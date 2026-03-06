from astropy.io import fits
from astropy.table import Table, vstack
import numpy as np
from tqdm import tqdm
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor

from paws.filepaths import PathManager
from paws.definitions import phase_param_name
from paws.io import make_dir, get_spacing 
from .clustering import clustering  

class ResultAnalysisManager:
    """
    Manages the collection, filtering, and storage of search results.
    """
    def __init__(self, config, target):
        """
        Initialize the ResultManager.

        Parameters:
            target (dict): Target object containing target information.
            config (dict): Configuration dictionary.
        """
        self.config = config
        self.target = target
        self.paths = PathManager(config, target)
        
    def make_outlier_table(self, data, mean2f_th, num_toplist=1000):    
        """Filters data to create an outlier table."""
        # Read and limit the data to the top entries
        data = data[:num_toplist]
        
        # Mask data with mean 2F values greater than the threshold
        mask = data['mean2F'] >= mean2f_th
        data = Table(data[mask])       
        data.add_column(mean2f_th * np.ones(len(data)), name='mean2F threshold')

        return data

    def make_injection_table(self, inj_param, search_param):   
        """Creates a table comparing injections with search results."""
        inj_param = Table(inj_param)   
        
        # Calculate h0 from aPlus and aCross
        aplus, across = inj_param['aPlus'], inj_param['aCross']
        h0 = aplus + np.sqrt(aplus**2 - across**2)
        inj_param.add_column(h0 * np.ones(len(inj_param)), name='h0')
        
        # Rename reference time if exists
        if 'refTime_s' in inj_param.colnames:
            inj_param.rename_column('refTime_s', 'refTime')   

        search_param = Table(search_param)[:1] 

        return search_param, inj_param

    def _collect_outlier_data(self, taskname, freq, stage, job_indices, mean2f_th, 
                              num_toplist, freq_deriv_order, work_in_local_dir, 
                              max_workers, read_inj=False, check_saturation=False, 
                              separate_saturated=False, desc="Processing"):
        """Central engine for multithreaded FITS reading and outlier filtering."""
        # 1. Handle scalar vs array thresholds
        if np.isscalar(mean2f_th):
            thresholds = [mean2f_th] * len(job_indices)
        else:
            thresholds = mean2f_th # Follow-up provides an array of thresholds
            
        outlier_table_list = []
        sat_outlier_table_list = []
        inj_table_list = []
        info_list = [] # Stores (freq, job_idx, n_outliers, is_saturated)
        max_spacing = {}
        
        # 2. Universal Worker Function
        def _worker(args):
            i, job_idx, th = args
            file_path = self.paths.weave_output_file(freq, taskname, job_idx, stage)
            if work_in_local_dir:
                file_path = Path(file_path).name
                
            try:
                weave_data = fits.getdata(file_path, 1)
                spacing = get_spacing(file_path, freq_deriv_order)

                _outlier = self.make_outlier_table(weave_data, th, num_toplist)
                is_sat = 1 if len(_outlier) >= num_toplist else 0
                
                # Handle injections
                _inj_param = None
                if read_inj and _outlier is not None:
                    inj_data = fits.getdata(file_path, extname='injection_info')
                    _outlier, _inj_param = self.make_injection_table(inj_data, _outlier)
                    
                return (i, job_idx, _outlier, _inj_param, spacing, is_sat)
            except FileNotFoundError:
                return (i, job_idx, None, None, None, 0)

        # 3. Multithreading Queue
        job_args = list(zip(range(len(job_indices)), job_indices, thresholds))
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            results = list(tqdm(executor.map(_worker, job_args), total=len(job_args), desc=f"{desc} {freq}Hz"))
            
        # 4. Safe Sequential Unpacking
        missing_files = 0
        for res in results:
            _, job_idx, _outlier, _inj_param, spacing, is_sat = res
            
            if _outlier is not None:
                info_list.append((freq, job_idx, len(_outlier), is_sat))
                
                if separate_saturated and is_sat == 1:
                    sat_outlier_table_list.append(_outlier[:1]) # Keep ONLY the loudest 1
                else:
                    outlier_table_list.append(_outlier)
                    
                if _inj_param is not None and len(_outlier) > 0:
                    inj_table_list.append(_inj_param)
                
                if spacing is not None:
                    if not max_spacing:
                        max_spacing = spacing.copy()
                    else:
                        for k, v in spacing.items():
                            max_spacing[k] = max(max_spacing.get(k, 0), v)
            else:
                if is_sat == 0: 
                    missing_files += 1
                info_list.append((freq, job_idx, 0, is_sat))
                
        if missing_files > 0:
            print(f"Warning: {missing_files} files missing for {desc} {freq}Hz")
            
        return outlier_table_list, sat_outlier_table_list, inj_table_list, info_list, max_spacing
    
    def _write_clustered_results(self, freq, taskname, stage, outlier_data, freq_deriv_order, 
                                 primary_hdu, inj_hdu=None, non_sat_hdu=None, work_in_local_dir=False):
        """Central engine for clustering outliers and writing the clustered FITS file."""
        
        cluster_hdul = fits.HDUList()
        
        if outlier_data is None or len(outlier_data) == 0:
            # Create an empty outlier table (preserves column schema if outlier_data is an empty recarray)
            if outlier_data is not None:
                cluster_hdu = fits.BinTableHDU(data=outlier_data, name=stage+'_outlier')
            else:
                cluster_hdu = fits.BinTableHDU(name=stage+'_outlier')
                
            # Create an empty info table
            dtypes = [(key, '>f8') for key in ['freq', 'clusterIndex', 'noOutliersWithin']]
            info_data = np.recarray((0,), dtype=dtypes) # 0 rows!
            info_clustered_hdu = fits.BinTableHDU(data=info_data, name='info_clustered')

        else:
            # 1. Run the clustering algorithm
            _, dfn = phase_param_name(freq_deriv_order)
            spacing = {key: primary_hdu.header[f'HIERARCH {key}'] for key in dfn}
            
            centers_idx, cluster_size, cluster_member = clustering(outlier_data, spacing, self.config['cluster_n_spacing']) 

            # 2. Map Data (Handle Injection vs Standard)
            if inj_hdu is not None:
                # Injection mode: map every outlier to a cluster center so we don't lose injection tracking
                center_idx_for_each_outlier = np.full(outlier_data.size, -1)
                processed_indices = set()
                for ci, members in zip(centers_idx, cluster_member):
                    idx = np.array([item for item in members if item not in processed_indices])
                    if len(idx) > 0:
                        center_idx_for_each_outlier[idx] = ci
                        processed_indices.update(members)
                cluster_data = outlier_data[center_idx_for_each_outlier]
            else:
                # Standard mode: just grab the centers
                cluster_data = outlier_data[centers_idx]

            cluster_hdu = fits.BinTableHDU(data=cluster_data, name=stage+'_outlier')

            # 3. Build Info Table 
            dtypes = [(key, '>f8') for key in ['freq', 'clusterIndex', 'noOutliersWithin']]
            info_data = np.recarray((len(cluster_size),), dtype=dtypes)
            for i in range(len(cluster_size)):
                info_data[i] = freq, i, cluster_size[i]
            
            info_clustered_hdu = fits.BinTableHDU(data=info_data, name='info_clustered')

        # 4. Assemble Final HDUList
        cluster_hdul.append(primary_hdu)
        cluster_hdul.append(cluster_hdu)
        if inj_hdu is not None:
            cluster_hdul.append(inj_hdu)
        if non_sat_hdu is not None:                 
            cluster_hdul.append(non_sat_hdu)
        cluster_hdul.append(info_clustered_hdu)

        # 5. File Path Logic
        outlier_file_path = self.paths.outlier_file(freq, taskname, stage, cluster=True)
            
        if work_in_local_dir:
            outlier_file_path = Path(outlier_file_path).name
            
        cluster_hdul.writeto(outlier_file_path, overwrite=True)
        return outlier_file_path

    def _make_search_outlier(self, taskname, freq, mean2f_th, n_jobs, num_toplist=1000, 
                             stage='search', freq_deriv_order=2, cluster=False, 
                             work_in_local_dir=False, separate_saturated=True, max_workers=32):
        
        job_indices = list(range(1, n_jobs + 1))
        
        normal_outliers, sat_outliers, _, info_list, max_spacing = self._collect_outlier_data(
            taskname, freq, stage, job_indices, mean2f_th, num_toplist, 
            freq_deriv_order, work_in_local_dir, max_workers, 
            read_inj=False, separate_saturated=True, desc="Search"
        )

        # 1. Build Info Array (NOW WITH 4 COLUMNS)
        info_data = np.recarray((n_jobs,), dtype=[(key, '>f8') for key in ['freq', 'jobIndex', 'outliers', 'isSaturated']])
        for i, (f, j, o, s) in enumerate(info_list):
            info_data[i] = (f, j, o, s)
            
        # 2. Build Non-Saturated Band Array
        f0_band = self.config['f0_band']
        sat = info_data['isSaturated'].reshape(int(1. / f0_band), int(n_jobs * f0_band))  

        idx = np.where(~sat.any(axis=1))[0]
        
        non_sat_data = np.recarray((len(idx),), dtype=[('non_sat_band', '>f8')])
        
        non_sat_data['non_sat_band'] = int(freq) + idx * f0_band
            
        # 3. Create HDUs
        primary_hdu = fits.PrimaryHDU()
        primary_hdu.header['HIERARCH mean2F_th'] = mean2f_th    

        if max_spacing:
            for key, val in max_spacing.items():
                primary_hdu.header[f'HIERARCH {key}'] = val    
            
        if normal_outliers:
            out_hdu = fits.BinTableHDU(data=vstack(normal_outliers), name=stage+'_outlier')
        else:
            out_hdu = fits.BinTableHDU(name=stage+'_outlier')

        if sat_outliers:
            sat_hdu = fits.BinTableHDU(data=vstack(sat_outliers), name=stage+'_sat_outlier')
        else:
            sat_hdu = fits.BinTableHDU(name=stage+'_sat_outlier')
            
        info_hdu = fits.BinTableHDU(data=info_data, name='info') 
        non_sat_hdu = fits.BinTableHDU(data=non_sat_data, name='non_sat_band')

        # 4. Assemble FITS with 5 HDUs
        outlier_hdul = fits.HDUList([primary_hdu, out_hdu, sat_hdu, info_hdu, non_sat_hdu])
        outlier_file_path = self.paths.outlier_file(freq, taskname, stage, cluster=False)
        if work_in_local_dir:
            outlier_file_path = Path(outlier_file_path).name
            
        make_dir([outlier_file_path])
        outlier_hdul.writeto(outlier_file_path, overwrite=True)  
        
        # 5. Clustering via central engine
        if cluster and out_hdu.data is not None:
            primary_hdu.header['HIERARCH cluster_n_spacing'] = self.config.get('cluster_n_spacing', 1)
            self._write_clustered_results(
                freq, taskname, stage, out_hdu.data, freq_deriv_order, 
                primary_hdu, inj_hdu=None, non_sat_hdu=non_sat_hdu,
                work_in_local_dir=work_in_local_dir
            )
            
        return outlier_file_path

    def make_search_outlier(self, taskname, freq, mean2f_th, n_jobs, num_toplist=1000, 
                            stage='search', freq_deriv_order=2, cluster=False, 
                            work_in_local_dir=False, separate_saturated=True, max_workers=16):
        """
        Public wrapper to write search results.

        Parameters:
        - taskname: str
            The name of the task for the search, used in naming and organizing output files.

        - freq: int
            The frequency value for the 1Hz band being processed.
            
        - mean2F_th: float
            The threshold value of the mean 2F statistic, which determines whether an outlier qualifies for follow-up or further analysis.

        - numTopList: int, optional (default=1000)
            Maximum number of top outliers to keep for each job's results.

        - stage: str, optional (default='search')
            The stage of the analysis. Determines the naming and organizational conventions for output files.

        - freqDerivOrder: int, optional (default=2)
            Specifies the order of frequency derivatives to consider (e.g., df1dot, df2dot) when calculating threshold and creating results.

        - cluster: bool, optional (default=False)
            If True, clusters outliers to consolidate similar results, saving computational costs and storage.

        - workInLocalDir: bool, optional (default=False)
            If True, stores output files in the local directory. This option might be useful for local testing.

        - separateSaturated: bool, optional (default=True)
            If True, separates saturated outliers from non-saturated ones.

        - max_worker: int 
            The numer of threads being used in the analysis.
        """ 
        
        outlier_file_path = self._make_search_outlier(taskname, freq, mean2f_th, n_jobs, num_toplist, 
                                                      stage, freq_deriv_order, cluster, work_in_local_dir, 
                                                      separate_saturated, max_workers)
        print('Finish writing search result for {0} Hz'.format(freq))
        return outlier_file_path 

    def _make_injection_outlier(self, taskname, freq, mean2f_th, n_jobs, num_toplist=1000, 
                                stage='injection', freq_deriv_order=2, cluster=False, 
                                work_in_local_dir=False, separate_saturated=False, max_workers=32):
        
        job_indices = list(range(1, n_jobs + 1))
        
        # read_inj=True: Read the injection files!
        outliers, _, inj_data_list, info_list, max_spacing = self._collect_outlier_data(
            taskname, freq, stage, job_indices, mean2f_th, num_toplist, 
            freq_deriv_order, work_in_local_dir, max_workers, 
            read_inj=True, separate_saturated=separate_saturated, desc="Injection"
        )
        
        info_data = np.recarray((n_jobs,), dtype=[(key, '>f8') for key in ['freq', 'jobIndex', 'outliers', 'isSaturated']])
        for i, (f, j, o, s) in enumerate(info_list):
            info_data[i] = (f, j, o, s)
            
        primary_hdu = fits.PrimaryHDU()
        primary_hdu.header['HIERARCH mean2F_th'] = mean2f_th        

        if max_spacing:
            for key, val in max_spacing.items():
                primary_hdu.header[f'HIERARCH {key}'] = val
            
        if outliers:
            out_hdu = fits.BinTableHDU(data=vstack(outliers), name=stage+'_outlier')
        else:
            out_hdu = fits.BinTableHDU(name=stage+'_outlier')
            
        info_hdu = fits.BinTableHDU(data=info_data, name='info') 
        
        # Build Injection HDU
        if inj_data_list:
            inj_hdu = fits.BinTableHDU(data=vstack(inj_data_list), name='injection_info')
        else:
            inj_hdu = fits.BinTableHDU(name='injection_info')

        # Assemble FITS with 4 HDUs (including injection_info)
        outlier_hdul = fits.HDUList([primary_hdu, out_hdu, inj_hdu, info_hdu])
        
        outlier_file_path = self.paths.outlier_file(freq, taskname, stage, cluster=False)
        if work_in_local_dir:
            outlier_file_path = Path(outlier_file_path).name
            
        make_dir([outlier_file_path])
        outlier_hdul.writeto(outlier_file_path, overwrite=True)  
        
        if cluster and out_hdu.data is not None:
            primary_hdu.header['HIERARCH cluster_n_spacing'] = self.config.get('cluster_n_spacing', 1)
            self._write_clustered_results(
                freq, taskname, stage, out_hdu.data, freq_deriv_order, 
                primary_hdu, inj_hdu=inj_hdu,
                work_in_local_dir=work_in_local_dir
            )
            
        return outlier_file_path

    def make_injection_outlier(self, taskname, freq, mean2f_th, n_jobs, num_toplist=1000, 
                               stage='search', freq_deriv_order=2, cluster=False, work_in_local_dir=False, 
                               separate_saturated=False, max_workers=32):
        """
        Public wrapper to write injection-search results.
        """

        outlier_file_path = self._make_injection_outlier(taskname, freq, mean2f_th, n_jobs, num_toplist, 
                                                         stage, freq_deriv_order, cluster, work_in_local_dir, 
                                                         separate_saturated=separate_saturated, max_workers=max_workers)
        print('Finish writing injection result for {0} Hz'.format(freq))
        return outlier_file_path

    def _make_followup_outlier(self, taskname, freq, mean2f_th, n_jobs, num_toplist=1000, 
                               stage='followup', freq_deriv_order=2, cluster=False, 
                               work_in_local_dir=False, separate_saturated=False, max_workers=32):
        
        job_indices = list(range(1, n_jobs + 1))
        
        # check_saturation=False: Keep all results for followup
        outliers, _, _, info_list, max_spacing = self._collect_outlier_data(
            taskname, freq, stage, job_indices, mean2f_th, num_toplist, 
            freq_deriv_order, work_in_local_dir, max_workers, 
            read_inj=False, separate_saturated=separate_saturated, desc="Followup"
        )
        
        info_data = np.recarray((n_jobs,), dtype=[(key, '>f8') for key in ['freq', 'jobIndex', 'outliers', 'isSaturated']])
        for i, (f, j, o, s) in enumerate(info_list):
            info_data[i] = (f, j, o, s)
            
        primary_hdu = fits.PrimaryHDU()

        if max_spacing:
            for key, val in max_spacing.items():
                primary_hdu.header[f'HIERARCH {key}'] = val
            
        if outliers:
            out_hdu = fits.BinTableHDU(data=vstack(outliers), name=stage+'_outlier')
        else:
            out_hdu = fits.BinTableHDU(name=stage+'_outlier')
            
        info_hdu = fits.BinTableHDU(data=info_data, name='info') 

        # Assemble FITS with 3 HDUs
        outlier_hdul = fits.HDUList([primary_hdu, out_hdu, info_hdu])
        
        outlier_file_path = self.paths.outlier_file(freq, taskname, stage, cluster=False)
        if work_in_local_dir:
            outlier_file_path = Path(outlier_file_path).name
            
        make_dir([outlier_file_path])
        outlier_hdul.writeto(outlier_file_path, overwrite=True)  
        
        if cluster and out_hdu.data is not None:
            primary_hdu.header['HIERARCH cluster_n_spacing'] = self.config.get('cluster_n_spacing', 1)
            self._write_clustered_results(
                freq, taskname, stage, out_hdu.data, freq_deriv_order, 
                primary_hdu, inj_hdu=None, 
                work_in_local_dir=work_in_local_dir
            )
            
        return outlier_file_path
    
    def make_followup_outlier(self, taskname, freq, mean2f_th, n_jobs, num_toplist=1000, 
                              new_stage='followUp-1', new_freq_deriv_order=2, 
                              cluster=False, work_in_local_dir=True, 
                              separate_saturated=False, max_workers=32):
        """
        Public wrapper to write follow-up results.
        """
        
        outlier_file_path = self._make_followup_outlier(taskname, freq, mean2f_th, n_jobs, num_toplist, 
                                                        new_stage, new_freq_deriv_order, cluster, work_in_local_dir, 
                                                        separate_saturated=separate_saturated, max_workers=max_workers)

        print(f'Finish writing followUp result for {freq} Hz')
        return outlier_file_path

    def ensemble_followup_result(self, freq, taskname, stage, inj_stage, outlier_file_path_list, inj_outlier_file_path_list, 
                                 mean2f_ratio_list, num_toplist_list,
                                 final_stage, cluster=False, work_in_local_dir=False):
        """Combines results from multiple follow-up stages into one summary FITS file."""
        n_inj_table = len(inj_outlier_file_path_list)
        n_out_table = len(outlier_file_path_list)
        
        primary_hdu = fits.PrimaryHDU()
        outlier_hdul = fits.HDUList()
        
        # Metadata
        try:
            # Try to get threshold from the first available file
            source_file = outlier_file_path_list[0] if n_out_table > 0 else (inj_outlier_file_path_list[0] if n_inj_table > 0 else None)
            if source_file:
                mean2f_th = fits.getheader(source_file)['HIERARCH mean2F_th']
                primary_hdu.header['HIERARCH mean2F_th'] = mean2f_th
        except (IndexError, KeyError, FileNotFoundError):
            print("Warning: Unable to retrieve mean2F_th from header.")
            pass

        primary_hdu.header['HIERARCH injection_test'] = (n_inj_table != 0)
            
        # Record ratios and top lists for every stage in the header
        # We iterate up to the max number of stages provided
        max_stages = max(n_inj_table, n_out_table)
        
        # Note: The loop index 'i' corresponds to the follow-up stage index. 
        # Typically stage lists include the initial search, so we might offset by 1 if 'stage' list includes 'search' at index 0.
        for i in range(max_stages):
            # Check bounds for ratio list
            if i < len(mean2f_ratio_list):
                # We use stage[i+1] assuming the lists passed in include the initial search stage name at 0
                stage_name = stage[i+1] if (i+1) < len(stage) else f"stage_{i+1}"
                primary_hdu.header[f'HIERARCH mean2F_ratio_{stage_name}'] = mean2f_ratio_list[i]
            
            # Check bounds for top list
            if i < len(num_toplist_list):
                stage_name = stage[i+1] if (i+1) < len(stage) else f"stage_{i+1}"
                primary_hdu.header[f'HIERARCH numTopList_{stage_name}'] = num_toplist_list[i]
                
        outlier_hdul.append(primary_hdu)
                 
        # 1. Append Injection Follow-up Stages
        for i in range(n_inj_table):           
            try:
                # Outliers
                data = fits.getdata(inj_outlier_file_path_list[i], extname=inj_stage[i]+'_outlier')
                outlier_hdul.append(fits.BinTableHDU(data=data, name=inj_stage[i]+'_outlier'))
                
                # Injections
                data = fits.getdata(inj_outlier_file_path_list[i], extname='inj') 
                outlier_hdul.append(fits.BinTableHDU(data=data, name=inj_stage[i]+'_inj'))

                # Info
                data = fits.getdata(inj_outlier_file_path_list[i], extname='info')
                outlier_hdul.append(fits.BinTableHDU(data=data, name=inj_stage[i]+'_info'))
            except FileNotFoundError:
                print(f"Warning: Missing injection file {inj_outlier_file_path_list[i]}")

        # 2. Append Search Follow-up Stages
        for i in range(n_out_table):
            try:
                data = fits.getdata(outlier_file_path_list[i], extname=stage[i]+'_outlier')
                outlier_hdul.append(fits.BinTableHDU(data=data, name=stage[i]+'_outlier'))

                data = fits.getdata(outlier_file_path_list[i], extname='info')
                outlier_hdul.append(fits.BinTableHDU(data=data, name=stage[i]+'_info'))
            except FileNotFoundError:
                print(f"Warning: Missing outlier file {outlier_file_path_list[i]}")
        
        # Write Final Ensemble File
        outlier_file_path = self.paths.outlier_file(freq, taskname, final_stage, cluster=cluster)    
        
        if work_in_local_dir:
            outlier_file_path = Path(outlier_file_path).name        
        else:
            make_dir([outlier_file_path])    
        
        outlier_hdul.writeto(outlier_file_path, overwrite=True)
        return outlier_file_path
