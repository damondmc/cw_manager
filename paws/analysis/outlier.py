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
        
    def make_outlier_table(self, data, spacing, mean2f_th, num_toplist=1000):    
        """Filters data to create an outlier table."""
        # Read and limit the data to the top entries
        data = data[:num_toplist]
        
        # Mask data with mean 2F values greater than the threshold
        mask = data['mean2F'] > mean2f_th
        data = Table(data[mask])       
        data.add_column(mean2f_th * np.ones(len(data)), name='mean2F threshold')
    
        # Get parameter names (e.g., ['f0', 'f1dot', ...])
        _, deriv_params = phase_param_name(len(spacing) - 1)
        
        for param in deriv_params:
            data.add_column(spacing[param] * np.ones(len(data)), name=param) 
        return data

    def make_injection_table(self, inj_param, search_param):   
        """Creates a table comparing injections with search results."""
        inj_param = Table(inj_param)   
        
        # Calculate h0 from aPlus and aCross
        aplus, across = inj_param['aPlus'], inj_param['aCross']
        h0 = 0.5 * (2. * aplus + 2. * np.sqrt(aplus**2 - across**2))
        inj_param.add_column(h0 * np.ones(len(inj_param)), name='h0')
        
        # Rename reference time if exists
        if 'refTime_s' in inj_param.colnames:
            inj_param.rename_column('refTime_s', 'refTime')   

        search_param = Table(search_param)[:1] 

        return search_param, inj_param

    def _collect_outlier_data(self, taskname, freq, stage, job_indices, mean2f_th, 
                              num_toplist, freq_deriv_order, work_in_local_dir, 
                              max_workers, read_inj=False, check_saturation=False, desc="Processing"):
        """Central engine for multithreaded FITS reading and outlier filtering."""
        # 1. Handle scalar vs array thresholds
        if np.isscalar(mean2f_th):
            thresholds = [mean2f_th] * len(job_indices)
        else:
            thresholds = mean2f_th # Follow-up provides an array of thresholds
            
        outlier_table_list = []
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
                _outlier = self.make_outlier_table(weave_data, spacing, th, num_toplist)
                
                # Handle saturation (If searching and too many outliers, discard the data)
                is_sat = 0
                if check_saturation and len(_outlier) >= num_toplist:
                    is_sat = 1
                    _outlier = None 
                
                # Handle injections
                _inj_param = None
                if read_inj and _outlier is not None:
                    inj_data = fits.getdata(file_path, 2)
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
            i, job_idx, _outlier, _inj_param, spacing, is_sat = res
            
            if _outlier is not None:
                outlier_table_list.append(_outlier)
                if _inj_param is not None and len(_outlier) > 0:
                    inj_table_list.append(_inj_param)
                    
                info_list.append((freq, job_idx, len(_outlier), is_sat))
                
                # Track max spacing
                if spacing is not None:
                    if not max_spacing:
                        max_spacing = spacing.copy()
                    else:
                        for k, v in spacing.items():
                            max_spacing[k] = max(max_spacing.get(k, 0), v)
            else:
                if is_sat == 0: 
                    missing_files += 1
                # Record 0 outliers if missing, or num_toplist if saturated
                n_out = num_toplist if is_sat else 0
                info_list.append((freq, job_idx, n_out, is_sat))
                
        if missing_files > 0:
            print(f"Warning: {missing_files} files missing for {desc} {freq}Hz")
            
        return outlier_table_list, inj_table_list, info_list, max_spacing
    
    def _write_clustered_results(self, freq, taskname, stage, outlier_data, freq_deriv_order, 
                                 primary_hdu, inj_hdu=None, chunk_index=None, work_in_local_dir=False):
        """Central engine for clustering outliers and writing the clustered FITS file."""
        if outlier_data is None or outlier_data.size <= 1:
            return None # Nothing to cluster

        cluster_hdul = fits.HDUList()
        
        # 1. Run the clustering algorithm
        centers_idx, cluster_size, cluster_member = clustering(outlier_data, freq_deriv_order) 

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

        # 3. Build Info Table (Handle chunking columns dynamically)
        if chunk_index is not None:
            dtypes = [(key, '>f8') for key in ['freq', 'chunkIndex', 'clusterIndex', 'noOutliersWithin']]
            info_data = np.recarray((cluster_size.size,), dtype=dtypes)
            for i in range(cluster_size.size):
                info_data[i] = freq, chunk_index, i, cluster_size[i]
        else:
            dtypes = [(key, '>f8') for key in ['freq', 'clusterIndex', 'noOutliersWithin']]
            info_data = np.recarray((cluster_size.size,), dtype=dtypes)
            for i in range(cluster_size.size):
                info_data[i] = freq, i, cluster_size[i]
        
        info_clustered_hdu = fits.BinTableHDU(data=info_data, name='info_clustered')

        # 4. Assemble Final HDUList
        cluster_hdul.append(primary_hdu)
        cluster_hdul.append(cluster_hdu)
        if inj_hdu is not None:
            cluster_hdul.append(inj_hdu)
        cluster_hdul.append(info_clustered_hdu)

        # 5. File Path Logic
        outlier_file_path = self.paths.outlier_file(freq, taskname, stage, cluster=True)
        if chunk_index is not None:
            outlier_file_path = outlier_file_path.replace('.fts', f'_chunk{chunk_index}.fts')
            
        if work_in_local_dir:
            outlier_file_path = Path(outlier_file_path).name
            
        cluster_hdul.writeto(outlier_file_path, overwrite=True)
        return outlier_file_path

    def _make_search_outlier(self, taskname, freq, mean2f_th, n_jobs, num_toplist=1000, 
                             stage='search', freq_deriv_order=2, cluster=False, 
                             work_in_local_dir=False, max_workers=32):
        
        job_indices = list(range(1, n_jobs + 1))
        
        # FIX 1: Expect 4 return values from engine (added trailing `_`)
        outliers, _, info_list, _ = self._collect_outlier_data(
            taskname, freq, stage, job_indices, mean2f_th, num_toplist, 
            freq_deriv_order, work_in_local_dir, max_workers, 
            read_inj=False, check_saturation=True, desc="Search"
        )
        
        info_data = np.recarray((n_jobs,), dtype=[(key, '>f8') for key in ['freq', 'jobIndex', 'outliers']])
        
        # FIX 2: Unpack 4 items from info_list (added `s` for is_sat)
        for i, (f, j, o, s) in enumerate(info_list):
            info_data[i] = (f, j, o)
            
        primary_hdu = fits.PrimaryHDU()
        primary_hdu.header['HIERARCH mean2F_th'] = mean2f_th        
            
        if outliers:
            out_hdu = fits.BinTableHDU(data=vstack(outliers), name=stage+'_outlier')
        else:
            out_hdu = fits.BinTableHDU(name=stage+'_outlier')
            
        info_hdu = fits.BinTableHDU(data=info_data, name='info') 

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

    def make_search_outlier(self, taskname, freq, mean2f_th, n_jobs, num_toplist=1000, 
                            stage='search', freq_deriv_order=2, cluster=False, work_in_local_dir=False, max_workers=16):
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

        - max_worker: int 
            The numer of threads being used in the analysis.
        """ 
        
        outlier_file_path = self._make_search_outlier(taskname, freq, mean2f_th, n_jobs, num_toplist, 
                                                      stage, freq_deriv_order, cluster, work_in_local_dir, max_workers)
        print('Finish writing search result for {0} Hz'.format(freq))
        return outlier_file_path 
    
    def _make_search_outlier_from_saturated_band(self, taskname, freq, mean2f_th, job_index, 
                                                 num_toplist=1, stage='search', freq_deriv_order=2, 
                                                 work_in_local_dir=False, max_workers=32):
        
        # FIX: Expect 4 return values
        outliers, _, _, _ = self._collect_outlier_data(
            taskname, freq, stage, job_index, mean2f_th, num_toplist, 
            freq_deriv_order, work_in_local_dir, max_workers, 
            read_inj=False, desc="Processing Sat Bands"
        )
           
        primary_hdu = fits.PrimaryHDU()
        primary_hdu.header['HIERARCH mean2F_th'] = mean2f_th

        if len(outliers) == 0:
            outlier_hdu = fits.BinTableHDU(name=stage+'SatBand_outlier')
        else:
            outlier_hdu = fits.BinTableHDU(data=vstack(outliers), name=stage+'SatBand_outlier')
        
        outlier_hdul = fits.HDUList([primary_hdu, outlier_hdu])
        
        taskname = taskname + '_satband'
        outlier_file_path = self.paths.outlier_file(freq, taskname, stage, cluster=False)
        
        if work_in_local_dir:
            outlier_file_path = Path(outlier_file_path).name
            
        make_dir([outlier_file_path])
        outlier_hdul.writeto(outlier_file_path, overwrite=True)  
       
        return outlier_file_path
   
    def make_search_outlier_from_saturated_band(self, taskname, freq, mean2f_th, job_index, 
                                                num_toplist=1, stage='search', freq_deriv_order=2, 
                                                work_in_local_dir=False, max_workers=32):
        """
        Public wrapper for saturated band results.

        Parameters:
        - taskname: str
            The name of the task, used for naming and organizing output files.

        - freq: int
            The frequency value for the 1Hz band being processed.
        
        - mean2F_th: float
            The threshold value of the mean 2F statistic, which determines whether an outlier qualifies for follow-up or further analysis.
            
        - jobIndex: int array
            Array of indices identifying which jobs within the 1Hz band are saturated, used for tracking and organizing job-specific results.

        - numTopList: int, optional (default=1)
            Maximum number of top outliers to keep for each job's results.

        - stage: str, optional (default='search')
            The stage of the analysis. Determines the naming and organizational conventions for output files.

        - freqDerivOrder: int, optional (default=2)
            Specifies the order of frequency derivatives to consider (e.g., df1dot, df2dot) when calculating threshold and creating results.

        - workInLocalDir: bool, optional (default=False)
            If True, stores output files in the local directory. This option might be useful for local testing.

        - max_worker: int
            The number of threads being used in the analysis.
        """ 

        outlier_file_path = self._make_search_outlier_from_saturated_band(taskname, freq, mean2f_th, job_index, 
                                                                          num_toplist, stage, freq_deriv_order, 
                                                                          work_in_local_dir, max_workers)
        print('Finish writing search result for {0} Hz'.format(freq))
        return outlier_file_path

    def _make_injection_outlier(self, taskname, freq, mean2f_th, n_jobs, num_toplist=1000, 
                                stage='search', freq_deriv_order=2, cluster=False, 
                                work_in_local_dir=False, max_workers=32):
                                
        job_indices = list(range(1, n_jobs + 1))
        
        # FIX 1: Expect 4 return values
        outliers, injs, info_list, _ = self._collect_outlier_data(
            taskname, freq, stage, job_indices, mean2f_th, num_toplist, 
            freq_deriv_order, work_in_local_dir, max_workers, 
            read_inj=True, desc="Inj Collection"
        )

        info_data = np.recarray((n_jobs,), dtype=[(key, '>f8') for key in ['freq', 'jobIndex', 'outliers']])
        
        # FIX 2: Unpack 4 items
        for i, (f, j, o, s) in enumerate(info_list):
            info_data[i] = (f, j, o)

        primary_hdu = fits.PrimaryHDU()
        primary_hdu.header['HIERARCH mean2F_th'] = mean2f_th
        primary_hdu.header['HIERARCH cluster_n_spacing'] = ''
        
        if outliers:
            outlier_hdu = fits.BinTableHDU(data=vstack(outliers), name=stage+'_outlier')
        else:
            outlier_hdu = fits.BinTableHDU(Table(), name=stage+'_outlier')
            
        info_hdu = fits.BinTableHDU(data=info_data, name='info') 
      
        if injs:
            inj_hdu = fits.BinTableHDU(data=vstack(injs), name='inj')
        else:
            inj_hdu = fits.BinTableHDU(name='inj')
            print('No outliers found overlapping with injections.')
        
        outlier_hdul = fits.HDUList([primary_hdu, outlier_hdu, inj_hdu, info_hdu])
        
        outlier_file_path = self.paths.outlier_file(freq, taskname, stage, cluster=False)
        if work_in_local_dir:
            outlier_file_path = Path(outlier_file_path).name
            
        make_dir([outlier_file_path])
        outlier_hdul.writeto(outlier_file_path, overwrite=True) 
                
        if cluster and outlier_hdu.data is not None:
            primary_hdu.header['HIERARCH cluster_n_spacing'] = self.config.get('cluster_n_spacing', 1)
            self._write_clustered_results(
                freq, taskname, stage, outlier_hdu.data, freq_deriv_order, 
                primary_hdu, inj_hdu=inj_hdu if injs else None, 
                work_in_local_dir=work_in_local_dir
            )
            
        return outlier_file_path

    def make_injection_outlier(self, taskname, freq, mean2f_th, n_jobs, num_toplist=1000, 
                               stage='search', freq_deriv_order=2, cluster=False, work_in_local_dir=False, max_workers=32):
        """
        Public wrapper to write injection-search results.

        Parameters:
        - taskname: str
            The name of the task for the search, used in naming and organizing output files.
        - freq: float
            The frequency in Hz for which results are being written.
        - mean2f_th: numpy.ndarray
            The mean 2F threshold value for the injections.
        - n_jobs: int
            The number of jobs processed.
        - num_top_list: int, optional (default=1000)
            Maximum number of top outliers to keep for each job's results.
        - stage: str, optional (default='search')
            The stage of the analysis.
        - freq_deriv_order: int, optional (default=2)
            The order of frequency derivative.
        - cluster: bool, optional (default=False)
            If True, clusters outliers to consolidate similar results.
        - work_in_local_dir: bool, optional (default=False)
            If True, stores output files in the local directory.
        - max_workers: int, optional (default=32)
            The number of threads to use for parallel processing when reading and filtering FITS files.
        """

        outlier_file_path = self._make_injection_outlier(taskname, freq, mean2f_th, n_jobs, num_toplist, 
                                                         stage, freq_deriv_order, cluster, work_in_local_dir, 
                                                         max_workers=max_workers)
        print('Finish writing injection result for {0} Hz'.format(freq))
        return outlier_file_path

    def _make_followup_outlier(self, taskname, freq, mean2f_th, n_jobs, num_toplist=1000, 
                               stage='search', freq_deriv_order=2, 
                               cluster=False, work_in_local_dir=True, inj=False,
                               chunk_index=0, chunk_size=1, max_workers=32):
                               
        start_job = chunk_index * chunk_size + 1
        end_job = chunk_index * chunk_size + n_jobs + 1
        job_indices = list(range(start_job, end_job))
        
        # FIX 1: Expect 4 return values
        outliers, injs, info_list, _ = self._collect_outlier_data(
            taskname, freq, stage, job_indices, mean2f_th, num_toplist, 
            freq_deriv_order, work_in_local_dir, max_workers, 
            read_inj=inj, desc=f"Follow-up Chunk {chunk_index}"
        )

        info_data = np.recarray((n_jobs,), dtype=[(key, '>f8') for key in ['freq', 'jobIndex', 'outliers']])
        
        # FIX 2: Unpack 4 items
        for i, (f, j, o, s) in enumerate(info_list):
            info_data[i] = (f, j, o)

        primary_hdu = fits.PrimaryHDU()
        
        if outliers:
            outlier_hdu = fits.BinTableHDU(data=vstack(outliers), name=stage+'_outlier')
        else:
            outlier_hdu = fits.BinTableHDU(name=stage+'_outlier')
            print('No outlier in follow-up chunk.')
            
        info_hdu = fits.BinTableHDU(data=info_data, name='info')
        outlier_hdul = fits.HDUList([primary_hdu, outlier_hdu])
        
        if inj and injs:
            inj_hdu = fits.BinTableHDU(data=vstack(injs), name='inj')
            outlier_hdul.append(inj_hdu)
        elif inj:
            inj_hdu = fits.BinTableHDU(name='inj')
            outlier_hdul.append(inj_hdu)
        else:
            inj_hdu = None
            
        outlier_hdul.append(info_hdu)
        
        outlier_file_path = self.paths.outlier_file(freq, taskname, stage, cluster=(cluster and not inj))
        if chunk_size != 1:
            outlier_file_path = outlier_file_path.replace('.fts', f'_chunk{chunk_index}.fts')
        if work_in_local_dir:
            outlier_file_path = Path(outlier_file_path).name
            
        make_dir([outlier_file_path])
        outlier_hdul.writeto(outlier_file_path, overwrite=True)     
        
        if cluster and outlier_hdu.data is not None:
            active_chunk_index = chunk_index if chunk_size != 1 else None
            self._write_clustered_results(
                freq, taskname, stage, outlier_hdu.data, freq_deriv_order, 
                primary_hdu, inj_hdu=inj_hdu, chunk_index=active_chunk_index,
                work_in_local_dir=work_in_local_dir
            )
            
        return outlier_file_path
    def make_followup_outlier(self, taskname, freq, mean2f_th, num_toplist=1000, 
                              new_stage='followUp-1', new_freq_deriv_order=2, 
                              cluster=False, work_in_local_dir=True, inj=False,
                              chunk_index=0, chunk_size=1, chunk_count=None, max_workers=32):
        """
        Public wrapper to write follow-up results.

        Parameters:
        - new_cohDay: int
            The new coherence day for the analysis.

        - freq: float
            The frequency in Hz for which results are being written.
            
        - old_mean2F: numpy.ndarray
            The mean 2F value of the outliers at previous stage (shorter coherence time).

        - numTopList: int, optional
            The maximum number of top results to consider. Default is 1000.

        - new_stage: str, optional
            The stage of the current analysis (e.g., 'followUp-1'). Default is 'followUp-1'.

        - new_freqDerivOrder: int, optional
            The order of frequency derivative used in the new analysis. Default is 2.

        - ratio: float, optional
            The ratio to adjust the mean2F threshold for the follow-up analysis. Default is None.

        - workInLocalDir: bool, optional
            If True, work with local directory paths. Default is True.

        - inj: bool, optional
            If True, includes injections in the follow-up result. Default is False
        
        - chunk_index: int, optional
            The index of the current chunk being processed (0-based). Default is 0.
        
        - chunk_size: int, optional
            The number of jobs included in each chunk. Default is 1 (no chunking).
        
        - chunk_count: int, optional
            The total number of chunks for the frequency band. If provided, it will be used to

        - max_workers: int, optional
            The number of threads to use for parallel processing when reading and filtering FITS files. Default is 32.
        """
        #print(f'Follow-up F-statistic threshold: {mean2f_th}')
        
        if chunk_count is not None:
            # Slice the threshold array for this specific chunk
            mean2f_th = mean2f_th[chunk_index*chunk_size : (chunk_index+1)*chunk_size]
            
        n_jobs = mean2f_th.size
        
        outlier_file_path = self._make_followup_outlier(taskname, freq, mean2f_th, n_jobs, num_toplist, 
                                                        new_stage, new_freq_deriv_order, 
                                                        cluster, work_in_local_dir, inj, 
                                                        chunk_index=chunk_index, chunk_size=chunk_size, max_workers=max_workers)

        print(f'Finish writing followUp result for {freq} Hz')
        return outlier_file_path


    def ensemble_outlier_chunk(self, chunk_count, taskname, freq, stage, cluster, work_in_local_dir):
        """
        Combines outlier results from multiple chunks into a single output file.

        Parameters:
        - chunk_count: int
            The number of chunks to process.
        - taskname: str
            The name of the task.
        - freq: float
            The frequency in Hz.
        - stage: str
            The current stage of the analysis.
        - cluster: bool
            If True, indicates that clustering results should be included.
        - work_in_local_dir: bool
            If True, indicates that paths should be treated as local directory paths.

        Returns:
        - outlier_file_path: str
            The path to the output file containing the combined outlier results.
        """
        outlier_file_path = self.paths.outlier_file(freq, taskname, stage, cluster=cluster)
        
        if work_in_local_dir:
            outlier_file_path = Path(outlier_file_path).name
            
        outlier_table_list = []
        info_table_list = []
        
        # Iterate through each chunk
        for i in range(chunk_count):
            # Construct chunk filename manually based on convention
            _outlier_file_path = outlier_file_path.replace('.fts', f'_chunk{i}.fts')
            
            try:
                # Read tables 
                outlier_table_list.append(Table(fits.getdata(_outlier_file_path, extname=stage+'_outlier')))  
                info_table_list.append(fits.getdata(_outlier_file_path, extname='info'))
            except FileNotFoundError:
                print(f"Warning: Chunk file {_outlier_file_path} missing.")

        # Stack and Write
        outlier_hdul = fits.HDUList()
        # Grab header from first available chunk if possible
        if chunk_count > 0:
             # Re-construct first chunk path
            first_chunk = outlier_file_path.replace('.fts', '_chunk0.fts')
            try:
                primary_hdu = fits.PrimaryHDU(header=fits.getheader(first_chunk))
            except FileNotFoundError:
                primary_hdu = fits.PrimaryHDU()
        else:
            primary_hdu = fits.PrimaryHDU()
            
        outlier_hdul.append(primary_hdu)
        
        if outlier_table_list:
            outlier_hdu = fits.BinTableHDU(data=vstack(outlier_table_list), name=stage+'_outlier')
            outlier_hdul.append(outlier_hdu)
        
        if info_table_list:
            info_hdu = fits.BinTableHDU(data=np.hstack(info_table_list), name='info')
            outlier_hdul.append(info_hdu) 
            
        outlier_hdul.writeto(outlier_file_path, overwrite=True)
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
