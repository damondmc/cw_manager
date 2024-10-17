from . import readFile as rf
from . import tools as tools
from astropy.io import fits
from astropy.table import Table, vstack
import numpy as np
from ..utils import filePath as fp
from ..utils import setup_parameter as setup
from ..genParam import frequencyRange as fr
from tqdm import tqdm
from ..utils import utils as utils
from pathlib import Path
import warnings
    
class resultManager():
    def __init__(self, target, obsDay):
        self.obsDay = obsDay
        self.setup = setup
        self.target = target
        
    def readMean2F(self, freq, threshFilePath):
        _freq, _mean2F = np.loadtxt(threshFilePath).T
        
        if len(_freq[_freq==freq]) == 1:
            mean2F_th = _mean2F[_freq==freq][0]
        else:
            idx = (np.abs(_freq - freq)).argmin()
            mean2F_th = _mean2F[idx]
        return mean2F_th

    def calMean2F_threshold(self, cohDay, freq, nJobs):
        nTemp = self._readTemplateCount(cohDay, freq, nJobs)
        mean2F_th = tools.mean2F_threshold(sum(nTemp), self.nSeg)            
        return mean2F_th

    def readFollowUpMean2F(self, dataFilePath):
        data = fits.getdata(dataFilePath, 1)
        mean2F = data['mean2F threshold']
        return mean2F, data.size

    def makeOutlierTable(self, dataFilePath, mean2F_th, toplistLimit=1000, freqDerivOrder=2):    
        data = fits.getdata(dataFilePath, 1)[:toplistLimit]
        mask = data['mean2F'] > mean2F_th
        data = Table(data[mask])       
        data.add_column(mean2F_th*np.ones(len(data)), name='mean2F threshold')
        
        # get spacing         
        spacing = utils.getSpacing(dataFilePath, freqDerivOrder)
        _, name = utils.phaseParamName(freqDerivOrder)
        for i in range(len(name)):
            data.add_column(spacing[name[i]]*np.ones(len(data)), name=name[i]) 
        return data
    
    #def makeOutlierTable(self, dataFilePath, mean2F_th, toplistLimit=1000, freqDerivOrder=2):    
    #    data = fits.getdata(dataFilePath, 1)[:toplistLimit]
    #    mask = data['mean2F'] > mean2F_th
    #    data = Table(data[mask])       
    #    data.add_column(mean2F_th*np.ones(len(data)), name='mean2F threshold') 
    #    return data
       
    def makeInjectionTable(self, dataFilePath, sr, freqDerivOrder):      
        inj = fits.getdata(dataFilePath, 2)
        inj = Table(inj)   
        aplus, across = inj['aPlus'], inj['aCross']
        h0 = 0.5*(2.*aplus+2.*np.sqrt(aplus**2-across**2) )
        inj.add_column(h0*np.ones(len(inj)), name='h0')
        inj.rename_column('refTime_s', 'refTime')   
        fn, dfn = utils.phaseParamName(freqDerivOrder)
        #_d = fits.getheader(dataFilePath)
        #spacing = {key: _d['HIERARCH ' + key] for key in dfn}  
        mask = np.full(sr['freq'].shape, True)
        #mask *= ( (sr['freq']-setup.injOutlier_nSpacing*sr['df']) < inj['Freq'] )
        #mask *= ( (sr['freq']+setup.injOutlier_nSpacing*sr['df']) > inj['Freq'] )
        #for i in range(1, len(fn)): # just match up to 2nd order
        for i in range(1, 2):
            mask *= ( (sr[fn[i]]-setup.injOutlier_nSpacing*sr[dfn[i]]) < inj[fn[i]] )
            mask *= ( (sr[fn[i]]+setup.injOutlier_nSpacing*sr[dfn[i]]) > inj[fn[i]] )
        sr = Table(sr[mask])[:1] # only follow up the loudest one which covering the injection to save the cost 
        #_inj = inj.copy()
        #if len(sr) !=0:
        #    inj = vstack([inj for _ in range(len(sr))])
        return sr, inj
    
    # function to write result from weave output in each 1Hz band
    def _writeSearchResult(self, freq, mean2F_th, nJobs, numTopListLimit=1000, stage='search', freqDerivOrder=2, cluster=False, workInLocalDir=False):
        taskName = utils.taskName(self.target, stage, self.cohDay, freqDerivOrder, freq)
             
        outlierTableList = []
        info_data = np.recarray((nJobs,), dtype=[(key, '>f8') for key in ['freq', 'jobIndex', 'outliers', 'saturated']]) 

        for i, jobIndex in enumerate(tqdm(range(1, nJobs+1))):
            weaveFilePath = fp.weaveOutputFilePath(self.target, freq, taskName, jobIndex, stage)
            if workInLocalDir:
                weaveFilePath = Path(weaveFilePath).name
            _outlier = self.makeOutlierTable(weaveFilePath, mean2F_th, numTopListLimit, freqDerivOrder)  

            if len(_outlier) == numTopListLimit:
                info_data[i] = freq, jobIndex, 0, 1  
            else:
                info_data[i] = freq, jobIndex, len(_outlier), 0
                outlierTableList.append( _outlier )
        
        #compute non-saturated bands
        sat = info_data['saturated'].reshape(10, int(nJobs/10)).sum(axis=1)
        idx = np.where(sat == 0)[0] 
        nonSatBand = np.recarray((idx.size,), dtype=[(key, '>f8') for key in ['nonSatBand']])
        nonSatBand['nonSatBand'] = freq + idx * setup.fBand
           
        # append all tables in the file into one
        # Create a PrimaryHDU object
        primary_hdu = fits.PrimaryHDU()

        # Add 'mean2Fth' to the header with a value of 1.0
        primary_hdu.header['HIERARCH mean2F_th'] = mean2F_th
        # Add grid size to the header (i.e., df, df1dot, ...)
        spacing = utils.getSpacing(weaveFilePath, freqDerivOrder)
        _, name = utils.phaseParamName(freqDerivOrder)
        for i in range(len(name)):
            primary_hdu.header['HIERARCH '+ name[i]] = spacing[name[i]]

        outlier_hdu =  fits.BinTableHDU(data=vstack(outlierTableList))
        info_hdu =  fits.BinTableHDU(data=info_data) 
        nsb_hdu =  fits.BinTableHDU(data=nonSatBand)

        outlier_hdul = fits.HDUList()
        outlier_hdul.append(primary_hdu)
        outlier_hdul.append(outlier_hdu)
        outlier_hdul.append(info_hdu)
        outlier_hdul.append(nsb_hdu)

        outlierFilePath = fp.outlierFilePath(self.target, freq, taskName, stage, cluster=False)
        if workInLocalDir:
            outlierFilePath = Path(outlierFilePath).name
        utils.makeDir([outlierFilePath])
        outlier_hdul.writeto(outlierFilePath, overwrite=True)    
       
        if cluster:
            if outlier_hdu.data.size > 1:
                cluster_hdul = fits.HDUList()
                centers_idx, cluster_size = self.clustering(outlier_hdu.data, freqDerivOrder) ##### under development
                cluster_data = outlier_hdu.data[centers_idx]
                cluster_hdu = fits.BinTableHDU(data=cluster_data)
     
                info_data = np.recarray((cluster_size.size,), dtype=[(key, '>f8') for key in ['freq', 'clusterIndex', 'noOutliersWithin']]) 
                for i in range(cluster_size.size):
                    info_data[i] = freq, i, cluster_size[i]
                
                info_hdu =  fits.BinTableHDU(data=info_data)
                cluster_hdul.append(primary_hdu)
                cluster_hdul.append(cluster_hdu)
                cluster_hdul.append(info_hdu)
                cluster_hdul.append(nsb_hdu)
            else:
                cluster_hdul = outlier_hdul
            outlierClusterFilePath = fp.outlierFilePath(self.target, freq, taskName, stage, cluster=cluster)
            if workInLocalDir:
                outlierClusterFilePath = Path(outlierClusterFilePath).name
            cluster_hdul.writeto(outlierClusterFilePath, overwrite=True)
        return outlierFilePath 
    
    #work flow to write search result for each 1Hz band in (fmin, fmax) Hz
    def writeSearchResult(self, cohDay, freq, df1dot, df2dot, numTopList=1000, stage='search', freqDerivOrder=2, fmin=20, fmax=475, cluster=False, workInLocalDir=False):
        self.cohDay, self.cohTime, self.nSeg, self.obsTime, self.refTime = utils.getTimeSetup(self.target.name, self.obsDay, cohDay)    

        nf1dots = fr.getNf1dot(freq, self.setup.fBand, self.target.tau, df1dot=df1dot)
        nf2dots = fr.getNf2dot(freq, self.setup.fBand, self.target.tau, df2dot=df2dot)
        nJobs = int(nf1dots*nf2dots/self.setup.fBand)
        try:
            threshFilePath = fp.threshFilePath(self.target, fmin, fmax, 'search')
            mean2F_th = self.readMean2F(freq, threshFilePath)
        except:
            print('calculating mean2F threshold...')
            mean2F_th = self.calMean2F_threshold(cohDay, freq, nJobs)           
            
        print('mean2F threshold = ', mean2F_th)
        #mean2F_th = mean2F_th * np.ones(nJobs)
        outlierFilePath = self._writeSearchResult(freq, mean2F_th, nJobs, numTopList, stage, freqDerivOrder, cluster, workInLocalDir)
        print('Finish writing search result for {0} Hz'.format(freq))
        return outlierFilePath

    # function to write result from weave output in each 1Hz band
    def _writeInjectionResult(self, freq, mean2F_th, nJobs, numTopListLimit=1000, stage='search', freqDerivOrder=2, workInLocalDir=False, cluster=False):
        taskName = utils.taskName(self.target, stage, self.cohDay, freqDerivOrder, freq)
        outlierTableList = []
        injTableList = []
        info_data =np.recarray((nJobs,), dtype=[(key, '>f8') for key in ['freq', 'jobIndex', 'outliers']]) 
  
        weaveFilePath = fp.weaveOutputFilePath(self.target, freq, taskName, 1, stage)
        if workInLocalDir:
                weaveFilePath = Path(weaveFilePath).name
        #spacing = utils.getSpacing(weaveFilePath, freqDerivOrder)
        for i, jobIndex in enumerate(range(1, nJobs+1)):
            weaveFilePath = fp.weaveOutputFilePath(self.target, freq, taskName, jobIndex, stage)
            if workInLocalDir:
                weaveFilePath = Path(weaveFilePath).name
            _outlier = self.makeOutlierTable(weaveFilePath, mean2F_th, numTopListLimit, freqDerivOrder)  
            #_outlier, _inj = self.makeInjectionTable(weaveFilePath, _outlier, spacing, freqDerivOrder)
            _outlier, _inj = self.makeInjectionTable(weaveFilePath, _outlier, freqDerivOrder)
            
            if len(_outlier) == 0:
                outlierTableList.append( _outlier )
            else:
                outlierTableList.append( _outlier )
                injTableList.append( _inj )
            info_data[i] = freq, jobIndex, len(_outlier)  

        # append all tables in the file into one
        # Create a PrimaryHDU object
        primary_hdu = fits.PrimaryHDU()

        # Add 'mean2Fth' to the header with a value of 1.0
        #primary_hdu.header['HIERARCH mean2F_th'] = mean2F_th
        # Add grid size to the header (i.e., df, df1dot, ...)
        #spacing = utils.getSpacing(weaveFilePath, freqDerivOrder)
        #_, name = utils.phaseParamName(freqDerivOrder)
        #for i in range(len(name)):
        #    primary_hdu.header['HIERARCH '+ name[i]] = spacing[name[i]]

        outlier_hdu =  fits.BinTableHDU(data=vstack(outlierTableList))
        info_hdu =  fits.BinTableHDU(data=info_data) 
      
        outlier_hdul = fits.HDUList()
        outlier_hdul.append(primary_hdu)
        outlier_hdul.append(outlier_hdu)
        #outlier_hdul.append(info_hdu)


        if len(injTableList) != 0:
            inj_hdu =  fits.BinTableHDU(data=vstack(injTableList))
            outlier_hdul.append(inj_hdu)
        else:
            inj_hdu =  fits.BinTableHDU()
            outlier_hdul.append(inj_hdu)
            print('No outlier.')
        
        outlier_hdul.append(info_hdu)
        outlierFilePath = fp.outlierFilePath(self.target, freq, taskName, stage, cluster=False)
        if workInLocalDir:
            outlierFilePath = Path(outlierFilePath).name
 
        utils.makeDir([outlierFilePath])
        outlier_hdul.writeto(outlierFilePath, overwrite=True)    
        return outlierFilePath

    #work flow to write injection-search result for each 1Hz band in (fmin, fmax) Hz
    def writeInjectionResult(self, cohDay, fmin, fmax, nBands, nInj, numTopList=1000, stage='search', freqDerivOrder=2, cluster=False):
        self.cohDay, self.cohTime, self.nSeg, self.obsTime, self.refTime = utils.getTimeSetup(self.target.name, self.obsDay, cohDay)    
        freqList = utils.loadNonSaturatedBand(self.target, fmin, fmax, nBands)
        for freq in tqdm(freqList):
            try:
                threshFilePath = fp.threshFilePath(self.target, fmin, fmax, 'search')
                mean2F_th = self.readMean2F(freq, threshFilePath)
            except:
                taskName = utils.taskName(self.target, 'search', self.cohDay, freqDerivOrder, int(freq))
                outlierFilePath = fp.outlierFilePath(self.target, int(freq), taskName, 'search', cluster=cluster)
                data = fits.getdata(outlierFilePath)
                mean2F_th = data['mean2F'][0]
                
                print('No mean2F threshold file')
                           
            mean2F_th = mean2F_th * np.ones(nInj)
            self._writeInjectionResult(freq, mean2F_th, nInj, numTopList, stage, freqDerivOrder, cluster)
        print('Finish writing injection result for freq: {0}-{1} Hz'.format(fmin, fmax))
        return 0
    
    def writeInjectionResult1Hz(self, cohDay, freq, nInj, numTopList=1000, stage='search', freqDerivOrder=2, fmin=20, fmax=475, workInLocalDir=False, cluster=False):
        self.cohDay, self.cohTime, self.nSeg, self.obsTime, self.refTime = utils.getTimeSetup(self.target.name, self.obsDay, cohDay)    
        
        try:
            taskName = utils.taskName(self.target, 'search', self.cohDay, freqDerivOrder, int(freq))
            outlierFilePath = fp.outlierFilePath(self.target, int(freq), taskName, 'search', cluster=cluster)
            if workInLocalDir:
                outlierFilePath = Path(outlierFilePath).name
            mean2F_th = fits.getheader(outlierFilePath)['HIERARCH mean2F_th']
        except:
            print('No mean2F threshold file')
                
        outlierFilePath = self._writeInjectionResult(freq, mean2F_th, nInj, numTopList, stage, freqDerivOrder, workInLocalDir, cluster)
        print('Finish writing injection result for {0} Hz'.format(freq))
        return outlierFilePath

    # function to read template count from weave output in each 1Hz band
    def _readTemplateCount(self, cohDay, freq, nJobs, stage='search', freqDerivOrder=2):
        templateList = []
        taskName = utils.taskName(self.target, stage, cohDay, freqDerivOrder, freq)
        crfiles = fp.condorRecordFilePath(freq, self.target, taskName, stage)
        for jobIndex in range(1, nJobs+1):
            outFilePath = crfiles[0][:-8] + '{0}'.format(jobIndex)
            templateList.append(rf.readTemplateCount(outFilePath))
        return templateList
    
    # function to read runtime from weave output in each 1Hz band
    def _readJobStat(self, cohDay, freq, nJobs, stage='search', freqDerivOrder=2):
        timeList, memoryList = [], []
    
    # function to read runtime from weave output in each 1Hz band
    def _readJobStat(self, cohDay, freq, nJobs, stage='search', freqDerivOrder=2):
        timeList, memoryList = [], []
        taskName = utils.taskName(self.target, stage, cohDay, freqDerivOrder, freq)
        crfiles = fp.condorRecordFilePath(freq, self.target, taskName, stage)
        for jobIndex in range(1, nJobs+1):
            outFilePath = crfiles[0][:-8] + '{0}'.format(jobIndex)
            timeList.append(rf.readRunTime(outFilePath))
            memoryList.append(rf.readMemoryUsage(outFilePath))
        return timeList, memoryList
                
    def writeJobStat(self, cohDay, fmin, fmax, df1dot=None, df2dot=None, nJobs=None, stage='search', freqDerivOrder=2):          
        self.cohDay, self.cohTime, self.nSeg, self.obsTime, self.refTime = utils.getTimeSetup(self.target.name, self.obsDay, cohDay)    
        if stage == 'search':
            print('Working on mean2F threshold.')
            mean2FfilePath = fp.threshFilePath(self.target, fmin, fmax, stage)
            if not Path(mean2FfilePath).is_file():
                with open(mean2FfilePath, 'wt') as file:
                    file.write('#{0}\t{1}\n'.format('freq', 'mean2F threshold'))
                    for freq in tqdm(range(fmin, fmax)):
                        nf1dots = fr.getNf1dot(freq, self.setup.fBand, self.target.tau, df1dot)
                        nf2dots = fr.getNf2dot(freq, self.setup.fBand, self.target.tau, df2dot)
                        mean2F = self.calMean2F_threshold(cohDay, freq, nf1dots*nf2dots)
                        file.write('{0}\t{1}\n'.format(freq, mean2F))
                        
            freqList = range(fmin, fmax)
            filePath = fp.jobStatFilePath(self.target, fmin, fmax, stage)
            with open(filePath, 'wt') as file:
                file.write('#{0}\t{1}\t{2}\t{3}\t{4}\n'.format('freq', 'totalTimeUsage(s)', 'minMemoryUsage(MB)', 'maxMemoryUsage(MB)', 'no.Templates'))
                for freq in tqdm(freqList):
                    nf1dots = fr.getNf1dot(freq, self.setup.fBand, self.target.tau, df1dot)
                    nf2dots = fr.getNf2dot(freq, self.setup.fBand, self.target.tau, df1dot)
                    nJobs = int(round(nf1dots*nf2dots/self.setup.fBand, 1))
                    times, memorys = self._readJobStat(cohDay, freq, nJobs)
                    nTemp = self._readTemplateCount(cohDay, freq, nJobs)
                    file.write('{0}\t{1}\t{2}\t{3}\t{4}\n'.format(freq, sum(times), min(memorys), max(memorys), sum(nTemp)))                      
        elif 'followUp' in stage:
            print('Working on follow up.')
            freqList = [int(f) for f in nJobs.keys()]
            filePath = fp.jobStatFilePath(self.target, fmin, fmax, stage)
            with open(filePath, 'wt') as file:
                file.write('#{0}\t{1}\t{2}\t{3}\t{4}\n'.format('freq', 'totalTimeUsage(s)', 'minMemoryUsage(MB)', 'maxMemoryUsage(MB)', 'no.Templates'))
                for i, freq in enumerate(tqdm(freqList)):
                    times, memorys = self._readJobStat(cohDay, freq, nJobs[str(freq)], stage=stage, freqDerivOrder=freqDerivOrder)
                    nTemp = self._readTemplateCount(cohDay, freq, nJobs[str(freq)], stage=stage, freqDerivOrder=freqDerivOrder)
                    if nJobs[str(freq)]>0:
                        file.write('{0}\t{1}\t{2}\t{3}\t{4}\n'.format(freq, sum(times), min(memorys), max(memorys), sum(nTemp)))
                    else:
                        print('{0} Hz has 0 outlier to follow up.'.format(freq))
                        #file.write('{0}\t{1}\t{2}\t{3}\t{4}\n'.format(freq, 'nan', 'nan', 'nan', 'nan'))
                        
        elif 'injection' in stage:
            print('Working on injection test')
            #freqList = utils.loadNonSaturatedBand(self.target, fmin, fmax, nBands)
            freqList = range(fmin, fmax)
            filePath = fp.jobStatFilePath(self.target, fmin, fmax, stage)
            with open(filePath, 'wt') as file:
                file.write('#{0}\t{1}\t{2}\t{3}\t{4}\n'.format('freq', 'totalTimeUsage(s)', 'minMemoryUsage(MB)', 'maxMemoryUsage(MB)', 'no.Templates'))
                for freq in tqdm(freqList):
                    times, memorys = self._readJobStat(cohDay, freq, nJobs, stage=stage, freqDerivOrder=freqDerivOrder)
                    nTemp = self._readTemplateCount(cohDay, freq, nJobs, stage=stage, freqDerivOrder=freqDerivOrder)
                    file.write('{0}\t{1}\t{2}\t{3}\t{4}\n'.format(freq, sum(times), min(memorys), max(memorys), sum(nTemp)))

        return filePath

    def infoSummary(self, cohDay, fmin, fmax, df1dot=None, df2dot=None, nJobs=None, stage='search', freqDerivOrder=2, nBands=None):              
        filePath = fp.infoSummaryFilePath(self.target, fmin, fmax, stage)
        if 'search' in stage:
            print('Working on first stage of the search')
            freqList = range(fmin, fmax)
            with open(filePath, 'wt') as file:
                file.write('#{0}\t{1}\t{2}\n'.format('freq', 'outliers', 'saturated'))
                for freq in tqdm(freqList):
                    taskName = utils.taskName(self.target, stage, cohDay, freqDerivOrder, freq)
                    infoFilePath  = fp.outlierInfoFilePath(self.target, freq, taskName, stage, cluster=False)
                    _freq, jobIndex, noutlier, nsat = np.loadtxt(infoFilePath).T
                    nf1dots = fr.getNf1dot(freq, self.setup.fBand, self.target.tau, df1dot)
                    nf2dots = fr.getNf2dot(freq, self.setup.fBand, self.target.tau, df2dot)
                    
                    totalOutliers = noutlier.reshape(-1, nf1dots*nf2dots).sum(axis=-1)
                    total_saturated = nsat.reshape(-1, nf1dots*nf2dots).sum(axis=-1)
                    for i in range(len(totalOutliers)):
                        file.write('{0}\t{1}\t{2}\n'.format(freq+i*self.setup.fBand, totalOutliers[i], total_saturated[i]))
        
        elif 'followUp' in stage:
            print('Working on follow up.')
            freqList = freqList = np.array([int(f) for f in nJobs.keys()])
            with open(filePath, 'wt') as file:
                file.write('#{0}\t{1}\t{2}\n'.format('freq', 'outliers', 'saturated'))

                for freq in tqdm(freqList):
                    taskName = utils.taskName(self.target, stage, cohDay, freqDerivOrder, freq)
                    infoFilePath  = fp.outlierInfoFilePath(self.target, freq, taskName, stage, cluster=False)
                    if nJobs[str(freq)]>0:
                        _, _, noutlier, nsat = np.loadtxt(infoFilePath).T
                        file.write('{0}\t{1}\t{2}\n'.format(freq, noutlier.sum(), nsat.sum()))
                    else:
                        print('{0} Hz has 0 outlier to follow up.'.format(freq))
                        file.write('{0}\t{1}\t{2}\n'.format(freq, 'nan', 'nan'))
        elif 'injection' in stage:
            print('Working on injection test')
            freqList = utils.loadNonSaturatedBand(self.target, fmin, fmax, nBands)
            with open(filePath, 'wt') as file:
                file.write('#{0}\t{1}\t{2}\n'.format('freq', 'outliers', 'saturated'))

                for freq in tqdm(freqList):
                    taskName = utils.taskName(self.target, stage, cohDay, freqDerivOrder, freq)
                    infoFilePath  = fp.outlierInfoFilePath(self.target, freq, taskName, stage, cluster=False)
                    _, _, noutlier, nsat = np.loadtxt(infoFilePath).T
                    file.write('{0}\t{1}\t{2}\n'.format(freq, noutlier.sum(), nsat.sum()))
       
        return filePath

    
    def calInjectionEfficiency(self, cohDay, fmin, fmax, nBands=None, nInj=1000, nAmp=8, stage='search', freqDerivOrder=2, cluster=False):
        nInjPerAmp = int(nInj/nAmp)
        freqList = utils.loadNonSaturatedBand(self.target, fmin, fmax, nBands)
        saveFilePath = fp.efficiency_nonSaturatedBandFilePath(self.target, fmin, fmax, nBands=nBands)
        with open(saveFilePath, 'wt') as file:
            file.write('#{0}'.format('freq'))
            for i in range(nAmp):
                file.write('\t{0}'.format(i))
            file.write('\n')
            
            for freq in tqdm(freqList):
                file.write('{0}'.format(freq))

                taskName = utils.taskName(self.target, stage, cohDay, freqDerivOrder, freq)
                infoFilePath  = fp.outlierInfoFilePath(self.target, freq, taskName, 'injectionUpperLimit', cluster=False)
                _freq, _, nout, nsat = np.loadtxt(infoFilePath).T
                
                nout[(nout!=0)+(nsat!=0)] = 1 # if no. outlier > 0, we will follow up and find that signal
                nout = nout.reshape(nAmp, nInjPerAmp)
                nout = np.sum(nout, axis=1)
                p = nout/float(nInjPerAmp)
                
                if p.size != nAmp:
                    print('Error: size of pdet is not consistent with number of h0 points.')
                
                for i in range(nAmp):
                    file.write('\t{0}'.format(p[i]))
                file.write('\n') 
                
        return saveFilePath
 
    def calInjectionEfficiency1Hz(self, cohDay, fmin, fmax, nInj=1000, nAmp=8, stage='search', freqDerivOrder=2):
        nInjPerAmp = int(nInj/nAmp)
        saveFilePath = fp.efficiencyIn1HzFilePath(self.target, fmin, fmax, stage)
        with open(saveFilePath, 'wt') as file:
            file.write('#{0}'.format('freq'))
            for i in range(nAmp):
                file.write('\t{0}'.format(i))
            file.write('\n')
            
            for freq in tqdm(range(fmin,fmax)):
                file.write('{0}'.format(freq))
                taskName = utils.taskName(self.target, stage, cohDay, freqDerivOrder, freq)
                infoFilePath  = fp.outlierInfoFilePath(self.target, freq, taskName, stage, cluster=False)
                _freq, _, nout, nsat = np.loadtxt(infoFilePath).T
                
                nout[(nout!=0)+(nsat!=0)] = 1 # if no. outlier > 0, we will follow up and find that signal
                nout = nout.reshape(nAmp, nInjPerAmp)
                nout = np.sum(nout, axis=1)
                p = nout/float(nInjPerAmp)
                
                if p.size != nAmp:
                    print('Error: size of pdet is not consistent with number of h0 points.')
                
                for i in range(nAmp):
                    file.write('\t{0}'.format(p[i]))
                file.write('\n') 
        return saveFilePath
 
    def _writeFollowUpResult(self, cohDay, freq, mean2F_th, nJobs, numTopListLimit=1000, stage='search', freqDerivOrder=2, workInLocalDir=True, inj=False, cluster=False):
        taskName = utils.taskName(self.target, stage, cohDay, freqDerivOrder, freq)
    
        outlierTableList = []
        injTableList = []
        info_data =np.recarray((nJobs,), dtype=[(key, '>f8') for key in ['freq', 'jobIndex', 'outliers']]) 
 
        weaveFilePath = fp.weaveOutputFilePath(self.target, freq, taskName, 1, stage)
        if workInLocalDir:
            weaveFilePath = Path(weaveFilePath).name
            
        #spacing = utils.getSpacing(weaveFilePath, freqDerivOrder)
        for i, jobIndex in enumerate(range(1, nJobs+1)):
            weaveFilePath = fp.weaveOutputFilePath(self.target, freq, taskName, jobIndex, stage)
            if workInLocalDir:
                weaveFilePath = Path(weaveFilePath).name
            _outlier = self.makeOutlierTable(weaveFilePath, mean2F_th[i], numTopListLimit, freqDerivOrder)  
            if inj:
                #_outlier, _inj = self.makeInjectionTable(weaveFilePath, _outlier, spacing, freqDerivOrder)
                _outlier, _inj = self.makeInjectionTable(weaveFilePath, _outlier, freqDerivOrder)
                
            if len(_outlier) == 0:
                outlierTableList.append( _outlier )
            else:
                outlierTableList.append( _outlier )
                if inj:
                    injTableList.append( _inj )
                    
            info_data[i] = freq, jobIndex, len(_outlier)  

           
        # append all tables in the file into one
        # Create a PrimaryHDU object
        primary_hdu = fits.PrimaryHDU()

        # Add 'mean2Fth' to the header with a value of 1.0
        #primary_hdu.header['HIERARCH mean2F_th'] = mean2F_th
        # Add grid size to the header (i.e., df, df1dot, ...)
        #spacing = utils.getSpacing(weaveFilePath, freqDerivOrder)
        #_, name = utils.phaseParamName(freqDerivOrder)
        #for i in range(len(name)):
        #    primary_hdu.header['HIERARCH '+ name[i]] = spacing[name[i]]

        outlier_hdul = fits.HDUList()
        outlier_hdul.append(primary_hdu)
       
        if len(outlierTableList) !=0:
            outlier_hdu =  fits.BinTableHDU(data=vstack(outlierTableList))
        else:
            outlier_hdu =  fits.BinTableHDU()
            print('No outlier.')
        outlier_hdul.append(outlier_hdu)
        
        # if software injection is included 
        if inj:
            if len(injTableList) != 0:
                inj_hdu =  fits.BinTableHDU(data=vstack(injTableList))
                outlier_hdul.append(inj_hdu)
            else:
                inj_hdu =  fits.BinTableHDU()
        
        # summary information of the outliers
        info_hdu =  fits.BinTableHDU(data=info_data)
        outlier_hdul.append(info_hdu)

        if inj:
            outlierFilePath = fp.outlierFilePath(self.target, freq, taskName, stage, cluster=False)
        else:
            outlierFilePath = fp.outlierFilePath(self.target, freq, taskName, stage, cluster=cluster)
        if workInLocalDir:
            outlierFilePath = Path(outlierFilePath).name
        utils.makeDir([outlierFilePath])
        outlier_hdul.writeto(outlierFilePath, overwrite=True)    
        
        
        if not inj and cluster:
            if outlier_hdu.data.size > 1:
                cluster_hdul = fits.HDUList()
                centers_idx, cluster_size = self.clustering(outlier_hdu.data, freqDerivOrder) ##### under development
                cluster_data = outlier_hdu.data[centers_idx]
                cluster_hdu = fits.BinTableHDU(data=cluster_data)

                #if inj:
                #    inj_data = inj_hdu.data[centers_idx]
                #    inj_hdu = fits.BinTableHDU(data=inj_data) 

                info_data = np.recarray((cluster_size.size,), dtype=[(key, '>f8') for key in ['freq', 'clusterIndex', 'noOutliersWithin']]) 
                for i in range(cluster_size.size):
                    info_data[i] = freq, i, cluster_size[i]
                info_hdu =  fits.BinTableHDU(data=info_data)
                cluster_hdul.append(primary_hdu)
                cluster_hdul.append(cluster_hdu)
                #if inj:
                #    cluster_hdul.append(inj_hdu)
                #cluster_hdul.append(info_hdu)
            else:
                cluster_hdul = outlier_hdul
            outlierClusterFilePath = fp.outlierFilePath(self.target, freq, taskName, stage, cluster=cluster)
            if workInLocalDir:
                outlierClusterFilePath = Path(outlierClusterFilePath).name
            cluster_hdul.writeto(outlierClusterFilePath, overwrite=True)
        return outlierFilePath 
 
    #work flow to write injection-search result for each 1Hz band in (fmin, fmax) Hz
    def writeFollowUpResult(self, old_cohDay, new_cohDay, freq, numTopList=1000, old_stage='search', new_stage='followUp-1', old_freqDerivOrder=2, new_freqDerivOrder=2, ratio=None, workInLocalDir=True, inj=False, cluster=False):

        taskName = utils.taskName(self.target, old_stage, old_cohDay, old_freqDerivOrder, freq)
        outlierFilePath = fp.outlierFilePath(self.target, freq, taskName, old_stage, cluster=cluster)
        if workInLocalDir:
            outlierFilePath = Path(outlierFilePath).name
        
        data = fits.getdata(outlierFilePath, 1)
            
        if inj:
            mean2F_th = data['mean2F'] # no actual ratio increase threshold for inj test, keep all outlier for investigation
            nJobs = mean2F_th.size
            mean2F_th = np.zeros(nJobs)
        else:
            if ratio is None:
                ratio = utils.followUpMean2F_ratio(self.target, new_stage)
            mean2F_th = data['mean2F'] * ratio
            nJobs = mean2F_th.size
            print('ratio=',ratio)

        outlierFilePath = self._writeFollowUpResult(new_cohDay, freq, mean2F_th, nJobs, numTopList, new_stage, new_freqDerivOrder, workInLocalDir, inj, cluster)
        print('Finish writing followUp result for {0} Hz'.format(freq))
        return outlierFilePath

    def emsembleFollowUpResult(self, outlierFilePathList, inj_outlierFilePathList, mean2F_ratio_list,
                               freq, final_stage, taskName, workInLocalDir=False, cluster=False):
    
        stage = ['followUp-{}'.format(i+1) for i in range(len(outlierFilePathList))]

        primary_hdu = fits.PrimaryHDU()
        outlier_hdul = fits.HDUList()
        #inj_outlier_hdul = fits.HDUList()
          
        for i in range(len(outlierFilePathList)):
            primary_hdu.header['HIERARCH mean2F_ratio_{}'.format(stage[i])] = mean2F_ratio_list[i]
        outlier_hdul.append(primary_hdu)
        
        if len(inj_outlierFilePathList) !=0:
            data = fits.getdata(inj_outlierFilePathList[0], 1) 
            outlier_hdul.append(fits.BinTableHDU(data=data, name='inj_search'))
         
            for i in range(len(inj_outlierFilePathList)-1):           
                data = fits.getdata(inj_outlierFilePathList[i+1], 1)
                outlier_hdul.append(fits.BinTableHDU(data=data, name='inj_'+stage[i]))
            
        if len(outlierFilePathList) != 0:
            for i in range(len(outlierFilePathList)):
                data = fits.getdata(outlierFilePathList[i], 1)
                outlier_hdul.append(fits.BinTableHDU(data=data, name=stage[i]))
   
        outlierFilePath = fp.outlierFilePath(self.target, freq, taskName, final_stage, cluster=cluster)    
        if workInLocalDir:
            outlierFilePath = Path(outlierFilePath).name        
        else:
            utils.makeDir([outlierFilePath])    
        outlier_hdul.writeto(outlierFilePath, overwrite=True)
        return outlierFilePath


    def clustering(self, data, freqDerivOrder):
        # convert table data into normalized array
        fn, dfn = utils.phaseParamName(freqDerivOrder)
        _data = [data[key] for key in fn]
        _data = np.column_stack(_data)

        _spacing = [data[key] for key in dfn]
        _spacing = np.column_stack(_spacing)

        loudness = data['mean2F']

        # Sort the data by loudness in descending order
        sorted_indices = np.argsort(-loudness)
        sorted_coords = _data[sorted_indices]
        sorted_loudness = loudness[sorted_indices]
        sorted_spacing = _spacing[sorted_indices]
        # Initialize list for keeping track of loudest samples and their respective spheres
        centers_idx = []
        cluster_size = []
        # Set to keep track of processed indices
        processed_indices = set()
        # Loop over sorted samples (from loudest to least loud)
        for i, (center, gridsize) in enumerate(zip(sorted_coords, sorted_spacing)):
            if i in processed_indices:
                continue  # Skip already processed samples

            # Initialize list for dimension-wise indices
            within_dim_indices = []
            # Compute distances in each dimension separately and find indices within radius r0
            for dim in range(freqDerivOrder+1):
                r0 = setup.injOutlier_nSpacing * gridsize[dim]
                distances_dim = np.abs(sorted_coords[:, dim] - center[dim])
                within_dim = np.where(distances_dim <= r0)[0]
                within_dim_indices.append(within_dim)

            # Combine the conditions (intersection of all dimensions)
            within_r0_indices = within_dim_indices[0]
            for dim_indices in within_dim_indices[1:]:
                within_r0_indices = np.intersect1d(within_r0_indices, dim_indices)
            
            processed_indices.update(within_r0_indices)  # Mark all samples within this sphere as processed
            centers_idx.append(sorted_indices[i])
            cluster_size.append(len(within_r0_indices))
        # Convert lists to arrays
        centers_idx = np.array(centers_idx)
        cluster_size = np.array(cluster_size)
        return centers_idx, cluster_size


"""       
        if cluster:
            cluster_hdul = fits.HDUList()
            centers_idx, cluster_size = self.clustering(outlier_hdu.data, freqDerivOrder) ##### under development
            cluster_data = outlier_hdu.data[centers_idx]
            cluster_hdu = fits.BinTableHDU(data=cluster_data)

            inj_data = inj_hdu.data[centers_idx]
            inj_hdu = fits.BinTableHDU(data=inj_data) 

            info_data = np.recarray((cluster_size.size,), dtype=[(key, '>f8') for key in ['freq', 'clusterIndex', 'noOutliersWithin']]) 
            for i in range(cluster_size.size):
                info_data[i] = freq, i, cluster_size[i]
            info_hdu =  fits.BinTableHDU(data=info_data)
            cluster_hdul.append(primary_hdu)
            cluster_hdul.append(cluster_hdu)
            cluster_hdul.append(inj_hdu)
            cluster_hdul.append(info_hdu)
           
            outlierClusterFilePath = fp.outlierFilePath(self.target, freq, taskName, stage, cluster=cluster)
            if workInLocalDir:
                outlierClusterFilePath = Path(outlierClusterFilePath).name
            cluster_hdul.writeto(outlierClusterFilePath, overwrite=True)
"""   


