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
    
    def makeInjectionTable(self, dataFilePath, sr, freqDerivOrder):      
        inj = fits.getdata(dataFilePath, 2)
        inj = Table(inj)  
        aplus, across = inj['aPlus'], inj['aCross']
        h0 = 0.5*(2.*aplus+2.*np.sqrt(aplus**2-across**2) )
        inj.add_column(h0*np.ones(len(inj)), name='h0')
        inj.rename_column('refTime_s', 'refTime')   
        fn, dfn = utils.phaseParamName(freqDerivOrder)
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
    def _writeSearchResult(self, freq, mean2F_th, nJob, numTopListLimit=1000, stage='search', freqDerivOrder=2, cluster=False):
        taskName = utils.taskName(self.target, stage, self.cohDay, freqDerivOrder, freq)
    
        infoFilePath  = fp.outlierInfoFilePath(self.target, freq, taskName, stage, cluster=False)
        utils.makeDir([infoFilePath])

        outlierTableList = []
        with open(infoFilePath, 'w') as infoFile:    
            infoFile.write('#{0}\t{1}\t{2}\t{3}\n'.format('freq', 'jobIndex', 'outliers', 'saturated'))
            
            for i, jobIndex in enumerate(tqdm(range(1, nJob+1))):
                weaveFilePath = fp.weaveOutputFilePath(self.target, freq, taskName, jobIndex, stage)
                _outlier = self.makeOutlierTable(weaveFilePath, mean2F_th[i], numTopListLimit, freqDerivOrder)  

                if len(_outlier) == numTopListLimit:
                    infoFile.write('{0}\t{1}\t{2}\t{3}\n'.format(freq, jobIndex, 0, 1))
                else:
                    infoFile.write('{0}\t{1}\t{2}\t{3}\n'.format(freq, jobIndex, len(_outlier), 0))
                    outlierTableList.append( _outlier )
                    
        # append all tables in the file into one
        outlier_hdu =  fits.BinTableHDU(data=vstack(outlierTableList))
        
        outlier_hdul = fits.HDUList()
        outlier_hdul.append(outlier_hdu)
        
        outlierFilePath = fp.outlierFilePath(self.target, freq, taskName, stage, cluster=False)
        utils.makeDir([outlierFilePath])
        outlier_hdul.writeto(outlierFilePath, overwrite=True)    
        
        #if cluster:
        #    cluster_hdul = tools.cluster(outlier_hdul) ##### under development
        #    outlierClusterFilePath = fp.outlierFilePath(self.target, freq, taskName, stage, cluster=cluster)
        #    cluster_hdul.writeto(outlierClusterFilePath, overwrite=True)
        return 0 
    
    #work flow to write search result for each 1Hz band in (fmin, fmax) Hz
    def writeSearchResult(self, cohDay, freq, df1dot, df2dot, numTopList=1000, stage='search', freqDerivOrder=2, fmin=20, fmax=475, cluster=False):
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
        mean2F_th = mean2F_th * np.ones(nJobs)
        self._writeSearchResult(freq, mean2F_th, nJobs, numTopList, stage, freqDerivOrder, cluster)
        print('Finish writing search result for {0} Hz'.format(freq))
        return 0

    # function to write result from weave output in each 1Hz band
    def _writeInjectionResult(self, freq, mean2F_th, nJob, numTopListLimit=1000, stage='search', freqDerivOrder=2, cluster=False):
        taskName = utils.taskName(self.target, stage, self.cohDay, freqDerivOrder, freq)
    
        infoFilePath  = fp.outlierInfoFilePath(self.target, freq, taskName, stage, cluster=False)
        utils.makeDir([infoFilePath])

        outlierTableList = []
        injTableList = []
        with open(infoFilePath, 'w') as infoFile:    
            infoFile.write('#{0}\t{1}\t{2}\t{3}\n'.format('freq', 'jobIndex', 'outliers', 'saturated'))
            
            for i, jobIndex in enumerate(range(1, nJob+1)):
                weaveFilePath = fp.weaveOutputFilePath(self.target, freq, taskName, jobIndex, stage)
                _outlier = self.makeOutlierTable(weaveFilePath, mean2F_th[i], numTopListLimit, freqDerivOrder)  
                _outlier, _inj = self.makeInjectionTable(weaveFilePath, _outlier, freqDerivOrder)
                if len(_outlier) == 0:
                    infoFile.write('{0}\t{1}\t{2}\t{3}\n'.format(freq, jobIndex, 0, 0))
                    outlierTableList.append( _outlier )
                else:
                    infoFile.write('{0}\t{1}\t{2}\t{3}\n'.format(freq, jobIndex, len(_outlier), 0))
                    outlierTableList.append( _outlier )
                    injTableList.append( _inj )

        # append all tables in the file into one
        outlier_hdu =  fits.BinTableHDU(data=vstack(outlierTableList))
         
        outlier_hdul = fits.HDUList()
        outlier_hdul.append(outlier_hdu)

        if len(injTableList) != 0:
            inj_hdu =  fits.BinTableHDU(data=vstack(injTableList))
            outlier_hdul.append(inj_hdu)
        else:
            print('No outlier.')
        
        outlierFilePath = fp.outlierFilePath(self.target, freq, taskName, stage, cluster=False)
        utils.makeDir([outlierFilePath])
        outlier_hdul.writeto(outlierFilePath, overwrite=True)    
        return 0 

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
    
    def writeInjectionResult1Hz(self, cohDay, freq, nInj, numTopList=1000, stage='search', freqDerivOrder=2, fmin=20, fmax=475, cluster=False):
        self.cohDay, self.cohTime, self.nSeg, self.obsTime, self.refTime = utils.getTimeSetup(self.target.name, self.obsDay, cohDay)    
        
        try:
            threshFilePath = fp.threshFilePath(self.target, fmin, fmax, 'search')
            mean2F_th = self.readMean2F(freq, threshFilePath)
        except:
            print('No mean2F threshold file')
               
        mean2F_th = mean2F_th * np.ones(nInj)
        self._writeInjectionResult(freq, mean2F_th, nInj, numTopList, stage, freqDerivOrder, cluster)
        print('Finish writing injection result for {0} Hz'.format(freq))


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
    
    # function to write result from weave output in each 1Hz band
    def _writeFollowUpResult(self, cohDay,  freq, mean2F_th, nJobs, numTopListLimit=1000, stage='search', freqDerivOrder=2, inj=False, cluster=False):
        taskName = utils.taskName(self.target, stage, cohDay, freqDerivOrder, freq)
    
        infoFilePath  = fp.outlierInfoFilePath(self.target, freq, taskName, stage, cluster=False)
        utils.makeDir([infoFilePath])

        outlierTableList = []
        injTableList = []
        with open(infoFilePath, 'w') as infoFile:    
            infoFile.write('#{0}\t{1}\t{2}\t{3}\n'.format('freq', 'jobIndex', 'outliers', 'saturated'))
            for i, jobIndex in enumerate(range(1, nJobs+1)):
                weaveFilePath = fp.weaveOutputFilePath(self.target, freq, taskName, jobIndex, stage)
                _outlier = self.makeOutlierTable(weaveFilePath, mean2F_th[i], numTopListLimit, freqDerivOrder) 
                #_outlier.add_column(mean2F_th*np.ones(len(data)), name='mean2F threshold')
                if inj:
                    _outlier, _inj = self.makeInjectionTable(weaveFilePath, _outlier, freqDerivOrder)
                if len(_outlier) == 0:
                    infoFile.write('{0}\t{1}\t{2}\t{3}\n'.format(freq, jobIndex, 0, 0))
                    outlierTableList.append( _outlier )
                else:
                    infoFile.write('{0}\t{1}\t{2}\t{3}\n'.format(freq, jobIndex, len(_outlier), 0))
                    outlierTableList.append( _outlier )
                    if inj:
                        injTableList.append( _inj )
        
        outlier_hdul = fits.HDUList()
        # append all tables in the file into one
        if len(outlierTableList) !=0:
            outlier_hdu =  fits.BinTableHDU(data=vstack(outlierTableList))
        else:
            outlier_hdu =  fits.BinTableHDU()
        outlier_hdul.append(outlier_hdu)

        if inj:
            if len(injTableList) != 0:
                inj_hdu =  fits.BinTableHDU(data=vstack(injTableList))
                outlier_hdul.append(inj_hdu)
            else:
                print('No outlier.')
                
        outlierFilePath = fp.outlierFilePath(self.target, freq, taskName, stage, cluster=False)
        utils.makeDir([outlierFilePath])
        outlier_hdul.writeto(outlierFilePath, overwrite=True)    
        return 0 

    #work flow to write injection-search result for each 1Hz band in (fmin, fmax) Hz
    def writeFollowUpResult(self, old_cohDay, new_cohDay, freq, numTopList=1000, old_stage='search', new_stage='followUp-1', old_freqDerivOrder=2, new_freqDerivOrder=2, inj=False, cluster=False):

        taskName = utils.taskName(self.target, old_stage, old_cohDay, old_freqDerivOrder, freq)
        outlierFilePath = fp.outlierFilePath(self.target, freq, taskName, old_stage, cluster=cluster)
        data = fits.getdata(outlierFilePath, 1)
            
        if inj:
            mean2F_th = data['mean2F'] # no actual ratio increase threshold for inj test, keep all outlier for investigation
        else:        
            ratio = utils.followUpMean2F_ratio(self.target, new_stage)
            mean2F_th = data['mean2F'] * ratio
            print('ratio=',ratio)
        nJobs = mean2F_th.size
        if inj:
            mean2F_th = np.zeros(nJobs)
        self._writeFollowUpResult(new_cohDay, freq, mean2F_th, nJobs, numTopList, new_stage, new_freqDerivOrder, inj, cluster)
        print('Finish writing followUp result for {0} Hz'.format(freq))
