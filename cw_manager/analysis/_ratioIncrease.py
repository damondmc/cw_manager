import analysis.ReadFile as rf
import analysis.tools as tools
from astropy.io import fits
import numpy as np
import setup.setup_parameter as setup
import utils.utils as utils
from tqdm import tqdm
from condor.condorManager import condorManager


class resultManager(resultManager):    
    def __init__(self, target, newCohDay, newFreqDerivOrder, oldCohDay, oldFreqDerivOrder, obsDay, stage='Search', numTopList=1000):
        super().__init__(target, obsDay, oldCohDay, stage, oldFreqDerivOrder, numTopList)
        self.freqDerivOrder, self.newFreqDerivOrder = oldFreqDerivOrder, newFreqDerivOrder
        self.newCohDay, self.newCohTime, self.newNSeg, self.TObs, self.newRefTime = utils.getTimeSetup(target.name, obsDay, newCohDay)  
        
    def getNonSaturatedBand(self, fmin, fmax, nBand, saveFreqList=True):
        FilePath = self.setup.homeDir + 'Results/{0}/{0}_{1}_infoSummary_{2}-{3}Hz.txt'.format('Search', self.target.name, fmin, fmax)
        freq, _, nsat = np.loadtxt(FilePath).T
        freqList = freq[nsat==0]
        if len(freqList) > nBand:
            np.random.shuffle(freqList)
            nonSatBands = freqList[:nBand]
        else:
            nonSatBands = np.array([freqList[np.random.randint(0, len(freqList), 1)] for _ in range(nBand)])
            
        saveFilePath = self.setup.homeDir + 'Results/{0}/{0}_{1}_{2}-{3}Hz_{4}Bands.txt'.format('Search', self.target.name, fmin, fmax, nBand)
        if saveFreqList:
            with open(saveFilePath ,'w') as file:
                file.write('#freq(Hz)\n')
                for nonSatFreq in nonSatBands:
                    file.write('{0}\n'.format(nonSatFreq))
        return nonSatBands 
    
    def loadNonSaturatedBand(self, fmin, fmax, nBand):
        saveFilePath = self.setup.homeDir + 'Results/{0}/{0}_{1}_{2}-{3}Hz_{4}Bands.txt'.format('Search', self.target.name, fmin, fmax, nBand)
        nonSatBands = np.loadtxt(saveFilePath).T
        return nonSatBands
        
    def outlierOutput(self, freq, taskName):
        loudDataFile = setup.homeDir + 'Results/{0}/{1}/{2}/{3}/Outliers/{4}_outlier.fts'.format(self.stage, self.target.name, setup.sftSource, freq, taskName)
        infoFile = setup.homeDir + 'Results/{0}/{1}/{2}/{3}/Outliers/{4}_info.txt'.format(self.stage, self.target.name, setup.sftSource, freq, taskName)
        return loudDataFile, infoFile

    def injectionOutput(self, freq, taskName):
        injDataFile = setup.homeDir + 'Results/{0}/{1}/{2}/{3}/Outliers/{4}_injections.fts'.format(self.stage, self.target.name, setup.sftSource, freq, taskName)
        return injDataFile
    
    def readOutlierData(self, freq, cluster = False):
        taskName = self.taskName(freq)
        if cluster:
            dataFilePath, infoFilePath  = self.clusterOutput(freq, taskName)
        else:
            dataFilePath, infoFilePath  = self.outlierOutput(freq, taskName)
            
        hdul = fits.open(dataFilePath)
        data = hdul[1].data
        hdul.close()
        return data
 
        # function to write result from weave output in each 1Hz band
    def writeInjection(self, freq, mean2F_th, nJob, numTopListLimit=1000):
        taskName = self.taskName(freq)
        outlierFilePath, infoFilePath  = self.outlierOutput(freq, taskName)
        injectionFilePath = self.injectionOutput(freq, taskName)
        self.makeDir([infoFilePath])

        with open(infoFilePath, 'w') as infoFile:    
            infoFile.write('#{0}\t{1}\t{2}\t{3}\n'.format('freq', 'jobIndex', 'outliers', 'saturated'))
            
            injection_hdul = fits.HDUList()
            outlier_hdul = fits.HDUList()

            for jobIndex in range(1, nJob+1):
                weaveFilePath, _, _ = self.weaveOutput(freq, taskName, jobIndex)
                
                injection = tools.getInjection(weaveFilePath)   
                injection_hdul.append(injection)
                
                outlier = tools.filterOutlier(weaveFilePath, mean2F_th[jobIndex-1], numTopListLimit)   
                if len(outlier.data) == numTopListLimit:
                    empty_table = fits.BinTableHDU(outlier.data[:0])
                    outlier_hdul.append(empty_table)
                    infoFile.write('{0}\t{1}\t{2}\t{3}\n'.format(freq, jobIndex, 0, 1))
                else:
                    outlier_hdul.append(outlier)
                    infoFile.write('{0}\t{1}\t{2}\t{3}\n'.format(freq, jobIndex, len(outlier.data), 0))
        
        injection_hdul.writeto(injectionFilePath, overwrite=True)
        outlier_hdul.writeto(outlierFilePath, overwrite=True)
        return 0

    
    def writeInjectionResult(self, fmin, fmax, nBand, nJob, numTopListLimit=1000):
        freqList = self.loadNonSaturatedBand(fmin, fmax, nBand)
        for freq in tqdm(freqList):
            try:
                threshFilePath = self.threshFilePath(fmin, fmax)
                mean2F_th = self.readMean2F(freq, threshFilePath)
            except:
                print("Mean2F threshold file doesn't exist, calculate the mean2F threshold first...")
                threshFilePath = self.calMean2F_threshold(fmin, fmax)
                mean2F_th = self.readMean2F(freq, threshFilePath)
            mean2F_th = mean2F_th * np.ones(nJob)
            self.writeInjection(freq, mean2F_th, nJob, numTopListLimit)
        print('Finish writing result for freq: {0}-{1} Hz'.format(fmin, fmax))
        

    def writeFollowUpResult(self, params_list, numTopListLimit=1000, cluster_subband=0.01, clustering=True):
        freq_list = [int(f) for f in params_list.keys()]
        fmin, fmax = freq_list[0], freq_list[-1] 
        for freq in tqdm(freq_list):
            data=self.readOutlierData(freq, clustering=False)
            mean2F_th = data['mean2F'] * tools.mean2FIncRatio(self.stage)
            self.writeResult(freq, mean2F_th, nJob, numTopListLimit, cluster_subband, clustering)
        print('Finish writing result for freq: {0}-{1} Hz'.format(fmin, fmax))
