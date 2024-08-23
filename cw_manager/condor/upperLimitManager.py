#################### condor manager for estimate upper limit
from ..analysis.resultManager import resultManager
from . import writeCondor as wc
import numpy as np
from pathlib import Path
from ..utils import filePath as fp
from ..utils import utils as utils
from tqdm import tqdm

class upperLimitManager(resultManager):
    def __init__(self, target, obsDay):
        super().__init__(target, obsDay)    
    
    def estimateUpperLimitArgsStr(self): 
        argStr = ["alpha", "delta", "freq", "freq-band", "loudest-2F", "sft-patt", "output-file"]
        
        argListString = ""
        for s in argStr:
            argListString += "--{0}=$({1}) ".format(s, s.replace('-', '').upper())
        return argListString
    
    def estimateUpperLimitArgs(self, freq, fband, mean2F, sftFiles, resultFile): 
        inputFiles = ', '.join([s for s in sftFiles])
        sft = ';'.join([Path(s).name for s in sftFiles])
        argListString = 'ALPHA="{0}" DELTA="{1}" FREQ="{2}" FREQBAND="{3}" LOUDEST2F="{4}" SFTPATT="{5}" OUTPUTFILE="{6}" REMAPOUTPUTFILE="{7}" TRANSFERFILES="{8}"'.format(
                self.target.alpha, self.target.delta, freq, fband, mean2F, sft, Path(resultFile).name, resultFile, inputFiles)
        return argListString
    
    def makeEstimateUpperLimitDag(self, fmin=20, fmax=475, fBand=0.1, stage='ULEstimation'):
        taskName = 'ULEstimation_{0}Days'.format(self.obsDay)
        
        nonSatBandsFile = fp.nonSaturatedBandFilePath(self.target, fmin, fmax, stage='search')
        
        if Path(nonSatBandsFile ).is_file():
            nonSatBands = utils.loadNonSaturatedBand(self.target, fmin, fmax) # load non saturated sub-bands file for all 1Hz bands
        else:
            nonSatBands = utils.getNonSaturatedBand(self.target, fmin, fmax) # get non saturated bands  
        
        threshFilePath = fp.threshFilePath(self.target, fmin, fmax, stage='search') # file storing mean2F threshold
        
        subFileName = fp.condorSubFilePath(self.target, '', taskName, stage)
        argsStr = self.estimateUpperLimitArgsStr()
        exe = fp.estimateUpperLimitExcutable()
        
        wc.writeSearchSub(subFileName, exe, False, '/dev/null', '/dev/null', '/dev/null', argsStr, request_memory='1GB')

        dagFileName = fp.dagFilePath('', self.target, taskName+'_{0}-{1}Hz'.format(fmin, fmax), stage)
        Path(dagFileName).unlink(missing_ok=True)

        for jobIndex, freq in tqdm(enumerate(nonSatBands, 1)):
            
            sftFiles = utils.sftEnsemble(freq, self.obsDay)
            resultFile = fp.estimateUpperLimitFilePath(self.target, freq, taskName, stage)
            Path(Path(resultFile).resolve().parent).mkdir(parents=True, exist_ok=True)   

            mean2F = self.readMean2F(freq, threshFilePath)
            argListString = self.estimateUpperLimitArgs(freq, fBand, mean2F, sftFiles, resultFile)
            wc.writeSearchDag(dagFileName, taskName, subFileName, jobIndex, argListString)  

        print('Finish writing dag files for estimating upper strain limit for {0}-{1}Hz'.format(fmin, fmax))
