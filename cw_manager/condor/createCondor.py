from . import writeCondor as wc
import numpy as np
from pathlib import Path
from ..utils import filePath as fp
from ..utils import setup_parameter as setup
from ..utils import utils as utils
from tqdm import tqdm

class followupManager:
    def __init__(self, target, obsDay):
        self.obsDay = obsDay
        self.setup = setup
        self.target = target
    
    def analyzeResultArgs(self, cohDay, freq, stage, freqDerivOrder,  numTopList, cluster inj, workInLocalDir): 
        argListString = "argList=\" --target {0} --obsDay {1} --cohDay {2} --freq {3} --stage {4} --freqDerivOrder {5} --numTopList {6} --sftFiles {7}".format(
                self.target.name, self.obsDay, cohDay, freq, stage, freqDerivOrder, self.numTopList, sftFiles)
        for arg in [cluster, inj, workInLocalDir]
            if cluster:
                argListString += " --{0}".format(arg)
        argListString += "\""
        return argListString

    def weaveArgs(self, freq, params, taskName, sftFiles, jobIndex, OSG=True):
        # using OSG computing resources (different format for .sub file)
        argList+="OUTPUTFILE=\"{0}\" REMAPOUTPUTFILE=\"{1}\" TRANSFERFILES=\"{2}\" ".format(
            Path(resultFile).name, resultFile, inputFiles)
        return argList
    
###### main use function 
    def makeAnalyzeDag(self, cohDay, freq, stage, freqDerivOrder=2,  numTopList=1000, cluster=False, inj=False, workInLocalDir=False, OSG=False, OSDF=False):
        
        for jobIndex, freq in tqdm(enumerate(range(fmin, fmax), 1)):
            taskName = taskName = utils.taskName(self.target, new_stage, new_cohDay, new_freqDerivOrder, freq)
            exe = fp.
            subFileName = fp.condorSubFilePath(self.target, freq, taskName, self.stage)
            Path(subFileName).unlink(missing_ok=True)
            
            crFiles = fp.condorRecordFilePath(taskName, self.target, taskName, self.stage)
            utils.makeDir(crFiles)
            wc.writeSearchSub(subFileName, exe, True, crFiles[0], crFiles[1], crFiles[2], '$(argList)', request_memory='100MB', request_disk='100MB', OSG=OSG, OSDF=OSDF)
            
            dagFileName = fp.dagFilePath(freq, self.target, taskName, self.stage)
            Path(dagFileName).unlink(missing_ok=True
            # call function to write .sub files for analyze result
            ######################## Argument string use to write to DAG  ########################
            argList = self.analyzeResultArgs(old_cohDay, new_cohDay, freq, old_stage, new_stage, old_freqDerivOrder, new_freqDerivOrder,  numTopList, cluster, inj, workInLocalDir)
            # Call function from WriteCondorFiles.py which will write DAG 
            wc.writeSearchDag(dagFileName, taskName, subFileName, jobIndex, argList)
        print('Finish writing {0} dag files for {1}-{2}Hz'.format(self.stage, fmin, fmax))

