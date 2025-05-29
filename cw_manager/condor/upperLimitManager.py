from . import writeCondor as wc
import numpy as np
from pathlib import Path
from ..utils import filePath as fp
from ..utils import setup_parameter as setup
from ..utils import utils as utils
from tqdm import tqdm

class upperLimitManager:
    def __init__(self, target, obsDay):
        self.obsDay = obsDay
        self.setup = setup
        self.target = target
    
    def upperLimitArgs(self, cohDay, freq, stage, freqDerivOrder,  numTopList, nInj, skyUncertainty, num_cpus, sftFiles, est_sftFiles, cluster, workInLocalDir, OSDF): 
        sftFiles = ';'.join([Path(s).name for s in sftFiles])
        est_sftFiles = ';'.join([Path(s).name for s in est_sftFiles])
        argListString = '--target {0} --obsDay {1} --cohDay {2} --freq {3} --stage {4} --freqDerivOrder {5} --numTopList {6} --nInj {7} --sftFiles {8} --est_sftFiles {9} --num_cpus {10} --skyUncertainty {11}'.format(
                self.target.name, self.obsDay, cohDay, freq, stage, freqDerivOrder, numTopList, nInj, sftFiles, est_sftFiles, num_cpus, skyUncertainty)
        if cluster:
            argListString += " --cluster"
        
        if workInLocalDir:
            argListString += " --workInLocalDir"
        
        if OSDF:
            argListString += " --OSDF"
        return argListString

    def transferFileArgs(self, cohDay, freq, freqDerivOrder, stage, sftFiles, est_sftFiles, cluster=False, OSG=True):
    
        taskName = utils.taskName(self.target, 'search', cohDay, freqDerivOrder, freq)
        searchResultFile = fp.outlierFilePath(self.target, freq, taskName, 'search', cluster=cluster) 
        
        exe = fp.upperLimitExecutableFilePath()
        image = fp.imageFilePath()
        inputFiles = "{}, {}, {}".format(exe, image, searchResultFile)
        for sft in sftFiles:
            inputFiles += ", {}".format(sft)
        for sft in est_sftFiles:
            inputFiles += ", {}".format(sft)

        _, cohTime, nSeg, _, _ = utils.getTimeSetup(self.target.name, self.obsDay, cohDay)
        metric = fp.weaveSetupFilePath(cohTime, nSeg, freqDerivOrder)
        inputFiles += ", {}".format(metric)
        
        # using OSG computing resources (different format for .sub file)
        taskName = utils.taskName(self.target, stage, cohDay, freqDerivOrder, freq)
        outlierFilePath = fp.outlierFilePath(self.target, freq, taskName, stage, cluster=cluster) 
        utils.makeDir([outlierFilePath])
        argList="OUTPUTFILE=\"{0}\" REMAPOUTPUTFILE=\"{1}\" TRANSFERFILES=\"{2}\" ".format(
            Path(outlierFilePath).name, outlierFilePath, inputFiles)
        return argList
    
    def makeUpperLimitDag(self, fmin, fmax, cohDay, freqDerivOrder=2, stage='upperLimit', skyUncertainty=1e-4, nInj=100, numTopList=1000, num_cpus=4, request_memory='4GB', request_disk='4GB', cluster=False, workInLocalDir=False, OSG=False, OSDF=False):
              
        taskName = utils.taskName(self.target, stage, cohDay, freqDerivOrder, str(fmin)+'-'+str(fmax)) 
        dagFileName = fp.dagFilePath('', self.target, taskName, stage)
        Path(dagFileName).unlink(missing_ok=True)
        
        for jobIndex, freq in tqdm(enumerate(range(fmin, fmax), 1)):
            taskName = utils.taskName(self.target, stage, cohDay, freqDerivOrder, freq)
            exe = fp.upperLimitExecutableFilePath()
            exe = Path(exe).name
            subFileName = fp.condorSubFilePath(self.target, freq, taskName, stage)
            Path(subFileName).unlink(missing_ok=True)
            
            crFiles = fp.condorRecordFilePath(freq, self.target, taskName, stage)
            utils.makeDir(crFiles)
            
            image = fp.imageFilePath()
            image = Path(image).name
            sftFiles = utils.sftEnsemble(freq, self.obsDay, OSDF=OSDF)
            est_sftFiles = utils.sftEnsemble(freq, cohDay, OSDF=OSDF)
            argList = self.upperLimitArgs(cohDay, freq, stage, freqDerivOrder, numTopList, nInj, skyUncertainty, num_cpus, sftFiles, est_sftFiles, cluster, workInLocalDir, OSDF)
            wc.writeSearchSub(subFileName, exe, True, crFiles[0], crFiles[1], crFiles[2], argList, request_memory=request_memory, request_disk=request_disk, OSG=OSG, OSDF=OSDF, image=image)
            
           # call function to write .sub files for analyze result
            ######################## Argument string use to write to DAG  ########################        
            argList = self.transferFileArgs(cohDay, freq, freqDerivOrder, stage, sftFiles, est_sftFiles, cluster, OSG)
                             
            # Call function from WriteCondorFiles.py which will write DAG 
            wc.writeSearchDag(dagFileName, taskName, subFileName, jobIndex, argList)
        print('Finish writing upper limit dag from {0} stage for {1}-{2}Hz'.format(stage, fmin, fmax))

