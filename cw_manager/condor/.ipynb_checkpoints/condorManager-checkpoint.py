import condor.writeCondor as wc
import numpy as np
from pathlib import Path
import filePath.filePath as fp
import setup.setup_parameter as setup
import utils.utils as utils
from tqdm import tqdm

class condorManager:
    def __init__(self, target, obsDay):
        self.obsDay = obsDay
        self.setup = setup
        self.target = target
    
    def weaveArgs(self, freq, params, taskName, sftFiles, jobIndex, OSG=True):
        exe = fp.weaveExecutableFilePath()
        metric = fp.weaveSetupFilePath(self.cohTime, self.nSeg, self.freqDerivOrder)
        
        resultFile = fp.weaveOutputFilePath(self.target, freq, taskName, jobIndex, self.stage)
        utils.makeDir([resultFile])
        
        extraStats = "coh2F_det,mean2F,coh2F_det,mean2F_det"
        if self.nSeg != 1:
            kwargs = {"semi-max-mismatch": self.setup.semiMM,
                      "coh-max-mismatch": self.setup.cohMM,
                      "toplist-limit": self.numTopList,
                      "extra-statistics": extraStats}
        else:
            kwargs = {"semi-max-mismatch": self.setup.semiMM,
                      "toplist-limit": self.numTopList,
                      "extra-statistics": extraStats}
        
        argList = ("")
        if not OSG:
            argList+="argList= \"--output-file={0} ".format(resultFile)
            argList+="--sft-files={0} ".format(';'.join([s for s in sftFiles]))
            argList+="--setup-file={0} ".format(metric)
            
            for key, value in kwargs.items():  
                argList+="--{0}={1} ".format(key,value)
            argList+="--alpha={0}/{1} ".format(params['alpha'], params['dalpha'])
            argList+="--delta={0}/{1} ".format(params['delta'], params['ddelta'])
        
            for i in range(self.freqDerivOrder+1):
                key1, key2 = self.freqParamName[i], self.freqDerivParamName[i]
                argList+="--{0}={1}/{2} ".format(key1, params[key1], params[key2])
            argList+="\""
        else: # using OSG computing resources (different format for .sub file)
            argList+="OUTPUTFILE=\"{0}\" ".format(Path(resultFile).name)
            argList+="REMAPOUTPUTFILE=\"{0}\" ".format(resultFile)
            argList+="SETUPFILE=\"{0}\" ".format(Path(metric).name)
            sft = ';'.join([Path(s).name for s in sftFiles])
            argList+="SFTFILES=\"{0}\" ".format(sft)
            inputFiles = ', '.join([s for s in sftFiles]) + ', ' + metric 
            argList+="TRANSFERFILES=\"{0}\" ".format(inputFiles)
            
            for key, value in kwargs.items():
                argList+="{0}=\"{1}\" ".format(key.replace('-', '').upper(),value)        
            argList+="ALPHA=\"{0}\" DALPHA=\"{1}\" ".format(params['alpha'], params['dalpha'])
            argList+="DELTA=\"{0}\" DDELTA=\"{1}\" ".format(params['delta'], params['ddelta'])
            for i in range(self.freqDerivOrder+1):
                key1, key2 = self.freqParamName[i], self.freqDerivParamName[i]
                argList+="{0}=\"{1}\" ".format(key1.upper(), params[key1])
                argList+="{0}=\"{1}\" ".format(key2.upper(), params[key2])
        return argList
        
    def weaveArgStr(self): 
        if self.nSeg != 1:
            argStr = ["output-file", "sft-files", "setup-file", "semi-max-mismatch", "coh-max-mismatch", "toplist-limit", "extra-statistics"]
        else:
            argStr = ["output-file", "sft-files", "setup-file", "semi-max-mismatch", "toplist-limit", "extra-statistics"]
            
        argListString = ""
        for s in argStr:
            argListString += "--{0}=$({1}) ".format(s, s.replace('-', '').upper())
            
        argListString += "--alpha=$(ALPHA)/$(DALPHA) --delta=$(DELTA)/$(DDELTA) "
        for i in range(len(self.freqParamName)):
            argListString += "--{0}=$({1})/$({2}) ".format(self.freqParamName[i], self.freqParamName[i].upper(), self.freqDerivParamName[i].upper())
        
        return argListString
    
    def analyzeResultArgStr(self): 
        argStr = ["target", "obsDay", "cohDay", "stage", "freq", "freqDerivOrder", "numTopList", "df1dot", "cluster"]
        argListString = ""
        for s in argStr:
            argListString += "--{0}=$({1}) ".format(s, s.replace('-', '').upper())

        return argListString
    
    def analyzeResultArgs(self, fmin, fmax, df1dot, cluster, OSG): 
        if OSG:
            argListString = 'TARGET="{0}" OBSDAY="{1}" COHDAY="{2}" STAGE="{3}" FREQ="{4}" FREQDERIVORDER="{5}" TOPLISTLIMIT="{6}" df1dot="{7}" CLUSTER="{8}"'.format(
                    self.target.name, self.obsDay, self.cohDay, self.stage, freq, self.freqDerivOrder, self.numTopList, df1dot, int(cluster))
        else:
            argListString = "argList=\" --targetList {0} --obsDay {1} --cohDay {2} --stage {3} --fmin {4} --fmax {5} --freqDerivOrder {6} --numTopList {7} --df1dot {8} --cluster {9}\"".format(
                self.target.name, self.obsDay, self.cohDay, self.stage, fmin, fmax, self.freqDerivOrder, self.numTopList, df1dot, int(cluster))
        return argListString
    
    def injectionArgStr(self):
        argStr = ["injections"]
        
        argListString = ""
        for s in argStr:
            argListString += "--{0}=$({1}) ".format(s, s.replace('-', '').upper())
     
        return argListString
    
    def injectionArg(self, colnames, injParam, OSG):
        #injParamStr = ""
        injParamStr = ";".join(["{0}={1}".format(col, injParam[col]) for col in colnames])
        #for col in colnames:
        #    injParamStr += "{0}={1};".format(col, injParam[col])
        if OSG:
            argList = ("INJECTIONS=\"{{{0}}}\"".format(injParamStr))
        else:
            argList = ("--injections={{{0}}}".format(injParamStr))
        return argList

    
    def writeSub(self, freq, taskName, crFiles, argStr, request_memory, OSG, OSDF):
        exe = fp.weaveExecutableFilePath()
        metric = fp.weaveSetupFilePath(self.cohTime, self.nSeg, self.freqDerivOrder)
        # call function to write .sub files for search
        subFileName = fp.condorSubFilePath(self.target, freq, taskName, self.stage)
        Path(subFileName).unlink(missing_ok=True)
        wc.writeSearchSub(subFileName, exe, crFiles[0], crFiles[1], crFiles[2], argStr, request_memory=request_memory, request_disk='2GB', OSG=OSG, OSDF=OSDF)
        return subFileName
    
# to do: move to utils
    def memoryUsage(self, stage):
        if 'search' in stage:
            if '1987'in self.target.name:
                memory = '20GB'
            else:
                memory = '15GB'
                
        elif 'follow' in stage:
            memory = '2GB'
        else:
            memory = '15GB'
        return memory
    
    
###### main use function

    def makeSearchDag(self, cohDay, paramList, numTopList, stage, freqDerivOrder, OSG=False, OSDF=False):
        if OSDF and not OSG:
            print('Are you sure you want to read SFTs from OSDF but not using OSG computing resources?')
        self.freqParamName, self.freqDerivParamName = utils.phaseParamName(freqDerivOrder)
        self.freqDerivOrder = freqDerivOrder
        self.numTopList = numTopList
        self.stage = stage
        self.cohDay, self.cohTime, self.nSeg, self.obsTime, self.refTime = utils.getTimeSetup(self.target.name, self.obsDay, cohDay)
        
        
        freqList = np.array([int(f) for f in paramList.keys()])
        nJobs = np.array([paramList[f].data.size for f in paramList.keys()])
        fmin, fmax = freqList[0], freqList[-1] 
        dagNameFileName = fp.dagGroupFilePath(self.target, fmin, fmax+1, self.stage)
        Path(dagNameFileName).unlink(missing_ok=True)
        
        request_memory = self.memoryUsage(self.stage)
        
        with open(dagNameFileName, 'wt') as dagNameFile:
            for freq in tqdm(freqList[nJobs>0]):
                # call function to write .sub files for search
                taskName = utils.taskName(self.target, self.stage, self.cohDay, self.freqDerivOrder, freq)
                sftFiles = utils.sftEnsemble(freq, self.obsDay, OSDF=OSDF)
                
                dagFileName = fp.dagFilePath(freq, self.target, taskName, self.stage)
                Path(dagFileName).unlink(missing_ok=True)
                    
                crFiles = fp.condorRecordFilePath(freq, self.target, taskName, self.stage)
                utils.makeDir(crFiles)

                argStr = self.weaveArgStr()
                subFileName = self.writeSub(freq, taskName, crFiles, argStr, request_memory=request_memory, OSG=OSG, OSDF=OSDF)
                
                for jobIndex, params in enumerate(paramList[str(freq)].data, 1):
                    ######################## Argument string use to write to DAG  ########################
                    argList = self.weaveArgs(freq, params, taskName, sftFiles, jobIndex, OSG)
                    # Call function from WriteCondorFiles.py which will write DAG 
                    wc.writeSearchDag(dagFileName, taskName, subFileName, jobIndex, argList)
                if len(paramList[str(freq)].data) != 0:
                    dagNameFile.write('{0}\n'.format(dagFileName))
        print('Finish writing {0} dag files for {1}-{2}Hz'.format(self.stage, fmin, fmax+1))
        
        
    def makeAnalyzeDag(self, cohDay, fmin, fmax, numTopList=1000, df1dot=1.5e-9, stage='search', freqDerivOrder=2, cluster=False):
        self.freqParamName, self.freqDerivParamName = utils.phaseParamName(freqDerivOrder)
        self.freqDerivOrder = freqDerivOrder
        self.numTopList = numTopList
        self.stage = stage
        self.cohDay, self.cohTime, self.nSeg, self.obsTime, self.refTime = utils.getTimeSetup(self.target.name, self.obsDay, cohDay)
        
        
        taskName = 'AnalyzeResult_{0}-{1}Hz'.format(fmin, fmax)
        request_memory = '1GB'
        argStr = self.analyzeResultArgStr()
        exe = fp.analyzeResultExecutableFilePath()
        subFileName = fp.condorSubFilePath(self.target, taskName, taskName, self.stage)
        Path(subFileName).unlink(missing_ok=True)
        
        crFiles = fp.condorRecordFilePath(taskName, self.target, taskName, self.stage)
        utils.makeDir(crFiles)
        wc.writeSearchSub(subFileName, exe, crFiles[0], crFiles[1], crFiles[2], argStr, request_memory=request_memory, request_disk='1GB', OSG=False)
        
        dagFileName = fp.dagFilePath(taskName, self.target, taskName, self.stage)
        Path(dagFileName).unlink(missing_ok=True)
        for jobIndex, freq in tqdm(enumerate(range(fmin, fmax), 1)):
            # call function to write .sub files for analyze result
            ######################## Argument string use to write to DAG  ########################
            argList = self.analyzeResultArgs(freq, freq+1, df1dot, cluster=cluster, OSG=False)
            # Call function from WriteCondorFiles.py which will write DAG 
            wc.writeSearchDag(dagFileName, taskName, subFileName, jobIndex, argList)
        print('Finish writing {0} dag files for {1}-{2}Hz'.format(self.stage, fmin, fmax))
    
    def makeInjectionDag(self, cohDay, fmin, fmax, paramList, injParamList, numTopList=1000, stage='search', freqDerivOrder=2, injFreqDerivOrder=4, OSG=False, OSDF=False):
        if OSDF and not OSG:
            print('Are you sure you want to read SFTs from OSDF but not using OSG computing resources?')
        self.freqParamName, self.freqDerivParamName = utils.phaseParamName(freqDerivOrder)
        self.freqDerivOrder = freqDerivOrder
        self.numTopList = numTopList
        self.stage = stage
        self.cohDay, self.cohTime, self.nSeg, self.obsTime, self.refTime = utils.getTimeSetup(self.target.name, self.obsDay, cohDay)
        injFreqParamName, _ = utils.phaseParamName(injFreqDerivOrder)
        self.injParamName = utils.injParamName() + injFreqParamName[1:]
        
        
        freqList = np.array([f for f in paramList.keys()])
        nJobs = np.array([paramList[f].data.size for f in paramList.keys()])
        dagNameFileName = fp.dagGroupFilePath(self.target, fmin, fmax, self.stage)
        Path(dagNameFileName).unlink(missing_ok=True)
        
        request_memory = '2GB'
        
        with open(dagNameFileName, 'wt') as dagNameFile:
            for freq in tqdm(freqList[nJobs>0]):
                # call function to write .sub files for search
                taskName = utils.taskName(self.target, self.stage, self.cohDay, self.freqDerivOrder, freq)
                sftFiles = utils.sftEnsemble(freq, self.obsDay, OSDF=OSDF)
                
                dagFileName = fp.dagFilePath(freq, self.target, taskName, self.stage)
                Path(dagFileName).unlink(missing_ok=True)
                
                crFiles = fp.condorRecordFilePath(freq, self.target, taskName, self.stage)
                utils.makeDir(crFiles)

                argStr = self.weaveArgStr() + self.injectionArgStr()
                subFileName = self.writeSub(freq, taskName, crFiles, argStr, request_memory=request_memory, OSG=OSG, OSDF=OSDF)
                
                #injParamName = injParamList[str(freq)].columns.names
                injParamName = self.injParamName
                for jobIndex, (searchParam, injParam) in enumerate(zip(paramList[str(freq)].data, injParamList[str(freq)].data), 1):
                    ######################## Argument string use to write to DAG  ########################
                    if not OSG:
                        argList = self.weaveArgs(freq, searchParam, taskName, sftFiles, jobIndex, OSG)[:-1]
                        argList += self.injectionArg(injParamName, injParam, OSG) + '\"'
                    else:
                        argList = self.weaveArgs(freq, searchParam, taskName, sftFiles, jobIndex, OSG) + self.injectionArg(injParamName, injParam, OSG)
                    # Call function from WriteCondorFiles.py which will write DAG 
                    wc.writeSearchDag(dagFileName, taskName, subFileName, jobIndex, argList)
                if len(paramList[str(freq)].data) != 0:
                    dagNameFile.write('{0}\n'.format(dagFileName))
        print('Finish writing {0} dag files for {1}-{2}Hz'.format(self.stage, fmin, fmax))
        return dagNameFileName
        