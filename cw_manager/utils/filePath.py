# This file contains all file path used in the Weave-based SNR search pipeline
# Please put all file path in this file for better maintanace in the future
from . import setup_parameter as setup


############################################ Core file
# path for Weave main program in igwn grid (OSG accessible)
def weaveExecutableFilePath():
    filePath = '/cvmfs/software.igwn.org/conda/envs/igwn-py39-20231212/bin/lalpulsar_Weave'
    return filePath

# path for python main program for the followUp process
def followUpExecutableFilePath():
    filePath = setup.homeDir + 'followUp.py'
    return filePath

# path for python main program for the injectionFollow process
def injFollowUpExecutableFilePath():
    filePath = setup.homeDir + 'followInjection.py'
    return filePath

# path for python main program for the upper limit determination process
def upperLimitExecutableFilePath():
    filePath = setup.homeDir + 'upperLimit.py'
    return filePath

# metric file for Weave 
def weaveSetupFilePath(cohTime, nSeg, freqOrder):
    filePath = setup.homeDir+'metricSetup/Start{0}_TCoh{1}_N{2}_Spin{3}.fts'.format(setup.startTime, cohTime, nSeg, freqOrder)
    return filePath

def analyzeResultExecutableFilePath():
    filePath = setup.homeDir+'/scripts/new/analyze.py'
    return filePath

# sft files
def sftFilePath(obsDay, freq, detector='H1', OSDF=False):
    if OSDF:
        rootDir = setup.OSDFDir + '/SFTs/o4a_data/'
    else:
        #rootDir = setup.homeDir
        rootDir = '/home/'+setup.user+'/SFTs/o4a_data/'
        
    if detector == 'H1':
        filePath = rootDir+'SFTs/narrowBand_age300yr/{0}days/H1/{1}/'.format(obsDay, int(freq))
    elif detector == 'L1':
        filePath = rootDir+'SFTs/narrowBand_age300yr/{0}days/L1/{1}/'.format(obsDay, int(freq))
    return filePath

def estimateUpperLimitExcutable():
    filePath = '/cvmfs/software.igwn.org/conda/envs/igwn-py39-20231212/bin/lalpulsar_ComputeFstatMCUpperLimit'
    return filePath
    
############################################ Condor-related

# file that store all dag filenames
def dagGroupFilePath(target, fmin, fmax, stage):
    filePath = setup.homeDir + 'dagJob/{0}_{1}_{2}-{3}Hz_dagFiles.txt'.format(target.name, stage, fmin, fmax)
    return filePath
    
# dag file that store the variable to submit condor jobs
def dagFilePath(freq, target, taskName, stage):
    filePath = setup.homeDir + 'condorFiles/{0}/{1}/{2}/{3}.dag'.format(stage, target.name, freq, taskName) 
    return filePath
    
# path for condor submit file
def condorSubFilePath(target, freq, taskName, stage):
    filePath = setup.homeDir + 'condorFiles/{0}/{1}/{2}/{3}.sub'.format(stage, target.name, freq, taskName)
    return filePath

# path for condor to submit another DAG
def SubmitCondorSubFilePath(target, freq, stage):
    filePath = setup.homeDir + 'condorFiles/{0}/{1}/{2}/{3}.sub'.format(stage, target.name, freq, 'submit')
    return filePath

# condor job records 
def condorRecordFilePath(freq, target, taskName, stage):
    outputFile = setup.homeDir + 'results/{0}/{1}/{2}/{3}/OUT/{4}.out.$(JobID)'.format(stage, target.name, setup.sftSource, freq, taskName)
    errorFile = setup.homeDir + 'results/{0}/{1}/{2}/{3}/ERR/{4}.err.$(JobID)'.format(stage, target.name, setup.sftSource, freq, taskName)
    logFile = setup.homeDir + 'results/{0}/{1}/{2}/{3}/LOG/{4}_Log.txt.$(JobID)'.format(stage, target.name, setup.sftSource, freq, taskName)
    return [outputFile, errorFile, logFile]

############################################ Weave output file

# path to save the Weave output
def weaveOutputFilePath(target, freq, taskName, jobIndex, stage):
    filePath = setup.homeDir + 'results/{0}/{1}/{2}/{3}/Result/{4}.fts.{5}'.format(stage, target.name, setup.sftSource, freq, taskName, jobIndex)
    return filePath
    
#  checkpoint file path for Weave
def checkPointFilePath(target, freq, jobIndex, stage):
    filePath = setup.homeDir + 'results/{0}/{1}/{2}/{3}/Result/Checkpoint.fts.{4}'.format(stage, target.name, setup.sftSource, freq, jobIndex)
    return filePath

# file to save the outlier after analyzing the weave result file
def outlierFilePath(target, freq, taskName, stage, cluster=False):
    if not cluster:
         filePath = setup.homeDir + 'results/{0}/{1}/{2}/{3}/Outliers/{4}_outlier.fts'.format(stage, target.name, setup.sftSource, freq, taskName)
    else:
         filePath = setup.homeDir + 'results/{0}/{1}/{2}/{3}/Outliers/{4}_outlier_clustered.fts'.format(stage, target.name, setup.sftSource, freq, taskName)
    return filePath

# file to save the outlier after analyzing the weave result file
def outlierFromSaturatedFilePath(target, freq, taskName, stage):
    filePath = setup.homeDir + 'results/{0}/{1}/{2}/{3}/Outliers/{4}_LoudestOutlierFromSaturated.fts'.format(stage, target.name, setup.sftSource, freq, taskName)
    return filePath

############################################ Estimate upper limit output file
def estimateUpperLimitFilePath(target, freq, taskName, stage):
    filePath = setup.homeDir + 'results/{0}/{1}/{2}_{3}Hz.txt'.format(stage, target.name, taskName, freq)
    return filePath

############################################ Info & summary file 
# file to save the the meta info of outlier after analyzing the weave result file (no. of outlier per sub-band)
def outlierInfoFilePath(target, freq, taskName, stage, cluster=False):
    if not cluster:
        filePath = setup.homeDir + 'results/{0}/{1}/{2}/{3}/Outliers/{4}_info.txt'.format(stage, target.name, setup.sftSource, freq, taskName)
    else:
        filePath = setup.homeDir + 'results/{0}/{1}/{2}/{3}/Outliers/{4}_clustered_info.txt'.format(stage, target.name, setup.sftSource, freq, taskName)      
    return filePath 

def imageFilePath(OSDF=False):
    if OSDF:    
        filePath = 'osdf:///igwn/cit/staging/hoitim.cheung/images/'
    else:
        filePath = '/home/hoitim.cheung/gc/cw_manager/'
    filePath += 'cw_manager_gc.sif'
    return filePath
