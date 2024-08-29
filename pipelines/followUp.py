#!/opt/conda/bin/python
import numpy as np
import argparse
import importlib
from tqdm import tqdm
import subprocess
from pathlib import Path
import time

import cw_manager 
from cw_manager.utils import setup_parameter as setup
from cw_manager.utils import utils as utils
from cw_manager.utils import filePath as fp
from cw_manager.analysis.resultManager import resultManager
from cw_manager.genParam.followUpParam1Hz import followUpParams

parser = argparse.ArgumentParser(description='Plot upper strain limit for a certain target base on injection test.')
parser.add_argument('--target', type=str, help='Target to be analyzed', default='CassA')
parser.add_argument('--obsDay',type=float,help='Total observation time',default=240)  
parser.add_argument('--cohDay',type=float,help='Coherence time for the previous stage.',default=5) 
parser.add_argument('--freq',type=int,help='Frequency band to be analyzed.',default=20)
parser.add_argument('--stage',type=str,help='Previous stage of the pipeline.',default='search') 
parser.add_argument('--sftFiles',type=str,help='SFT filepaths being used.',default='./') 
parser.add_argument('--freqDerivOrder',type=int,help='Highest order used in the phase evolution for the previous stage.',default=2)
parser.add_argument('--numTopList',type=int,help='Number of loudest templates returned by Weave',default=10)
parser.add_argument('--workInLocalDir',action='store_true',help='Work in local directory.')
parser.add_argument('--inj',action='store_true',help='To indicate whether injection is included.')
parser.add_argument('--cluster',action='store_true',help='Enable clustering for the outliers.')
parser.add_argument('--OSG',action='store_true',help='Enable the use of OSG.')
parser.add_argument('--OSDF',action='store_true',help='Enable the use of OSDF.')

args = parser.parse_args()
obsDay = args.obsDay
old_cohDay = args.cohDay
if old_cohDay.is_integer(): 
    old_cohDay = int(old_cohDay)
freq = args.freq
old_stage = args.stage
old_freqDerivOrder = args.freqDerivOrder
sftFiles = args.sftFiles
numTopList = args.numTopList
workInLocalDir = args.workInLocalDir
cluster = args.cluster
inj = args.inj
OSG = args.OSG
OSDF = args.OSDF

t0 = time.time()

if inj:
    print('Analyzing follow-up with software injection.')

target = importlib.import_module('cw_manager.target.'+ args.target)
print('Working on {0} ...'.format(target.name))

# define the follow-up parameter for each stage in a list
new_cohDay = [15, 30, 60, 120, 240]
new_freqDerivOrder = [2, 2, 3, 4, 4]
new_stage = ['followUp-{}'.format(i+1) for i in range(len(new_cohDay))]

# Weave main program path
weave_exe = fp.weaveExecutableFilePath()
#weave_exe = Path(weave_exe).name
# Weave extra output setup
extraStats = "coh2F_det,mean2F,coh2F_det,mean2F_det"

# initialize the manager class
fm = followUpParams(target, obsDay)
rm = resultManager(target, obsDay=obsDay)

##############
oc = old_cohDay
of = old_freqDerivOrder
os = old_stage
print(sftFiles)
# save all the outlier.fts files
outlierFilePathList = []
# start the follow-up process for-loop
for nc, nf, ns in zip(new_cohDay, new_freqDerivOrder, new_stage):
    print('{} Hz {}:'.format(freq, os))
    sp = fm.genFollowUpParam(
                oc, freq, stage=os, oldFreqDerivOrder=of, newFreqDerivOrder=nf, workInLocalDir=workInLocalDir)

    if sp[str(freq)].data.size == 0:
        print('Go to final process: combine all outlier files into one.')
        break
    else:
        print('{} outliers to follow-up...'.format(sp[str(freq)].data.size))
    
    taskName = utils.taskName(target, ns, nc, nf, freq)
    cohDay, cohTime, nSeg, _, _ = utils.getTimeSetup(target.name, obsDay, nc)
    metric = fp.weaveSetupFilePath(cohTime, nSeg, nf)
    metric = Path(metric).name
    
    # loop over every jobs in each Hz
    for jobIndex, params in enumerate(sp[str(freq)].data, 1):
        resultFile = fp.weaveOutputFilePath(target, freq, taskName, jobIndex, ns)
        resultFile = Path(resultFile).name
        if Path(resultFile).exists():
            continue
        command = '{} --output-file={} --sft-files=\"{}\" --setup-file={} --semi-max-mismatch={} --coh-max-mismatch={} --toplist-limit={} --extra-statistics={} --alpha={} --delta={}'.format(
            weave_exe, resultFile, sftFiles, metric, setup.semiMM, setup.cohMM, numTopList, extraStats, target.alpha, target.delta)
        
        newFreqParamName, newFreqDerivParamName = utils.phaseParamName(nf)
        for _f, _df in zip(newFreqParamName, newFreqDerivParamName):
            command += ' --{}={}/{}'.format(_f, params[_f], params[_df])
            
        print(command)
        # Run the command
        result = subprocess.run(command, shell=True, capture_output=True, text=True)
        
        # Print the standard output
        print(result.stdout)
        print(result.stderr)
        
    # analyze the result 
    _outlierFilePath = rm.writeFollowUpResult(oc, nc, freq, numTopList=numTopList, old_stage=os, new_stage=ns, old_freqDerivOrder=of, new_freqDerivOrder=nf, workInLocalDir=workInLocalDir, inj=inj, cluster=cluster)
    outlierFilePathList.append(_outlierFilePath)
    
    # update parameters
    oc = nc
    of = nf
    os = ns
    print('Finish {} {} for {} Hz'.format(target.name, ns, freq))
    

# combine all outlier files into one file
if len(outlierFilePathList) != 0:
    taskName = utils.taskName(target, 'followUp', old_cohDay, old_freqDerivOrder, freq)
    finalOutlierPath = rm.emsembleFollowUpResult(outlierFilePathList, freq, taskName, old_stage, workInLocalDir=workInLocalDir, cluster=cluster)
    print('Final result saved as', finalOutlierPath)
else:
    print('No outlier starting from the beginning.')
    
print('Finish all follow-up process for {}'.format(target.name))
print('Time used = {} s'.format(time.time()-t0))
