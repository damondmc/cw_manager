from astropy.io import fits
import numpy as np
from ..utils import filePath as fp
from ..utils import setup_parameter as setup
from ..utils import utils as utils
from pathlib import Path
import warnings
import subprocess
import time
from multiprocessing import Pool, cpu_count    

# Function to process each search job in parallel
def searchJob(params, sftFiles, metric, semiMM, cohMM, numTopList, extraStats, ra, dec, nc, nf, obsDay):
    
    # Weave main program path
    weave_exe = fp.weaveExecutableFilePath()
    
    resultFile, param = params
    utils.makeDir([resultFile])
    #print(resultFile)
    if Path(resultFile).exists():
        return resultFile
    
    command = '{} --output-file={} --sft-files=\"{}\" --setup-file={} --semi-max-mismatch={} --coh-max-mismatch={} --toplist-limit={} --extra-statistics={} --alpha={} --delta={}'.format(
        weave_exe, resultFile, sftFiles, metric, semiMM, cohMM, numTopList, extraStats, ra, dec)
    if nc == obsDay:
        command = '{} --output-file={} --sft-files=\"{}\" --setup-file={} --semi-max-mismatch={} --toplist-limit={} --extra-statistics={} --alpha={} --delta={}'.format(
            weave_exe, resultFile, sftFiles, metric, semiMM, numTopList, extraStats, ra, dec)

    newFreqParamName, newFreqDerivParamName = utils.phaseParamName(nf)
    for _f, _df in zip(newFreqParamName, newFreqDerivParamName):
        command += ' --{}={}/{}'.format(_f, param[_f], param[_df])

    # Run the command
    result = subprocess.run(command, shell=True, capture_output=True, text=True)

    # Print the standard output and errors for debugging
    #print(result.stdout)
    #print(result.stderr)
    return resultFile


# Function to process each injection job in parallel
def injectionJob(params, inj, sftFiles, metric, semiMM, cohMM, numTopList, extraStats, ra, dec, nc, nf, obsDay):

    # Weave main program path
    weave_exe = fp.weaveExecutableFilePath()
    
    resultFile, param = params
    utils.makeDir([resultFile])
    #print(resultFile)
    if Path(resultFile).exists():
        return resultFile
    
    command = '{} --output-file={} --sft-files=\"{}\" --setup-file={} --semi-max-mismatch={} --coh-max-mismatch={} --toplist-limit={} --extra-statistics={} --alpha={} --delta={}'.format(
        weave_exe, resultFile, sftFiles, metric, semiMM, cohMM, numTopList, extraStats, ra, dec)
    if nc == obsDay:
        command = '{} --output-file={} --sft-files=\"{}\" --setup-file={} --semi-max-mismatch={} --toplist-limit={} --extra-statistics={} --alpha={} --delta={}'.format(
            weave_exe, resultFile, sftFiles, metric, semiMM, numTopList, extraStats, ra, dec)

    newFreqParamName, newFreqDerivParamName = utils.phaseParamName(nf)
    for _f, _df in zip(newFreqParamName, newFreqDerivParamName):
        command += ' --{}={}/{}'.format(_f, param[_f], param[_df])
           
    injection_command = ' --injections=\"{{Alpha={};Delta={};refTime={};aPlus={};aCross={};psi={};Freq={};f1dot={};f2dot={};f3dot={};f4dot={}}}\"'.format(
    inj['Alpha'], inj['Delta'], inj['refTime'], 
    inj['aPlus'], inj['aCross'], inj['psi'], 
    inj['Freq'], inj['f1dot'], inj['f2dot'], 
    inj['f3dot'], inj['f4dot'])
    
    command += injection_command
    #print(command)
    # Run the command
    result = subprocess.run(command, shell=True, capture_output=True, text=True)

    # Print the standard output and errors for debugging
    #print(result.stdout)
    #print(result.stderr)
    return resultFile

def determineEfficiency(sftFiles, setup, cohDay, obsDay, im, rm, h0est, target, taskName, freq, nInj, freqDerivOrder, stage, numTopList, extraStats, num_cpus, cluster, workInLocalDir, saveIntermediate=False):
    
    sp, ip = im._genParam(h0=h0est, freq=freq, nInj=nInj, injFreqDerivOrder=4, freqDerivOrder=freqDerivOrder, skyUncertainty=1e-5, workInLocalDir=workInLocalDir)

    if workInLocalDir:
        search_params = [(Path(fp.weaveOutputFilePath(target, freq, taskName, jobIndex, stage)).name, params) for jobIndex, params in enumerate(sp[str(freq)].data, 1)]
    else:
        search_params = [(fp.weaveOutputFilePath(target, freq, taskName, jobIndex, stage), params) for jobIndex, params in enumerate(sp[str(freq)].data, 1)]
    inj_params = ip[str(freq)].data

    print("Generated params for h0={}, running Weave...".format(h0est)) 
    injResultFileList = [] 
    
    _, cohTime, nSeg, _, _ = utils.getTimeSetup(target.name, obsDay, cohDay)
    metric = fp.weaveSetupFilePath(cohTime, nSeg, freqDerivOrder)
    if workInLocalDir:  
        metric = Path(metric).name

    with Pool(processes=num_cpus) as pool:
        results = pool.starmap(injectionJob, [(params, inj, sftFiles, metric, setup.semiMM, setup.cohMM, numTopList, extraStats, target.alpha, target.delta, cohDay, freqDerivOrder, obsDay) for params, inj in zip(search_params, inj_params)])
    # Collect the results
    injResultFileList.extend(results)

    outlierFilePath = rm.writeInjectionResult1Hz(cohDay, freq, nInj, numTopList=numTopList, stage=stage, freqDerivOrder=freqDerivOrder, cluster=cluster, workInLocalDir=workInLocalDir)

    if not saveIntermediate:
        # delete files to release disk storage
        for _f in injResultFileList:
            command = 'rm {}'.format(_f)
            _ = subprocess.run(command, shell=True, capture_output=True, text=True)
        print('Deleted weave result files to release disk storage.\n')

    nout = fits.getdata(outlierFilePath,1).size
    p = nout / nInj
    print('{}% ({}/{}) above mean2F threshold, saved to {}.'.format(round(p*100,2), nout, nInj, outlierFilePath))
    return p, outlierFilePath

def injectionFollowUp(fm, rm, target, obsDay, freq, sftFiles, 
                      old_cohDay, old_freqDerivOrder, old_stage, new_cohDay, 
                      new_freqDerivOrder, new_stage, nInj, numTopList, extraStats, num_cpus, setup, 
                      cluster, workInLocalDir, saveIntermediate=False):
    print('Doing injection follow-up...')    
    
    sp, ip = fm.genFollowUpParamFromInjection1Hz(old_cohDay, freq, stage='inj_'+old_stage, oldFreqDerivOrder=old_freqDerivOrder, newFreqDerivOrder=new_freqDerivOrder, cluster=cluster, workInLocalDir=workInLocalDir)

    if sp[str(freq)].data.size == 0:
        print('0 outliers from the previous injection stage: Error!')
        exit()
    else:
        print('{} injections to be carried out.'.format(sp[str(freq)].data.size))

    injResultFileList = []
    taskName = utils.taskName(target, 'inj_'+new_stage, new_cohDay, new_freqDerivOrder, freq)
    cohDay, cohTime, nSeg, _, _ = utils.getTimeSetup(target.name, obsDay, new_cohDay)
    metric = fp.weaveSetupFilePath(cohTime, nSeg, new_freqDerivOrder)
    if workInLocalDir:  
        metric = Path(metric).name

    # Prepare job parameters for parallel execution
    if workInLocalDir:
        search_params = [(Path(fp.weaveOutputFilePath(target, freq, taskName, jobIndex, 'inj_'+new_stage)).name, params) for jobIndex, params in enumerate(sp[str(freq)].data, 1)]
    else:
        search_params = [(fp.weaveOutputFilePath(target, freq, taskName, jobIndex, 'inj_'+new_stage), params) for jobIndex, params in enumerate(sp[str(freq)].data, 1)]

    inj_params = ip[str(freq)].data

    print("Generated params, running Weave...")

    # Use multiprocessing to process the jobs in parallel
    with Pool(processes=num_cpus) as pool:
        results = pool.starmap(injectionJob, [(params, inj, sftFiles, metric, setup.semiMM, setup.cohMM, 1000, extraStats, target.alpha, target.delta, new_cohDay, new_freqDerivOrder, obsDay) for params, inj in zip(search_params, inj_params)])
    #for result in results:
    #    print(result)        
    # Collect the results
    injResultFileList.extend(results)

    print('Analyzing injection result.')
    # analyze the result
    outlierFilePath = rm.writeFollowUpResult(old_cohDay, new_cohDay, freq, numTopList=numTopList, old_stage='inj_'+old_stage, new_stage='inj_'+new_stage, old_freqDerivOrder=old_freqDerivOrder, new_freqDerivOrder=new_freqDerivOrder, workInLocalDir=workInLocalDir, inj=True, cluster=cluster)
    #inj_outlierFilePathList.append(_outlierFilePath)

    # delete files to release disk storage
    if not saveIntermediate:
        for _f in injResultFileList:
            command = 'rm {}'.format(_f)
            _ = subprocess.run(command, shell=True, capture_output=True, text=True)
        print('Deleted weave result files to release disk storage.\n')
     
    return outlierFilePath


def realFollowUp(rm, sp, target, obsDay, freq, sftFiles, mean2F_ratio,
                 old_cohDay, old_freqDerivOrder, old_stage, 
                 new_cohDay, new_freqDerivOrder, new_stage, 
                 numTopList, extraStats, num_cpus, setup, 
                 cluster, workInLocalDir, saveIntermediate=False):
    print('Doing real follow-up...')
     
    searchResultFileList = []
    taskName = utils.taskName(target, new_stage, new_cohDay, new_freqDerivOrder, freq)
    cohDay, cohTime, nSeg, _, _ = utils.getTimeSetup(target.name, obsDay, new_cohDay)
    metric = fp.weaveSetupFilePath(cohTime, nSeg, new_freqDerivOrder)
    if workInLocalDir:  
        metric = Path(metric).name

    # Prepare job parameters for parallel execution
    if workInLocalDir:
        search_params = [(Path(fp.weaveOutputFilePath(target, freq, taskName, jobIndex, new_stage)).name, params) for jobIndex, params in enumerate(sp[str(freq)].data, 1)]
    else:
        search_params = [(fp.weaveOutputFilePath(target, freq, taskName, jobIndex, new_stage), params) for jobIndex, params in enumerate(sp[str(freq)].data, 1)]

    # Use multiprocessing to process the jobs in parallel

    print("Generated params, running Weave...")
    #resultFile = Path(resultFile).name
    with Pool(processes=num_cpus) as pool:
        results = pool.starmap(searchJob, [(params, sftFiles, metric, setup.semiMM, setup.cohMM, 1000, extraStats, target.alpha, target.delta, new_cohDay, new_freqDerivOrder, obsDay) for params in search_params])
    #for result in results:
    #    print(result)        
    # Collect the results
    searchResultFileList.extend(results)

    # analyze the result
    outlierFilePath = rm.writeFollowUpResult(old_cohDay, new_cohDay, freq, numTopList=numTopList, old_stage=old_stage, new_stage=new_stage, old_freqDerivOrder=old_freqDerivOrder, new_freqDerivOrder=new_freqDerivOrder, ratio=mean2F_ratio, workInLocalDir=workInLocalDir, inj=False, cluster=cluster)
    #outlierFilePathList.append(_outlierFilePath)

    print('Finished analysis for current stage.')
    if not saveIntermediate:
        # delete files to release disk storage
        for _f in searchResultFileList:
            command = 'rm {}'.format(_f)
            _ = subprocess.run(command, shell=True, capture_output=True, text=True)
        print('Deleted weave result files to release disk storage.')

    return outlierFilePath


def determineMean2FRatio(percentile, target, freq, 
                         old_cohDay, old_freqDerivOrder, old_stage, 
                         new_cohDay, new_freqDerivOrder, new_stage, 
                         cluster=False, workInLocalDir=False):
    taskName = utils.taskName(target=target, stage=old_stage, cohDay=old_cohDay, order=old_freqDerivOrder, freq=freq)
    filePath = fp.outlierFilePath(target, freq, taskName, old_stage, cluster=cluster)
    if workInLocalDir:
        filePath = Path(filePath).name
    olddata = fits.getdata(filePath, 1)

    taskName = utils.taskName(target=target, stage=new_stage, cohDay=new_cohDay, order=new_freqDerivOrder, freq=freq)
    filePath = fp.outlierFilePath(target, freq, taskName, new_stage, cluster=cluster)
    if workInLocalDir:
        filePath = Path(filePath).name
    
    newdata = fits.getdata(filePath, 1)
    ratio = np.sort(newdata['mean2F']/olddata['mean2F'])
    r = np.percentile(ratio, percentile) ## 99.5-99.6 percentile
    r = int(r*100.)/100.
    print('Ratio = {} for {}% percentile for {} injections.\n'.format(r, (1-percentile)*100, ratio.size))
    return r 
