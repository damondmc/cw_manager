############# basic functions 
import numpy as np
from . import setup_parameter as setup
from pathlib import Path
from . import filePath as fp
from astropy.io import fits

def sftEnsemble(freq, obsDay, OSDF=False):
    H1path = Path(fp.sftFilePath(obsDay, freq, detector='H1', OSDF=OSDF))
    L1path = Path(fp.sftFilePath(obsDay, freq, detector='L1', OSDF=OSDF))
    if not OSDF:
        lst = [str(s) for s in H1path.glob("*.sft")] + [str(s) for s in L1path.glob("*.sft")]
    else:
        lst = ['osdf://'+str(s)[5:] for s in H1path.glob("*.sft")] + ['osdf://'+str(s)[5:] for s in L1path.glob("*.sft")]     
    return lst

def makeDir(filenames):
    for name in filenames:
        Path(name).resolve().parent.mkdir(parents=True, exist_ok=True)
    
def taskName(target, stage, cohDay, order, freq):
    tn='{0}_{1}_TCoh{2}_O{3}_{4}Hz'.format(
        target.name, stage, cohDay, order, freq)
    return tn

def phaseParamName(order):
    freqParamName = ["freq", "f1dot", "f2dot", "f3dot", "f4dot"]
    freqDerivParamName = ["df", "df1dot", "df2dot", "df3dot", "df4dot"]        
    return freqParamName[:order+1], freqDerivParamName[:order+1]


def injParamName():
    #return ["Alpha", "Delta", "refTime", "h0", "cosi", "psi", "Freq"]
    return ["Alpha", "Delta", "refTime", "aPlus", "aCross", "psi", "Freq"]
    


def memoryUsage(self, stage):
    if 'search' in stage:
        memory = '10GB'
                
    elif 'follow' in stage:
        memory = '2GB'
    return memory
    
    
def getSpacing(dataFilePath, freqDerivOrder):
    
    file = fits.open(dataFilePath)
    metaData = file[0].header
    file.close()
    freqParamName, freqBandWidthName = phaseParamName(freqDerivOrder)
    n = len(freqParamName)

    nTemp_cum = [metaData['NSEMITMPL NU{0}DOT'.format(i)] for i in range(n)]
    nTemp = []
    nTemp.append(int(nTemp_cum[0]/nTemp_cum[-1]))           # avg. no. of templates for f0
    nTemp.append(nTemp_cum[1])                          # avg. no. of templates for f1
    for i in range(2,len(freqParamName)):
        nTemp.append(int(nTemp_cum[i]/nTemp_cum[i-1]))       # avg. no. of templates for f2 - fn

    df = {}
    for i in range(n):
        a,b = metaData['PROGARG {0}'.format(freqParamName[i]).upper()].split(',')
        bandwidth = float(b) - float(a)
        df[freqBandWidthName[i]] = bandwidth/nTemp[i]
    return df

def getTimeSetup(source, obsDay, cohDay):
    cohTime = int(cohDay*setup.secondsInDay)
    nSeg = int(obsDay/cohDay)
    obsTime = cohTime*nSeg
    refTime = int(setup.startTime + obsTime/2)
    return cohDay, cohTime, nSeg, obsTime, refTime


def genh0Points(i, h0, nInj, nAmp):
    i = np.atleast_1d(i)
    if nAmp != 1:
        injPerPoint = int (nInj/nAmp)
        factor = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9] # Vela v1 and G347 (new)
        #factor = [0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8] # (new2)
        return h0 * np.array([factor[int(ii/injPerPoint)] for ii in i] )
    else:
        return h0

def getNonSaturatedBand(target, fmin, fmax, nBands=None, complete=False, saveFreqList=True):
    filePath = fp.infoSummaryFilePath(target, fmin, fmax, 'search')
    freq, _, nsat = np.loadtxt(filePath).T
    
    if nBands is not None:
        freqList = freq[nsat==0]
        # randomly sample nBands:
        if len(freqList) > nBands:
            np.random.shuffle(freqList)
            nonSatBands = freqList[:nBands]
        else:
            nonSatBands = np.array([freqList[np.random.randint(0, len(freqList), 1)] for _ in range(nBands)])

        saveFilePath = fp.nonSaturatedBandFilePath(target, fmin, fmax, nBands, 'search')
        if saveFreqList:
            with open(saveFilePath ,'w') as file:
                file.write('#freq(Hz)\n')
                for nonSatFreq in nonSatBands:
                    file.write('{0}\n'.format(nonSatFreq))

    elif complete:
        # return all band:
        saveFilePath = fp.nonSaturatedBandFilePath(target, fmin, fmax, nBands, 'search')
        nonSatBands = freq[nsat==0]
        if saveFreqList:
            with open(saveFilePath ,'w') as file:
                file.write('#freq(Hz)\n')
                for nonSatFreq in nonSatBands:
                    file.write('{0}\n'.format(nonSatFreq))  
   
    else:
        # return 1 sample per each 1Hz band:
        nonSatBands =np.zeros(fmax-fmin)

        for i, f0 in enumerate(range(fmin, fmax)):
            idx = np.where((freq.astype(int)==f0)*(nsat==0))[0]
            freqList = freq[idx]
            if len(freqList) != 0:
                nonSatBands[i] = freqList[0]

        saveFilePath = fp.nonSaturatedBandFilePath(target, fmin, fmax, nBands, 'search')
        if saveFreqList:
            with open(saveFilePath ,'w') as file:
                file.write('#freq(Hz)\n')
                for nonSatFreq in nonSatBands:
                    file.write('{0}\n'.format(nonSatFreq))

    return nonSatBands 
        

def loadNonSaturatedBand(target, fmin, fmax, nBands=None):
    saveFilePath = fp.nonSaturatedBandFilePath(target, fmin, fmax, nBands, 'search')
    nonSatBands = np.loadtxt(saveFilePath).T
    return nonSatBands

def getSaturatedBand(target, fmin, fmax, saveFreqList=True):
    filePath = fp.infoSummaryFilePath(target, fmin, fmax, 'search')
    freq, _, nsat = np.loadtxt(filePath).T
    
    saveFilePath = fp.saturatedBandFilePath(target, fmin, fmax, nBands, 'search')
    satBands = freq[nsat!=0]
    if saveFreqList:
        with open(saveFilePath ,'w') as file:
            file.write('#freq(Hz)\n')
            for satFreq in satBands:
                file.write('{0}\n'.format(satFreq))  
    return nonSatBands 
        

def loadSaturatedBand(target, fmin, fmax, nBands=None):
    saveFilePath = fp.saturatedBandFilePath(target, fmin, fmax, nBands, 'search')
    satBands = np.loadtxt(saveFilePath).T
    return satBands


def readOutlierData(target, freq, cohDay, freqDerivOrder, stage, cluster=False):
    taskName = utils.taskName(target, stage, cohDay, freqDerivOrder, freq)
    outlierFilePath = fp.outlierFilePath(target, freq, taskName, stage, cluster)            
    hdul = fits.getdata(outlierFilePath)
    data = hdul[1].data
    hdul.close()
    return data

def followUpMean2F_ratio(target, stage):
    
    if target.name == 'CassA':
        # 99.5% for 20000 injections
        followUpRatio = {
            "followUp-1": 1.3,
            "followUp-2": 1.35,
            "followUp-3": 1.6,
            "followUp-4": 1.65,
            #"followUp-5": 1.,
        }          
    elif target.name == 'VelaJr': 
        # 99.5% for 20000 injections
        followUpRatio = {
            "followUp-1": 1.35,
            "followUp-2": 1.5,
            "followUp-3": 1.6,
#            "followUp-4": 1.65,
            "followUp-4": 1.45,
#            "followUp-5": 1.75,
            "followUp-5": 1.4,
        }    
    elif target.name == 'G347.3-0.5':
        # 99.5% for 20000 injections
        followUpRatio = {
            "followUp-1": 1.35,
            "followUp-2": 1.5,
            "followUp-3": 1.6,
            #"followUp-4": 1.65,
            "followUp-4": 1.5,
            #"followUp-5": 1.75,
            "followUp-5": 1.8,
        }          
    
    for key in followUpRatio.keys():
        if key in stage:
            return followUpRatio[key]
        
    print('Mean2F ratio for {0} stage is not provided.'.format(stage))
    
