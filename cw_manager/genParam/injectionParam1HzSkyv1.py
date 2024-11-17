from ..utils import utils as utils
from tqdm import tqdm
from astropy.io import fits
from astropy.table import Table
import numpy as np
from ..utils import setup_parameter as setup
from astropy.io import fits
from . import frequencyRange as fr
from ..analysis import readFile as rf
from ..utils import filePath as fp
from pathlib import Path    

class injectionParams:    
    def __init__(self, target, obsDay, cohDay, fBand=0.1):
        self.target = target
        self.cohDay= cohDay
        _, _, _, _, self.refTime = utils.getTimeSetup(self.target.name, obsDay, cohDay)
        self.fBand = fBand
        self.injParamName = utils.injParamName()
        
    def getF0FromNonSatBands(self, nonSatBands, nInj):
        p = np.ones(nonSatBands.shape)
        p = p/p.sum() # Normalize to sum up to one
        f0 = [np.random.choice([np.random.uniform(freq, freq+self.fBand) for freq in nonSatBands], p=p) for _ in range(nInj)]
        return np.array(f0).reshape(nInj)

    def genInjParamTable(self, nonSatBands, h0, freq, nInj, nAmp, freqDerivOrder, skyUncertainty):
        freqParamName, freqDerivParamName = utils.phaseParamName(freqDerivOrder)
        injData =np.recarray((nInj,), dtype=[(key, '>f8') for key in (self.injParamName+freqParamName[1:])]) 

        for i in range(nInj):
            injData[i]['psi'] = np.random.uniform(0,np.pi/4.0)
            # alpha is uniformly distributed in [0, 2pi]
            injData[i]["Alpha"] = np.random.uniform(self.target.alpha-skyUncertainty, self.target.alpha+skyUncertainty)
            # sin(delta) is uniformly distributed in [-1, 1]
            sinDelta = np.random.uniform(np.sin(self.target.delta-skyUncertainty), np.sin(self.target.delta+skyUncertainty))
            injData[i]["Delta"] = np.arcsin(sinDelta)
            injData[i]["refTime"] = self.refTime
            cosi = np.random.uniform(-1,1)
            _h0 = utils.genh0Points(i, h0, nInj, nAmp) 
            injData[i]["aPlus"] = _h0*(1.+cosi**2)/2.
            injData[i]["aCross"] = _h0*cosi
            
            # draw injection params from defined search range
            f0 = self.getF0FromNonSatBands(nonSatBands, 1) # draw from non saturated bands in 1Hz
            injData[i]['Freq'] = f0
            
            f1min, f1max, _ = fr.f1BroadRange(f0, 0, self.target.tau)
            f1 = np.random.uniform(f1min, f1max, 1)
            injData[i]['f1dot'] = f1
            
            f2min, f2max, _ = fr.f2BroadRange(f0, 0, f1, f1)
            f2 = np.random.uniform(f2min, f2max, 1)
            injData[i]['f2dot'] = f2
            
            if freqDerivOrder >= 3:
                injData[i]['f3dot'] = fr.f3Value(f0, f1, f2)
            
            if freqDerivOrder >= 4:
                injData[i]['f4dot'] = fr.f4Value(f0, f1, f2)
            
        return fits.BinTableHDU(injData)
            
    def genSearchRangeTable(self, dataFilePath, freq, injData, stage, freqDerivOrder):
        freqParamName, freqDerivParamName = utils.phaseParamName(freqDerivOrder)
        nSpacing = setup.followUp_nSpacing
        n = injData.size

        _d = fits.getheader(dataFilePath)
        spacing = {key: _d['HIERARCH ' + key] for key in freqDerivParamName}        

        data =np.recarray((n,), dtype=[(key, '>f8') for key in (freqParamName+freqDerivParamName)]) 
        for i in range(n):            
            # get frequency evolution parameters' spacing 
            taskName = utils.taskName(self.target, stage, self.cohDay, freqDerivOrder, int(freq))

            # avoid the search arange to cross the sub-band bounday and hit the saturated band
            idx1, idx2 = freqParamName[0], freqDerivParamName[0]
            eps = spacing[idx2]
            injf0 =  injData[i][idx1.capitalize()] # Weave use "Freq" not "freq" for injection's f0
            
            f0_round = np.floor(injf0*10)/10
            f0min, f0max, _ = fr.f0BroadRange(f0_round, self.fBand)
            if (injf0-nSpacing*eps)<f0min:
                data[i][idx1], data[i][idx2] = f0min, 2*nSpacing*eps
            elif (injf0+nSpacing*eps)>f0max:
                data[i][idx1], data[i][idx2] = f0max - 2*nSpacing*eps, 2*nSpacing*eps
            else:
                data[i][idx1], data[i][idx2] = injf0 - nSpacing*eps, 2*nSpacing*eps
                
            # f1dot
            f1min, f1max, _ = fr.f1BroadRange(injf0, 0, self.target.tau)
            idx1, idx2 = freqParamName[1], freqDerivParamName[1]
            eps = spacing[idx2]
            injf1 =  injData[i][idx1] 
            #if (injf1 - nSpacing * eps) < f1min:
            #    data[i][idx1], data[i][idx2] = f1min, 2*nSpacing*eps
            #elif (injf1 + nSpacing * eps) > f1max:
            #    data[i][idx1], data[i][idx2] = f1max - 2*nSpacing*eps, 2*nSpacing*eps
            #else:
            #    data[i][idx1], data[i][idx2] = injf1 - nSpacing*eps, 2*nSpacing*eps
            data[i][idx1], data[i][idx2] = injf1 - nSpacing*eps, 2*nSpacing*eps            
   
            # f2dot
            f2min, f2max, _ = fr.f2BroadRange(injf0, 0, injf1, injf1)
            idx1, idx2 = freqParamName[2], freqDerivParamName[2]
            eps = spacing[idx2]
            injf2 =  injData[i][idx1] 
            #if (injf2 - nSpacing*eps)<f2min:
            #    data[i][idx1], data[i][idx2] = f2min, 2*nSpacing*eps
            #elif (injf2 + nSpacing*eps)>f2max:
            #    data[i][idx1], data[i][idx2] = f2max - 2*nSpacing*eps, 2*nSpacing*eps
            #else:
            #    data[i][idx1], data[i][idx2] = injf2 - nSpacing*eps, 2*nSpacing*eps
            data[i][idx1], data[i][idx2] = injf2 - nSpacing*eps, 2*nSpacing*eps            
                           
            # f3dot
            if freqDerivOrder >= 3:
                f3min, f3max, _ = fr.f3BroadRange(injf0, 0, injf1, injf1, injf2, injf2)
                idx1, idx2 = freqParamName[2], freqDerivParamName[2]
                eps = spacing[idx2]
                injf3 =  injData[i][idx1] 
                #if (injf3 - nSpacing*eps)<f3min:
                #    data[i][idx1], data[i][idx2] = f3min, 2*nSpacing*eps
                #elif (injf3 + nSpacing*eps)>f3max:
                #    data[i][idx1], data[i][idx2] = f3max - 2*nSpacing*eps, 2*nSpacing*eps
                #else:
                #    data[i][idx1], data[i][idx2] = injf3 - nSpacing*eps, 2*nSpacing*eps
                data[i][idx1], data[i][idx2] = injf3 - nSpacing*eps, 2*nSpacing*eps

            # f4dot
            if freqDerivOrder >= 4:
                f4min, f4max, _ = fr.f4BroadRange(injf0, 0, injf1, injf1, injf2, injf2)
                idx1, idx2 = freqParamName[2], freqDerivParamName[2]
                eps = spacing[idx2]
                injf4 =  injData[i][idx1] # Weave use "Freq" for injection
                #if (injf4 - nSpacing*eps)<f4min:
                #    data[i][idx1], data[i][idx2] = f4min, 2*nSpacing*eps
                #elif (injf4 + nSpacing*eps)>f4max:
                #    data[i][idx1], data[i][idx2] = f4max - 2*nSpacing*eps, 2*nSpacing*eps
                #else:
                #    data[i][idx1], data[i][idx2] = injf4 - nSpacing*eps, 2*nSpacing*eps
                data[i][idx1], data[i][idx2] = injf4 - nSpacing*eps, 2*nSpacing*eps

        data = Table(data)
        data.add_column(self.target.alpha*np.ones(n), name='alpha')
        data.add_column(self.target.dalpha*np.ones(n), name='dalpha')
        data.add_column(self.target.delta*np.ones(n), name='delta')
        data.add_column(self.target.ddelta*np.ones(n), name='ddelta')
 
        return fits.BinTableHDU(data)
    
    def _genParam(self, h0, freq, nBands=None, nInj=1, nAmp=1, injFreqDerivOrder=4, skyUncertainty=0, freqDerivOrder=2, stage='search', cluster=False, workInLocalDir=False):
        if freqDerivOrder > 4:
            print('Error: frequency derivative order larger than 4.')
        if injFreqDerivOrder > 4:
            print('Error: Injection frequency derivative order larger than 4.')
            
        searchParamDict, injParamDict = {}, {}
       
        taskName = utils.taskName(self.target, stage, self.cohDay, freqDerivOrder, freq)
        dataFilePath = fp.outlierFilePath(self.target, freq, taskName, stage, cluster=cluster)
        if workInLocalDir:
            dataFilePath = Path(dataFilePath).name
        nonSatBands = fits.getdata(dataFilePath, 3)['nonSatBand']  

        ip = self.genInjParamTable(nonSatBands, h0, freq, nInj, nAmp, injFreqDerivOrder, skyUncertainty)
        sp = self.genSearchRangeTable(dataFilePath, freq, ip.data, stage, freqDerivOrder)
        
        injParamDict[str(freq)] = ip
        searchParamDict[str(freq)] = sp
       
        return searchParamDict, injParamDict        
