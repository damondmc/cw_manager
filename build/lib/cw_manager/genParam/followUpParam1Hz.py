from ..utils import utils as utils
from tqdm import tqdm
import numpy as np
from astropy.io import fits
from astropy.table import Table
from . import frequencyRange as fr
from ..utils import setup_parameter as setup
from ..utils import filePath as fp
from pathlib import Path
                
class followUpParams():    
    def __init__(self, target, obsDay, fBand=0.1):
        self.obsDay = obsDay
        self.fBand = fBand
        self.target = target
    
    def makeFollowUpTable(self, cohDay, freq, stage, oldFreqDerivOrder, newFreqDerivOrder, cluster=False, workInLocalDir=False): 
        nSpacing = setup.followUp_nSpacing
        taskName = utils.taskName(self.target, stage, cohDay, oldFreqDerivOrder, freq)
        dataFilePath = fp.outlierFilePath(self.target, freq, taskName, stage, cluster=cluster)
        if workInLocalDir:
            dataFilePath = Path(dataFilePath).name
        
        data = Table.read(dataFilePath, hdu=1)  

        # check if there's outlier 
        if len(data) == 0:
            return fits.BinTableHDU(data=data)        
        
        freqParamName, freqDerivParamName = utils.phaseParamName(oldFreqDerivOrder)
        newFreqParamName, newFreqDerivParamName = utils.phaseParamName(newFreqDerivOrder)
        for _f, _df in zip(freqParamName, freqDerivParamName):
            data[_f] = data[_f] - nSpacing*data[_df] # it may change the sign of the parameter (from -ve to +ve or vice versa), due with it in the next step 
            data[_df] = 2*nSpacing*data[_df] 
            newFreqParamName.remove(_f)
            newFreqDerivParamName.remove(_df)
            
        # make the search range to be physical/consistent to the signal model (e.g. f1dot<=0, f2dot >=0)
        
        # f0
        f0min, f0max, _ = fr.f0BroadRange(data['freq'], data['df'])
        idx1, idx2 = freqParamName[0], freqDerivParamName[0]
        eps = data[idx2]
       
        # if the lower nominal freq search limit is lower than the broad range lower limit, relocate it to broad range lower limit. 
        mask = data[idx1]<f0min
        data[idx1][mask] = f0min[mask]
        mask = (data[idx1]+eps)>f0max
        data[idx1][mask] = f0max[mask] - eps[mask]
        
        # f1dot
        f1min, f1max, _ = fr.f1BroadRange(data['freq'], data['df'], self.target.tau)
        idx1, idx2 = freqParamName[1], freqDerivParamName[1]
        eps = data[idx2] 
        # if the lower f1dot search limit is lower than the broad range lower limit, relocate it to broad range lower limit. 
        mask = (data[idx1]+eps)>f1max
        data[idx1][mask] = f1max[mask] - eps[mask]
    
        # f2dot
        f2min, f2max, _ = fr.f2BroadRange(data['freq'], data['df'], data['f1dot'], data['f1dot'] + data['df1dot'])
        idx1, idx2 = freqParamName[2], freqDerivParamName[2]
        eps = data[idx2]
        mask = (data[idx1]+eps)>f2max
        data[idx1][mask] = f2max[mask] - eps[mask]
        #relocate mininum at the end to make sure f2dot >=0
        mask = data[idx1]<f2min
        data[idx1][mask] = f2min[mask]

        ## f3dot
        if oldFreqDerivOrder >=3:
            #print('f3')
            f3min, f3max, f3band = fr.f3BroadRange(data['freq'], data['df'], data['f1dot'], data['f1dot'] + data['df1dot'], data['f2dot'], data['f2dot']+data['df2dot'])
            idx1, idx2 = freqParamName[3], freqDerivParamName[3]
            data[idx1] = f3min 
            data[idx2] = f3band
    
        # f4dot
        if oldFreqDerivOrder >=4:
            #print('f4')
            f4min, f4max, _ = fr.f4BroadRange(data['freq'], data['df'], data['f1dot'], data['f1dot'] + data['df1dot'], data['f2dot'], data['f2dot']+data['df2dot'])
            idx1, idx2 = freqParamName[4], freqDerivParamName[4]
            eps = data[idx2]
            mask = (data[idx1]+eps)>f4max
            data[idx1][mask] = f4max[mask] - eps[mask]
            #relocate mininum at the end to make sure f4dot >=0
            mask = data[idx1]<f4min
            data[idx1][mask] = f4min[mask]

        # new added f3dot - f4dot (didnt search over in the previous stage)
        for _f, _df in zip(newFreqParamName, newFreqDerivParamName):
            if _f == 'f3dot':
                #print('new f3')
                f3min, f3max, f3band = fr.f3BroadRange(
                    data['freq'], data['df'], data['f1dot'], data['f1dot']+data['df1dot'], data['f2dot'], data['f2dot']+data['df2dot']
                )
                data.add_column(f3min, name=_f)
                data.add_column(f3band, name=_df)
                
            if _f == 'f4dot':
                #print('new f4')
                f4min, f4max, f4band = fr.f4BroadRange(
                    data['freq'], data['df'], data['f1dot'], data['f1dot']+data['df1dot'], data['f2dot'], data['f2dot']+data['df2dot']
                )
                data.add_column(f4min, name=_f)
                data.add_column(f4band, name=_df)
                
        return fits.BinTableHDU(data=data)
      
    def makeInjectionTable(self, cohDay, freq, stage, oldFreqDerivOrder, cluster=False):  
        taskName = utils.taskName(self.target, stage, cohDay, oldFreqDerivOrder, freq)
        dataFilePath = fp.outlierFilePath(self.target, freq, taskName, stage, cluster=cluster)
        injData = Table.read(dataFilePath, hdu=2)  
        
        return fits.BinTableHDU(data=injData)
    
    def genFollowUpParam(self, cohDay, freq, stage, oldFreqDerivOrder, newFreqDerivOrder, cluster=False, workInLocalDir=False):               
        if oldFreqDerivOrder > 4 or newFreqDerivOrder > 4:
            print('Error: frequency derivative order larger than 4.')
            
        params = {}
        dataset = self.makeFollowUpTable(cohDay, freq, stage, oldFreqDerivOrder, newFreqDerivOrder, cluster, workInLocalDir)
        if dataset.data.size == 0:
            print('{0} Hz has no outlier to follow up.'.format(freq))
            
        params[str(freq)] = dataset
        
        print('Done parameter generation of {0} from {1} for {2} Hz.'.format(self.target.name, stage, freq))
        return params
     
    def genFollowUpParamFromInjection1Hz(self, cohDay, freq, stage, oldFreqDerivOrder, newFreqDerivOrder, cluster=False):                 
        if oldFreqDerivOrder > 4 or newFreqDerivOrder > 4:
            print('Error: frequency derivative order larger than 4.')
            
            
        params, injParams = {}, {}
        
        params[str(freq)] = self.makeFollowUpTable(cohDay, freq, stage, oldFreqDerivOrder, newFreqDerivOrder, cluster)
        injParams[str(freq)] = self.makeInjectionTable(cohDay, freq, stage, oldFreqDerivOrder, cluster)
            
        print('Done generation of {0} follow-up parameters for {1} Hz.'.format(self.target.name, freq))
        return params, injParams
