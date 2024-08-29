#Author - Damon Cheung
from ..utils import utils as utils
from . import frequencyRange as fr
from tqdm import tqdm
import numpy as np
from astropy.io import fits
from astropy.table import Table

class initSearchParams:    
    def __init__(self, fBand=0.1, freqDerivOrder=2):
        self.fBand = fBand
        freqParamName = ["freq", "f1dot", "f2dot", "f3dot", "f4dot"]
        freqDerivParamName = ["df", "df1dot", "df2dot", "df3dot", "df4dot"]
        self.freqParamName = freqParamName[:freqDerivOrder+1]
        self.freqDerivParamName = freqDerivParamName[:freqDerivOrder+1]
        
    def genParamTable(self, freq, nf1dots, nf2dots):
        n = int(nf1dots*nf2dots/self.fBand)
        data =np.recarray((n,), dtype=[(key, '>f8') for key in (self.freqParamName + self.freqDerivParamName)]) 
        
        for i in range(int(1.0/self.fBand)):
            f0 = freq + i *self.fBand
            f0min, f0max, f0band = fr.f0BroadRange(f0, self.fBand)
            for j in range(nf1dots):
                _f1min, _, _f1Band = fr.f1BroadRange(f0, self.fBand, self.target.tau)
                f1Band = _f1Band/nf1dots  # divide f1dot into n segment
                f1min = _f1min + j*f1Band
                f1max = f1min + f1Band
                if j == nf1dots - 1:
                    f1max = 0.0             # to mannually set f1dot upper limit to 0 (numerical accuracy/error exits)
                    f1Band = 0.0 - f1min
                        
                for k in range(nf2dots):
                    _f2min, _, _f2Band = fr.f2BroadRange(freq, self.fBand, f1min, f1max)
                    f2Band = _f2Band/nf2dots
                    f2min = _f2min + k*f2Band
                    f2max = f2min + f2Band
                    if k == 0:
                        f2min = 0.0
                        f2Band = f2max    # to mannually set f2dot lower limit to 0 (numerical accuracy/error exits)
                    
                    
                    idx = i * nf1dots * nf2dots + j * nf2dots + k 
                    data[idx]['freq'], data[idx]['df'] = f0min, f0band
                    data[idx]['f1dot'], data[idx]['df1dot'] = f1min, f1Band
                    data[idx]['f2dot'], data[idx]['df2dot'] = f2min, f2Band
                    
                    # sky location
#                    data[idx]['alpha'], data[idx]['dalpha'] = self.target.alpha, self.target.dalpha
#                    data[idx]['delta'], data[idx]['ddelta'] = self.target.delta, self.target.ddelta
        data = Table(data)
        data.add_column(self.target.alpha*np.ones(n), name='alpha')
        data.add_column(self.target.dalpha*np.ones(n), name='dalpha')
        data.add_column(self.target.delta*np.ones(n), name='delta')
        data.add_column(self.target.ddelta*np.ones(n), name='ddelta')           
        return fits.BinTableHDU(data)
        
    # need to do, at search result stage, append the table.
    def genParam(self, target, fmin, fmax, df1dot=1e-9, df2dot=1e-19):
        self.target = target
        params = {}
        for freq in tqdm(range(fmin, fmax)):
            nf1dots = fr.getNf1dot(freq, self.fBand, target.tau, df1dot=df1dot) # number of segment for f1dot range
            nf2dots = fr.getNf2dot(freq, self.fBand, target.tau, df2dot=df2dot) # number of segment for f2dot range
            params[str(freq)] = self.genParamTable(freq, nf1dots, nf2dots)
        return params        
