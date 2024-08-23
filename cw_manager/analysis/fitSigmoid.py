import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import warnings

class fitSigmoid:    
    def __init__(self, target, nInj=1000, nAmp=8):
        self.target = target 
        self.injPerPoint = int(nInj/nAmp)         # number of injections per point

    def sigmoid(self, x, k, x0):
        return 1./(1.+np.exp(-k*(x-x0))) 

    def inv_sigmoid(self, y, k, x0):
        return -np.log(1./y-1.) /k +x0        

    def rescale_h0(self, h0, h0_list):
        self.h0_mean = h0_list.mean()
        self.h0_max = h0_list.max()
        return (h0-self.h0_mean)/self.h0_max

    def _invRescale_h0(self, x):
        return x*self.h0_max + self.h0_mean

    def binomialError(self, y, n):
        err =  np.sqrt(y*(1.-y)/n)
        err[err==0] = 1./n
        return err

    def h0_fromPercentile(self, per=0.95):    
        x = self.inv_sigmoid(per, *self.popt)
        a, b = self.popt.T
        dxda =  -1.0/a**2 * np.log(per/(1-per))
        dxdb = 1.0
        dx = np.sqrt(self.pcov[0,0]*dxda**2 + self.pcov[1,1]*dxdb**2 + 2*self.pcov[0,1]*dxda*dxdb)
        dh = dx*self.h0_max
        h = self._invRescale_h0(x)
        return h, dh/h
    def fit(self, freq, h0_list, p):
        if p[0] > 0.8 or p[-1] <0.95:
            warnings.warn("Freq={0} is problematic.".format(freq))
            
        err = self.binomialError(p, self.injPerPoint)
        x = self.rescale_h0(h0_list, h0_list)
        self.popt, self.pcov = curve_fit(self.sigmoid, x, p, p0=[5,0], sigma=err) 
        #perr = np.sqrt(np.diag(pcov))    
        return 0
    
    def plot(self, h0_list, p, savePath=None):
        fig, ax = plt.subplots()
        colors = plt.rcParams['axes.prop_cycle'].by_key()['color'] 
        
        h0_low, _ = self.h0_fromPercentile(per=0.01)
        h0_high, _ = self.h0_fromPercentile(per=0.99)
        h0_arr = np.linspace(h0_low, h0_high, 100)
        
        x = self.rescale_h0(h0_arr, h0_list)
        p_inter = self.sigmoid(x, *self.popt)
        h0, err = self.h0_fromPercentile(per=p_inter)
        ax.plot(h0_arr, p_inter, color=colors[0], label='fit')
        ax.errorbar(h0_arr, p_inter, xerr=h0*err, fmt='', color=colors[0], ls='none')
                 
        ax.scatter(h0_list, p, color='black')
        err = self.binomialError(p, self.injPerPoint)
        ax.errorbar(h0_list, p, yerr=err, fmt='', label='data', color='k', ls='none')
               
        h95, dx = self.h0_fromPercentile(per=0.95)
        p95 = 0.95
        x_highlighting = [0, h95, h95]
        y_highlighting = [p95, p95, -1]

        plt.plot(x_highlighting, y_highlighting, color='r', label=r'$h_{95}=$'+'{:.3e}'.format(h95)+r'$\pm$'+'{:.2f}%'.format(dx*100))

        ax.legend()
        ax.set_ylim(-0.05, 1.05)
        ax.set_xlim(0.9*min(h0_list[0], h0_arr[0]), 1.05*max(h0_list[-1], h0_arr[-1]))       
        ax.set_xlabel(r'$h_0$')
        ax.set_ylabel(r'$p_{\mathrm{det}}$')
        
        if savePath is not None:
            fig.savefig(savePath)
        return fig
