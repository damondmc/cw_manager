#!/cvmfs/software.igwn.org/conda/envs/igwn/bin/python
import importlib
from pathlib import Path
import argparse
import numpy as np
from analysis.resultManager import resultManager
import setup.setup_parameter as setup

parser = argparse.ArgumentParser(description='Analyze Weave search result.')
parser.add_argument('--target',type=str,help='Target to be analyzed',default='CassA')
parser.add_argument('--obsDay',type=float,help='Total observation time',default=240)  
parser.add_argument('--cohDay',type=float,help='Coherence time.',default=5)  
parser.add_argument('--stage',type=str,help='Stage of the pipeline.',default='search')  
parser.add_argument('--freq',type=int,help='Frequency band to be analyzed',default=20)
parser.add_argument('--freqDerivOrder',type=int,help='Highest order used in the phase evolution.',default=2)
parser.add_argument('--numTopList',type=int,help='Number of loudest templates returned by Weave',default=1000)
parser.add_argument('--dfdot',type=float,help='fdot bandwidth for the search',default=8e-10)
parser.add_argument('--cluster',type=int,help='To cluster the result or not.',default=0)
args = parser.parse_args()


obsDay = args.obsDay
cohDay = args.cohDay
stage = args.stage
freq = args.freq
freqDerivOrder = args.freqDerivOrder
numTopList = args.numTopList
dfdot = args.dfdot
cluster = bool(args.cluster)

target = importlib.import_from_string(args.target)

np.random.seed(0)
rm = resultManager(target, obsDay=obsDay, cohDay=cohDay, stage=stage, freqDerivOrder=freqDerivOrder, numTopList=numTopList)
rm.writeSearchResult(freq, freq+1,  dfdot=dfdot, numTopListLimit=numTopList, cluster=cluster)