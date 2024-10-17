user = 'hoitim.cheung'
homeDir = '/home/hoitim.cheung/testing/'
OSDFDir = '/osdf/igwn/cit/staging/hoitim.cheung/'
#startTime = 1238166483 #GPS time of first SFT in the data
startTime = 1368970000
ObsDay = 240
fBand = 0.1 #Resolution of search in frequency
semiMM = 0.2 #The maximum accepted mismatch of the semicoherent metric
cohMM = 0.1 #The maximum accepted mismatch of the N coherent metrics (N=number of segments)
OSG = True # argument: to use OSG data in CVMFS format (if True) or data in hdfs format

secondsInDay = 86400
#freqDerivOrder = 2

sftSource = 'C00_Gated_G01_Sub60Hz_1800s'
accGroup = 'ligo.prod.o4.cw.directedisolated.semicoherent' #argument for submitting jobs in LIGO cluster

nc_min, nc_max = 2.0, 7.0
injOutlier_nSpacing = 4.0 # for injection and followUp 1-3
followUp_nSpacing = 5.0 # for injection and followUp 1-3
#injOutlier_nSpacing = 3.0 # for injection and followUp 1-3
#followUp_nSpacing = 3.0 # for injection and followUp 1-3
#injOutlier_nSpacing = 3.0 # for injection and followUp 4
#followUp_nSpacing = 3.0 # for injection and followUp 4

