user = 'hoitim.cheung'
homeDir = '/home/hoitim.cheung/gc/'
OSDFDir = '/osdf/igwn/cit/staging/hoitim.cheung/'
startTime = 1368970000
fBand = 0.1 #Resolution of search in frequency
semiMM = 0.2 #The maximum accepted mismatch of the semicoherent metric
cohMM = 0.1 #The maximum accepted mismatch of the N coherent metrics (N=number of segments)
secondsInDay = 86400

sftSource = 'C00_Gated_G02_1800s'
accGroup = 'ligo.prod.o4.cw.directedisolated.semicoherent' #argument for submitting jobs in LIGO cluster

nc_min, nc_max = 2.0, 7.0
injOutlier_nSpacing = 4.0 # for injection and followUp 4
cluster_nSpacing = 3.0 # for injection and followUp 4
followUp_nSpacing = 5.0 # for injection and followUp 
