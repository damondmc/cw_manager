#!/bin/bash


#NarrowBand
STARTTIME=1368970000
TCOH=864000
NSEG=24
NSPINS=3
FILENAME='/home/hoitim.cheung/snrsearch/o4/metricSetup/Start'${STARTTIME}'_TCoh'${TCOH}'_N'${NSEG}'_Spin'${NSPINS}'.fts'
echo ${FILENAME}

/cvmfs/software.igwn.org/conda/envs/igwn-py39-20231212/bin/lalpulsar_WeaveSetup --output-file=${FILENAME} --detectors='H1,L1' --first-segment=${STARTTIME}/${TCOH} --segment-count=${NSEG} --spindowns=${NSPINS}



#search stage
#TCOH=432000
#NSEG=36
#NSPINS=2

#followup
#TCOH=864000
#NSEG=18
#NSPINS=2
#FILENAME='./MetricSetup/Start1238166483_TCoh'${TCOH}'_N'${NSEG}'_Spin'${NSPINS}'.fts'
#echo ${FILENAME}

#lalpulsar_WeaveSetup --output-file=${FILENAME} --detectors='H1,L1' --first-segment=1238166483/${TCOH} --segment-count=${NSEG} --spindowns=${NSPINS}
