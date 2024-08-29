#!/bin/bash
STARTTIME=$1
TCOH=$2
NSEG=$3
NSPINS=$4
SAVEPATH=$5
FILENAME=''${SAVEPATH}'/metricSetup/Start'${STARTTIME}'_TCoh'${TCOH}'_N'${NSEG}'_Spin'${NSPINS}'.fts'
echo ${FILENAME}

/cvmfs/software.igwn.org/conda/envs/igwn-py39-20231212/bin/lalpulsar_WeaveSetup --output-file=${FILENAME} --detectors='H1,L1' --first-segment=${STARTTIME}/${TCOH} --segment-count=${NSEG} --spindowns=${NSPINS}

