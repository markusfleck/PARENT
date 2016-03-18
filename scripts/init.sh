#!/bin/bash

# check if the parameter file was correctly supplied
PARAMS=$1

if [ "$PARAMS" == "" ]
then
  echo "USAGE: ${0} parameterfile"
  exit
fi
if [ ! -e $PARAMS ]
then
  echo "ERROR: no parameterfile named ${PARAMS}"
  exit
fi

# parameters are read in to get $OUTDIR
source $PARAMS

# Remove  the output directory, then recreate it. Generally write everything to a logfile called log.txt
rm -r ${OUTDIR} >> log.txt 2>&1
mkdir ${OUTDIR} >> log.txt 2>&1
mv log.txt ${OUTDIR}


# Write the hostname, the MPI-version, the date and the working directory to the logfile
hostname >> ${OUTDIR}/log.txt 2>&1
#~ mpirun --version >> ${OUTDIR}/log.txt 2>&1
date >> ${OUTDIR}/log.txt 2>&1
echo Working directory is `pwd` >> ${OUTDIR}/log.txt 2>&1


