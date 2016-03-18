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

rm -r exec obj ${OUTDIR}
# rm run.pbs.e* run.pbs.o*
