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

# parameters for the calculation are read in
source $PARAMS

# old folders are deleted, new ones are created, some "statistics" are written
./scripts/init.sh $PARAMS

# all programs are compiled
make >> ${OUTDIR}/log.txt 2>&1

if [ ! -d $OUTDIR ]
then
  mkdir $OUTDIR
fi
cd ${OUTDIR} >> ${OUTDIR}/log.txt 2>&1




# Run the conversion to bond-angle-torsion coordinates (program BAT_builder.x), specifying the .top input file, the .xtc input file, the name of the output file, the name of the backbone atoms (as in the .top file)
../exec/BAT_builder.x -t ../${TOP} -x ../${TRJ} -o ${NAME}.bat -bb "${BACKBONE_ATOMS}" >> log.txt 2>&1

# Run the calculation of the entropy terms (program PARENT.x) specifying the input .bat file, the output .par file and the number of bins for bonds, angles and dihedrals (torsions)  in 1D and 2D for building the histograms
# Depending on the architecture of your cluster/workstation as well as your Open MPI version, you might want to change the mpirun parameters. Especially if you are using InfiniBand you might want to change "tcp" to "openib"
mpirun --report-bindings --map-by ppr:1:socket:pe=$OMP_NUM_THREADS --use-hwthread-cpus --bind-to hwthread --mca btl vader,self,tcp ../exec/PARENT.x ${NAME}.bat ${NAME}.par ${BBINS1D} ${ABINS1D} ${DBINS1D} ${BBINS2D} ${ABINS2D} ${DBINS2D} >> log.txt 2>&1

# Extract the bond-angle-torsion topologyby using get_PAR_info.x on the .par file and redirecting stdout to the textfile ${NAME}_MIE_topology.txt (also redirect stderr to log.txt)
../exec/get_PAR_info.x ${NAME}.par > ${NAME}_MIE_topology.txt 2>> log.txt
# Extract the entropy and mutual information terms by using get_PAR_MIE.x on the .par file and redirecting stdout to the textfile ${NAME}_MIE.txt (also redirect stderr to log.txt)
../exec/get_PAR_MIE.x ${NAME}.par > ${NAME}_MIE.txt 2>> log.txt

# Run the calculation of the MIST approximation specifying the input .par file and the output .txt file
# Depending on the architecture of your cluster/workstation as well as your Open MPI version, you might want to change the mpirun parameters. Especially if you are using InfiniBand you might want to change "tcp" to "openib"
mpirun --report-bindings --map-by ppr:1:socket:pe=$OMP_NUM_THREADS --use-hwthread-cpus --bind-to hwthread --mca btl vader,self,tcp ../exec/get_PAR_MIST.x ${NAME}.par ${NAME}_MIST.txt >> log.txt 2>&1

#write the date of finishing
date >> log.txt 2>&1




































