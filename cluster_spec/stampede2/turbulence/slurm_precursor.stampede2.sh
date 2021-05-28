#!/bin/bash

#SBATCH -J turbulence_nograph          # Job name
#SBATCH -o nograph.o%j       # Name of stdout output file
#SBATCH -e nograph.e%j       # Name of stderr error file
#SBATCH -p normal          # Queue (partition) name
#SBATCH -N 8            # Total # of nodes 
#SBATCH -n 512             # Total # of mpi tasks
#SBATCH -t 24:00:00        # Run time (hh:mm:ss)
#SBATCH -A TG-OCE140023    # Allocation name


#The executable name
EXE=curved_fixREtau_precursor
#Parameter value
LEVEL=9
RE_tau=720 #Default 40000
BO=200
g=4
ak=0
TIME=500
emaxRATIO=0.3
alterMU=8


export ScratchDir="/scratch/06342/tg856416/turbulence/${EXE}_REtau${RE_tau}_ak${ak}_LEVEL${LEVEL}_emax${emaxRATIO}"
echo $ScratchDir
rm -rf $ScratchDir
mkdir -p $ScratchDir
cp /work/06342/tg856416/stampede2/f$EXE/$EXE $ScratchDir 
cp /home1/06342/tg856416/windwave/project_spec_source/turbulence/curved/${EXE}.c $ScratchDir
cd $ScratchDir
cp ../curved_fixREtau_precursor_precursor_REtau7200_ak0_LEVEL9_emax0.3/dump10 ./restart
mkdir ./field
mkdir ./eta
ibrun ./$EXE $RE_tau $BO $LEVEL $g $ak $TIME $emaxRATIO $alterMU > message 2>&1 

# ---------------------------------------------------
