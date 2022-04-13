#!/bin/bash

#SBATCH -J turbulence_nograph          # Job name
#SBATCH -o nograph.o%j       # Name of stdout output file
#SBATCH -e nograph.e%j       # Name of stderr error file
#SBATCH -p long          # Queue (partition) name
#SBATCH -N 8            # Total # of nodes 
#SBATCH -n 512             # Total # of mpi tasks
#SBATCH -t 120:00:00        # Run time (hh:mm:ss)
#SBATCH -A TG-OCE140023    # Allocation name


#The executable name
EXE=curved_fixREtau_boundary
#Parameter value
LEVEL=10
RE_tau=720 #Default 40000
BO=200
g=4
ak=0
TIME=50
emaxRATIO=0.3
alterMU=8

export ScratchDir="/scratch/06342/tg856416/turbulence/${EXE}_REtau${RE_tau}_BO${BO}_g${g}_ak${ak}_MU${alterMU}_LEVEL${LEVEL}_emax${emaxRATIO}"
echo $ScratchDir
#rm -rf $ScratchDir
#mkdir -p $ScratchDir
cp /work/06342/tg856416/stampede2/f$EXE/$EXE $ScratchDir 
cp /home1/06342/tg856416/windwave/project_spec_source/turbulence/curved/${EXE}.c $ScratchDir
cd $ScratchDir
cp ./dump172 ./restart
mkdir ./field
mkdir ./eta
ibrun ./$EXE $RE_tau $BO $LEVEL $g $ak $TIME $emaxRATIO $alterMU > message 2>&1 

# ---------------------------------------------------
