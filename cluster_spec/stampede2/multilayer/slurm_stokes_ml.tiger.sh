#!/bin/bash

#SBATCH -J turbulence_nograph          # Job name
#SBATCH -o nograph.o%j       # Name of stdout output file
#SBATCH -e nograph.e%j       # Name of stderr error file
#SBATCH -p normal          # Queue (partition) name
#SBATCH -N 1            # Total # of nodes 
#SBATCH -n 64             # Total # of mpi tasks
#SBATCH -t 2:00:00        # Run time (hh:mm:ss)
#SBATCH -A TG-OCE140023    # Allocation name

#The executable name
EXE=stokes_ml
#Parameter value
NLAYER=30
N=256
coeff=0
RE=40000
ak=0.4

export ScratchDir="/scratch/06342/tg856416/multilayer/${EXE}_thetaH0.8_RE${RE}_nl${NLAYER}_N${N}_coeff${coeff}_ak${ak}"
echo $ScratchDir
rm -rf $ScratchDir
mkdir -p $ScratchDir
#Copy the executable from executable directory
cp /work2/06342/tg856416/stampede2/f$EXE/$EXE $ScratchDir 
cp /home1/06342/tg856416/windwave/project_spec_source/multilayer/${EXE}.c $ScratchDir
cd $ScratchDir
cp ../${EXE}_RE${RE}_nl${NLAYER}_N${N}_coeff${coeff}_ak${ak}/dump1 ./restart
echo ibrun ./$EXE $RE $NLAYER $coeff $N $ak
ibrun ./$EXE $RE $NLAYER $coeff $N $ak

#To move the whole directory to /tigress
#cp -r $ScratchDir /tigress/jiarongw/
#rm -rf $ScratchDir
