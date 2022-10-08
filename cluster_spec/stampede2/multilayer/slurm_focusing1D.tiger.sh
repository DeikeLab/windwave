#!/bin/bash

#!/bin/bash

#SBATCH -J turbulence_nograph          # Job name
#SBATCH -o nograph.o%j       # Name of stdout output file
#SBATCH -e nograph.e%j       # Name of stderr error file
#SBATCH -p normal          # Queue (partition) name
#SBATCH -N 1            # Total # of nodes 
#SBATCH -n 64             # Total # of mpi tasks
#SBATCH -t 02:00:00        # Run time (hh:mm:ss)
#SBATCH -A TG-OCE140023    # Allocation name

#The executable name
EXE=field_focusing_1D
#Parameter value
NLAYER=30
LEVEL=8
TEND=60
nu=0.000025
#Initial field info
FIELD=P0010


export ScratchDir="/scratch/06342/tg856416/multilayer/${EXE}_${FIELD}_RE40000_256_30"
echo $ScratchDir
rm -rf $ScratchDir
mkdir -p $ScratchDir
#Copy the executable from executable directory
cp /work2/06342/tg856416/stampede2/f$EXE/$EXE $ScratchDir 
cp /home1/06342/tg856416/windwave/project_spec_source/multilayer/${EXE}.c $ScratchDir
#Copy the initial field
cp ${SCRATCH}/multilayer/F_kx_${FIELD} $ScratchDir/F_kx
cd $ScratchDir
mkdir ./surface
echo ibrun ./$EXE $NLAYER $LEVEL $TEND $nu
ibrun ./$EXE $NLAYER $LEVEL $TEND $nu

#To move the whole directory to /tigress
#cp -r $ScratchDir /tigress/jiarongw/
#rm -rf $ScratchDir
