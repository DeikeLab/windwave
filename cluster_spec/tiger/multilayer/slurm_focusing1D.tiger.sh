#!/bin/bash
#SBATCH --nodes=2
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
#SBATCH --time=00:30:00
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=jiarongw@princeton.edu

#The executable name
EXE=field_focusing_1D
#Parameter value
NLAYER=30
LEVEL=8
TEND=60
nu=0.000025
#Initial field info
FIELD=P0010


export ScratchDir="/scratch/gpfs/jiarongw/multilayer/${EXE}_${FIELD}_RE40000_256_30"
echo $ScratchDir
rm -rf $ScratchDir
mkdir -p $ScratchDir
#Copy the executable from executable directory
cp /scratch/gpfs/jiarongw/executable/f$EXE/$EXE $ScratchDir
cp ${SCRATCH}/multilayer/F_kx_${FIELD} $ScratchDir/F_kx 
#cp -r ${SCRATCH}/multilayer/${EXE}_${FIELD}/dump59 $ScratchDir/restart
cd $ScratchDir
mkdir ./surface
echo srun ./$EXE $NLAYER $LEVEL $TEND $nu
srun ./$EXE $NLAYER $LEVEL $TEND $nu

#To move the whole directory to /tigress
#cp -r $ScratchDir /tigress/jiarongw/
#rm -rf $ScratchDir
