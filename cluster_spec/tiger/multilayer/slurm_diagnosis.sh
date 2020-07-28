#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=42:00:00
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=jiarongw@princeton.edu

#The executable name
EXE=field
#Parameter value
NLAYER=60
LEVEL=8
MAXLEVEL=8
MINLEVEL=8
ETAE=0.01
TEND=60
#Initial field info
FIELD=focusing2D

export ScratchDir="/scratch/gpfs/jiarongw/multilayer/${FIELD}"
echo $ScratchDir
#rm -rf $ScratchDir
#mkdir -p $ScratchDir
#Copy the executable from executable directory
#cp /scratch/gpfs/jiarongw/executable/f$EXE/$EXE $ScratchDir
#cp -r ${SCRATCH}/multilayer/pre_${FIELD} $ScratchDir/pre 
cd $ScratchDir
#mkdir ./surface
echo srun ./$EXE $NLAYER $LEVEL $MAXLEVEL $MINLEVEL $ETAE $TEND
srun ./$EXE $NLAYER $LEVEL $MAXLEVEL $MINLEVEL $ETAE $TEND

module load python
python hello.py

#To move the whole directory to /tigress
#cp -r $ScratchDir /tigress/jiarongw/
#rm -rf $ScratchDir
