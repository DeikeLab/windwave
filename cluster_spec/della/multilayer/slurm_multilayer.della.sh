#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=8:00:00
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
TEND=200
#Initial field info
FIELD=P0015N10layer60

export ScratchDir="/scratch/gpfs/jiarongw/multilayer/${EXE}_${FIELD}"
echo $ScratchDir
rm -rf $ScratchDir
mkdir -p $ScratchDir
#Copy the executable from executable directory
cp /scratch/gpfs/jiarongw/executable/f$EXE/$EXE $ScratchDir
cp -r ./pre_${FIELD} $ScratchDir/pre 
cd $ScratchDir
mkdir ./surface
srun ./$EXE $NLAYER $LEVEL $MAXLEVEL $MINLEVEL $ETAE $TEND

#To move the whole directory to /tigress
#cp -r $ScratchDir /tigress/jiarongw/
#rm -rf $ScratchDir
