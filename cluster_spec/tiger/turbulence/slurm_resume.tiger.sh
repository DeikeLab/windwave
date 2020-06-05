#!/bin/bash
#SBATCH --nodes=16
#SBATCH --ntasks-per-node=40
#SBATCH --time=08:00:00
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=jiarongw@princeton.edu

#The executable name
EXE=twophase_uniform_forcing_weaker
#Parameter value
LEVEL1=8
LEVEL2=7
RE=40000 #Default 40000
#START=0
#END=50


export ScratchDir="/scratch/gpfs/jiarongw/turbulence/twophase_uniform_forcing_weaker"
echo $ScratchDir
#Copy the executable from executable directory
cp /scratch/gpfs/jiarongw/executable/f$EXE/$EXE $ScratchDir
cd $ScratchDir
srun ./$EXE $RE $LEVEL1 $LEVEL2  

#To move the whole directory to /tigress
#cp -r $ScratchDir /tigress/jiarongw/
#rm -rf $ScratchDir
