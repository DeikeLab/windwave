#!/bin/bash
#SBATCH --nodes=4
#SBATCH --ntasks=160
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=jiarongw@princeton.edu

#The executable name
EXE=curved_uniform_forcing
#Parameter value
LEVEL1=8
LEVEL2=7
RE=40000 #Default 40000
#START=0
#END=50


export ScratchDir="/scratch/gpfs/jiarongw/turbulence/${EXE}_RE${RE}"
echo $ScratchDir
rm -rf $ScratchDir
mkdir -p $ScratchDir
cp /scratch/gpfs/jiarongw/executable/f$EXE/$EXE $ScratchDir
cd $ScratchDir
srun ./$EXE $RE $LEVEL1 $LEVEL2 

#To move the whole directory to /tigress
#cp -r $ScratchDir /tigress/jiarongw/
#rm -rf $ScratchDir
