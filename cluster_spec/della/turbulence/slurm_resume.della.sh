#!/bin/bash
#SBATCH --nodes=16
#SBATCH --ntasks-per-node=20
#SBATCH --time=20:00:00
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=jiarongw@princeton.edu

#The executable name
EXE=turbulence_RE_limited
#Parameter value
LEVEL1=8
LEVEL2=6
RE=64000 #Default 40000
START=0
END=50


export ScratchDir="/scratch/gpfs/jiarongw/turbulence/adlim_RE${RE}LEVEL1${LEVEL1}LEVEL2${LEVEL2}START${START}END${END}"
echo $ScratchDir
#Copy the executable from executable directory
cp /scratch/gpfs/jiarongw/executable/f$EXE/$EXE $ScratchDir
cd $ScratchDir
cp ./dump20 ./restart
srun ./$EXE $RE $LEVEL1 $LEVEL2 $START $END 

#To move the whole directory to /tigress
#cp -r $ScratchDir /tigress/jiarongw/
#rm -rf $ScratchDir
