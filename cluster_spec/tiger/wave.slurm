#!/bin/bash
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=40
#SBATCH --time=0:30:00
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=jiarongw@princeton.edu

#The executable name
EXE=wavewind
#Parameter value
LEVEL=11
ak=0.05
BO=200
RE=4000 #Default 40000
m=8
B=0
UstarRATIO=1


export ScratchDir="/scratch/gpfs/jiarongw/parameter/m${m}B${B}Ustar${UstarRATIO}ak${ak}RE${RE}"
echo $ScratchDir
#rm -rf $ScratchDir
#mkdir -p $ScratchDir
#Copy the executable from executable directory
cp /scratch/gpfs/executable/f$EXE/EXE $ScratchDir
cd $ScratchDir
srun ./$EXE $LEVEL $ak $BO $RE $m $B $UstarRATIO

#To move the whole directory to /tigress
#cp -r $ScratchDir /tigress/jiarongw/
#rm -rf $ScratchDir
