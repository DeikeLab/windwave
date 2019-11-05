#!/bin/bash
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=40
#SBATCH --time=00:30:00
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=jiarongw@princeton.edu

#The executable name
EXE=wavewind
#Parameter value
LEVEL=10
ak=0.05
BO=3.45
RE=31000. #Default 40000
m=8
B=0
UstarRATIO=0.44


export ScratchDir="/scratch/gpfs/jiarongw/parameter/linear_m${m}B${B}Ustar${UstarRATIO}ak${ak}RE${RE}LEVEL${LEVEL}_4nodes160cores"
echo $ScratchDir
rm -rf $ScratchDir
mkdir -p $ScratchDir
#Copy the executable from executable directory
cp /scratch/gpfs/jiarongw/executable/f$EXE/$EXE $ScratchDir
cd $ScratchDir
srun ./$EXE $LEVEL $ak $BO $RE $m $B $UstarRATIO

#To move the whole directory to /tigress
#cp -r $ScratchDir /tigress/jiarongw/
#rm -rf $ScratchDir
