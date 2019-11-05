#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --time=02:00:00
#SBATCH --mail-type=end
#SBATCH --mail-user=jiarongw@princeton.edu

# This is a slurm script for Della
#The executable name
EXE=wavewind
#Parameter value
LEVEL=10
ak=0.05
BO=0.36
RE=31000. #Default 40000
m=8
B=0
UstarRATIO=0.528


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
