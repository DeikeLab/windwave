#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --time=04:00:00
#SBATCH --mail-type=end
#SBATCH --mail-user=jiarongw@princeton.edu

# This is a slurm script for Della
# The executable name
EXE=wavewind_LEVELmin6
#Parameter value
LEVEL=11
ak=0.05
BO=0.53
RE=5000. #Default 40000
m=5
B=0
UstarRATIO=1

export ScratchDir="/scratch/gpfs/jiarongw/test/linear_m${m}B${B}Ustar${UstarRATIO}ak${ak}Bo${BO}Re${RE}LEVEL${LEVEL}_min6"
echo $ScratchDir

#Copy the executable from executable directory
rm -rf $ScratchDir
mkdir -p $ScratchDir
#Copy the executable from executable directory
cp /scratch/gpfs/jiarongw/executable/f$EXE/$EXE $ScratchDir
cd $ScratchDir

srun ./$EXE $LEVEL $ak $BO $RE $m $B $UstarRATIO

#To move the whole directory to /tigress
#cp -r $ScratchDir /tigress/jiarongw/
#rm -rf $ScratchDir
