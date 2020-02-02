#!/bin/bash
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=16
#SBATCH --time=2:00:00
#SBATCH --mail-type=end
#SBATCH --mail-user=jiarongw@princeton.edu

# This is a slurm script for Della
#The executable name
EXE=wavewind_linear_adaptive
#Parameter value
LEVEL=10
ak=0.05
BO=3.34
RE=31000. #Default 40000
m=5
B=0
UstarRATIO=0.44

export ScratchDir="/scratch/gpfs/jiarongw/miscellaneous/linear_adaptive_varyslope_refineu001_m${m}B${B}Ustar${UstarRATIO}ak${ak}Bo${BO}Re${RE}LEVEL${LEVEL}"
echo $ScratchDir
rm -rf $ScratchDir
mkdir -p $ScratchDir
#Copy the executable from executable directory
cp /scratch/gpfs/jiarongw/executable/f$EXE/$EXE $ScratchDir
cd $ScratchDir
mkdir ./pressure
mkdir ./shear
srun ./$EXE $LEVEL $ak $BO $RE $m $B $UstarRATIO

#To move the whole directory to /tigress
#cp -r $ScratchDir /tigress/jiarongw/
#rm -rf $ScratchDir
