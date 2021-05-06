#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=8:00:00
#SBATCH --mail-type=end
#SBATCH --mail-user=jiarongw@princeton.edu

# This is a slurm script for Della
# The executable name
EXE=wavewind_rerun_test
#Parameter value
LEVEL=11
ak=0.05
BO=200
RE=10660. #Default 40000
m=5
B=0
UstarRATIO=1

export ScratchDir="/scratch/gpfs/jiarongw/rerun/linear_${EXE}_Ustar${UstarRATIO}ak${ak}Bo${BO}Re${RE}LEVEL${LEVEL}"
echo $ScratchDir
rm -rf $ScratchDir
mkdir -p $ScratchDir
#Copy the executable from executable directory
cp /scratch/gpfs/jiarongw/executable/f$EXE/$EXE $ScratchDir
cd $ScratchDir
mkdir ./pressure
mkdir ./shear
srun ./$EXE $LEVEL $ak $BO $RE $m $B $UstarRATIO 0 > message 2>&1

#To move the whole directory to /tigress
#cp -r $ScratchDir /tigress/jiarongw/
#rm -rf $ScratchDir
