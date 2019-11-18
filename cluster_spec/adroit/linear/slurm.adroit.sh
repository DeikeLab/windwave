#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=04:00:00
#SBATCH --mail-type=end
#SBATCH --mail-user=jiarongw@princeton.edu

EXE=wavewind_watercurrent_0

LEVEL=11
ak=0.05
BO=0.53
RE=5000. #40000 default
m=5
B=0.5
UstarRATIO=1

export ScratchDir="/scratch/network/jiarongw/watercurrent/linear_m${m}B${B}Ustar${UstarRATIO}ak${ak}Bo${BO}Re${RE}LEVEL${LEVEL}"
echo $ScratchDir
rm -rf $ScratchDir
mkdir -p $ScratchDir
#Copy the executable from executable directory
cp /scratch/network/jiarongw/executable/f$EXE/$EXE $ScratchDir
cd $ScratchDir
srun ./$EXE $LEVEL $ak $BO $RE $m $B $UstarRATIO >& output.log




