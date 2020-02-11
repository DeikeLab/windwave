#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=02:00:00
#SBATCH --mail-type=end
#SBATCH --mail-user=jiarongw@princeton.edu

EXE=wavewind_refineu

LEVEL=10
ak=0.05
BO=3.23
RE=31000. #40000 default
m=5
B=0.5
UstarRATIO=0.5

export ScratchDir="/scratch/network/jiarongw/miscellaneous/linear_refineu_m${m}B${B}Ustar${UstarRATIO}ak${ak}Bo${BO}Re${RE}LEVEL${LEVEL}"
echo $ScratchDir
rm -rf $ScratchDir
mkdir -p $ScratchDir
#Copy the executable from executable directory
cp /scratch/network/jiarongw/executable/f$EXE/$EXE $ScratchDir
cd $ScratchDir
srun ./$EXE $LEVEL $ak $BO $RE $m $B $UstarRATIO >& output.log




