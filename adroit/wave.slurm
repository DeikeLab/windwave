#!/bin/bash
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=8
#SBATCH --time=4:00:00
#SBATCH --mail-user=jiarongw@princeton.edu

LEVEL=11
ak=0.05
BO=200
RE=100000. #40000 default
m=5
B=0
UstarRATIO=0.8

export ScratchDir="/scratch/network/jiarongw/parameter/pressure_testRE_m${m}B${B}Ustar${UstarRATIO}ak${ak}Bo${BO}Re${RE}"
echo $ScratchDir
rm -rf $ScratchDir
mkdir -p $ScratchDir
cp ./wavewind.c $ScratchDir
cp ./wavewind_parallel $ScratchDir

#SBATCH --output=${ScratchDir}/slurm.out
cd $ScratchDir
#mv ./dump ./restart
srun ./wavewind_parallel $LEVEL $ak $BO $RE $m $B $UstarRATIO

