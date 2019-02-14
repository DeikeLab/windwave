#!/bin/bash
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=8
#SBATCH --time=2:00:00
#SBATCH --mail-user=jiarongw@princeton.edu

LEVEL=11
ak=0.05
BO=200
RE=40000. #40000 default
m=10
B=0
UstarRATIO=1

export ScratchDir="/scratch/network/jiarongw/parameter/m${m}B${B}Ustar${UstarRATIO}ak${ak}"
echo $ScratchDir
rm -rf $ScratchDir
mkdir -p $ScratchDir
cp ./wavewind_parallel $ScratchDir

#SBATCH --output=${ScratchDir}/slurm.out
cd $ScratchDir
#mv ./dump ./restart
srun ./wavewind_parallel $LEVEL $ak $BO $RE $m $B $UstarRATIO

