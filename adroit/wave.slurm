#!/bin/bash
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=8
#SBATCH --time=2:00:00
#SBATCH --mail-user=jiarongw@princeton.edu

LEVEL=11
ak=0.01
BO=3.45
RE=31000. #40000 default
m=5
B=0
UstarRATIO=0.97

export ScratchDir="/scratch/network/jiarongw/parameter/modified_m${m}B${B}Ustar${UstarRATIO}ak${ak}Bo${BO}Re${RE}"
echo $ScratchDir
rm -rf $ScratchDir
mkdir -p $ScratchDir
cp ./wavewind_parallel $ScratchDir

#SBATCH --output=${ScratchDir}/slurm.out
cd $ScratchDir
#mv ./dump ./restart
srun ./wavewind_parallel $LEVEL $ak $BO $RE $m $B $UstarRATIO

