#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=4:00:00
#SBATCH --mail-user=jiarongw@princeton.edu

LEVEL=10
ak=0.1
BO=200
RE=100000. #40000 default
m=5
B=0
UstarRATIO=0.8

export ScratchDir="/scratch/network/jiarongw/parameter/pressure_serial_m${m}B${B}Ustar${UstarRATIO}ak${ak}Bo${BO}Re${RE}"
echo $ScratchDir
rm -rf $ScratchDir
mkdir -p $ScratchDir
cp ./wavewind.c $ScratchDir
cp ./wavewind_serial $ScratchDir

#SBATCH --output=${ScratchDir}/slurm.out
cd $ScratchDir
srun ./wavewind_serial $LEVEL $ak $BO $RE $m $B $UstarRATIO
