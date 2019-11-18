#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=08:00:00
#SBATCH --mail-type=end
#SBATCH --mail-user=jiarongw@princeton.edu


LEVEL=11
ak=0.05
BO=0.53
RE=5000. #40000 default
m=5
B=0
UstarRATIO=1.2

export ScratchDir="/scratch/network/jiarongw/parameter/linear_m${m}B${B}Ustar${UstarRATIO}ak${ak}Bo${BO}Re${RE}LEVEL${LEVEL}"
echo $ScratchDir

#SBATCH --output=${ScratchDir}/slurm.out
cd $ScratchDir
cp ./dump1.34375 ./restart
srun ./wavewind_parallel $LEVEL $ak $BO $RE $m $B $UstarRATIO >& output.log
