#!/bin/bash
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=40
#SBATCH --time=20:00:00
#SBATCH --mail-type=end
#SBATCH --mail-user=jiarongw@princeton.edu


LEVEL=11
ak=0.05
BO=3.45
RE=31000. #40000 default
m=5
B=0
UstarRATIO=0.5
k=2k0

export ScratchDir="/scratch/gpfs/jiarongw/parameter/varying_k/m${m}B${B}Ustar${UstarRATIO}ak${ak}Bo${BO}Re${RE}LEVEL${LEVEL}k${k}"
echo $ScratchDir
rm -rf $ScratchDir
mkdir -p $ScratchDir
cp ./wavewind_parallel $ScratchDir

#SBATCH --output=${ScratchDir}/slurm.out
cd $ScratchDir
srun ./wavewind_parallel $LEVEL $ak $BO $RE $m $B $UstarRATIO

