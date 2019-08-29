#!/bin/bash
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=40
#SBATCH --time=40:00:00
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=jiarongw@princeton.edu

LEVEL=11
ak=0.05
BO=3.45
RE=31000. #40000 default
m=5
B=0
UstarRATIO=0.5

export ScratchDir="/scratch/gpfs/jiarongw/parameter/linear_m${m}B${B}Ustar${UstarRATIO}ak${ak}Bo${BO}Re${RE}LEVEL${LEVEL}"

#SBATCH --output=${ScratchDir}/slurm.out
cd $ScratchDir
cp ./dump30 ./restart
srun ./wavewind_parallel $LEVEL $ak $BO $RE $m $B $UstarRATIO
