#!/bin/bash
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=40
#SBATCH --time=0:30:00
#SBATCH --mail-type=end
#SBATCH --mail-user=jiarongw@princeton.edu


LEVEL=10
ak=0.05
BO=3.45
RE=31000. #40000 default
m=5
B=0
UstarRATIO=0.5

export ScratchDir="/scratch/gpfs/jiarongw/parameter/linear_nowind_fix_m${m}B${B}Ustar${UstarRATIO}ak${ak}Bo${BO}Re${RE}LEVEL${LEVEL}"
echo $ScratchDir
rm -rf $ScratchDir
mkdir -p $ScratchDir
cp ./wave_nowind_fix $ScratchDir

#SBATCH --output=${ScratchDir}/slurm.out
cd $ScratchDir
srun ./wave_nowind_fix $LEVEL $ak $BO $RE $m $B $UstarRATIO

