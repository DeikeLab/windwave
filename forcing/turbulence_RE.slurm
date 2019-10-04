#!/bin/bash
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=40
#SBATCH --time=20:00:00
#SBATCH --mail-type=end
#SBATCH --mail-user=jiarongw@princeton.edu

RE=132800
LEVEL=8
t_start=10
t_end=250

export ScratchDir="/scratch/gpfs/jiarongw/turbulence/RE_compare/RE${RE}LEVEL${LEVEL}start${t_start}end${t_end}"
echo $ScratchDir
rm -rf $ScratchDir
mkdir -p $ScratchDir
cp ./turbulence_RE.c $ScratchDir
cp ./turbulence_RE $ScratchDir

#SBATCH --output=${ScratchDir}/slurm.out
cd $ScratchDir
srun ./turbulence_RE ${RE} ${LEVEL} ${t_start} ${t_end}

