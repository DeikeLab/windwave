#!/bin/bash
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=40
#SBATCH --time=4:00:00
#SBATCH --mail-type=end
#SBATCH --mail-user=jiarongw@princeton.edu

RE=40000
LEVEL=8
t_start=10
t_end=250

export ScratchDir="/scratch/gpfs/jiarongw/turbulence/RE_compare/RE${RE}LEVEL${LEVEL}start${t_start}end${t_end}_resumetest"
echo $ScratchDir
mkdir $ScratchDir
cp ./turbulence_RE.c $ScratchDir
cp ./turbulence_RE $ScratchDir

#SBATCH --output=${ScratchDir}/slurm.out
cd $ScratchDir
cp ../../stopforcing/dump116 ./restart
srun ./turbulence_RE ${RE} ${LEVEL} ${t_start} ${t_end}

