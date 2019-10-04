#!/bin/bash
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=40
#SBATCH --time=40:00:00
#SBATCH --mail-type=end
#SBATCH --mail-user=jiarongw@princeton.edu

RE=40000
t_start=10
t_end=250

export ScratchDir="/scratch/gpfs/jiarongw/turbulence/RE_compare/RE${RE}start${t_start}end${t_end}"
echo $ScratchDir
cp ./turbulence_RE.c $ScratchDir
cp ./turbulence_RE $ScratchDir

#SBATCH --output=${ScratchDir}/slurm.out
cd $ScratchDir
cp ./dump106 ./restart
srun ./turbulence_RE ${RE} ${t_start} ${t_end}
