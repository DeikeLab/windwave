#!/bin/bash
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=40
#SBATCH --time=20:00:00
#SBATCH --mail-type=end
#SBATCH --mail-user=jiarongw@princeton.edu

export ScratchDir="/scratch/gpfs/jiarongw/turbulence/test2_stopforcing_bc2"
echo $ScratchDir
mkdir -p $ScratchDir
cp ./turbulence_stopforcing_bc2.c $ScratchDir
cp ./turbulence_stopforcing_bc2 $ScratchDir

#SBATCH --output=${ScratchDir}/slurm.out
cd $ScratchDir
cp ../test2_restore/dump115 ./restart
srun ./turbulence_stopforcing_bc2

