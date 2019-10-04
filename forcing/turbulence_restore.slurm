#!/bin/bash
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=40
#SBATCH --time=1:00:00
#SBATCH --mail-type=end
#SBATCH --mail-user=jiarongw@princeton.edu

RE=40000
LEVEL=7

export ScratchDir="/scratch/gpfs/jiarongw/turbulence/setfree_RE40000_test1"
echo $ScratchDir
rm -rf $ScratchDir
mkdir -p $ScratchDir
cp ./turbulence_setfree.c $ScratchDir
cp ./turbulence_setfree $ScratchDir

#SBATCH --output=${ScratchDir}/slurm.out
cd $ScratchDir
cp /scratch/gpfs/jiarongw/turbulence/stopforcing_restore_second/dump999 ./restart
srun ./turbulence_setfree 

