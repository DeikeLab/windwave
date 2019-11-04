#!/bin/bash
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=40
#SBATCH --time=10:00:00
#SBATCH --mail-type=end
#SBATCH --mail-user=jiarongw@princeton.edu

RE=4000

export ScratchDir="/scratch/gpfs/jiarongw/turbulence/RE_compare/RE${RE}_refine"
echo $ScratchDir
rm -rf $ScratchDir
mkdir -p $ScratchDir
mkdir -p $ScratchDir/movie
cp ./turbulence_refine_RE4000.c $ScratchDir
cp ./turbulence_refine_RE4000 $ScratchDir


#SBATCH --output=${ScratchDir}/slurm.out
cd $ScratchDir
cp ../RE4000start10end250/dump250 ./restart
srun ./turbulence_refine_RE4000 ${RE}
