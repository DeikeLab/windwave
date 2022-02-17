#!/bin/bash
#SBATCH --nodes=6
#SBATCH --ntasks=216
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=jiarongw@princeton.edu

#The executable name
EXE=onephase
#Parameter value

export ScratchDir="/scratch/gpfs/jiarongw/turbulence/channel/${EXE}_180_LEVEL10"
echo $ScratchDir
#rm -rf $ScratchDir
mkdir -p $ScratchDir
cp /scratch/gpfs/jiarongw/executable/f$EXE/$EXE $ScratchDir
cp /home/jiarongw/windwave/project_spec_source/turbulence/channel/${EXE}.c $ScratchDir
cd $ScratchDir
#cp ../curved_uniform_forcing_moving_RE10000_Ustar0.5_ak0.1_LEVEL9_02/dump4.6 ./restart
mkdir field
srun ./$EXE 10 180  > message 2>&1

#To move the whole directory to /tigress
#cp -r $ScratchDir /tigress/jiarongw/
#rm -rf $ScratchDir
