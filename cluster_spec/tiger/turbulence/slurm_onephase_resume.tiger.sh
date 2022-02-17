#!/bin/bash
#SBATCH --nodes=2
#SBATCH --ntasks=80
#SBATCH --cpus-per-task=1
#SBATCH --time=12:00:00
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=jiarongw@princeton.edu

#The executable name
EXE=square
#Parameter value
LEVEL=9
emax=0.3
RE=720

export ScratchDir="/scratch/gpfs/jiarongw/turbulence/channel/square720/${EXE}_LEVEL${LEVEL}_emax${emax}"
echo $ScratchDir
#rm -rf $ScratchDir
mkdir -p $ScratchDir
cp /scratch/gpfs/jiarongw/executable/f$EXE/$EXE $ScratchDir
cp /home/jiarongw/windwave/project_spec_source/turbulence/channel/${EXE}.c $ScratchDir
cd $ScratchDir
mkdir field
#cp ../dump180 ./restart
srun ./$EXE $LEVEL $RE $emax  > message 2>&1

#To move the whole directory to /tigress
#cp -r $ScratchDir /tigress/jiarongw/
#rm -rf $ScratchDir
