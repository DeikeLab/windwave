#!/bin/bash
#SBATCH --nodes=4
#SBATCH --ntasks=160
#SBATCH --cpus-per-task=1
#SBATCH --time=48:00:00
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=jiarongw@princeton.edu

#The executable name
EXE=curved_uniform_forcing_moving
#Parameter value
LEVEL=10
RE=10000 #Default 40000
BO=10
UstarRATIO=0.5
ak=0.1
TIME=50
emaxRATIO=0.5
#START=0
#END=50


export ScratchDir="/scratch/gpfs/jiarongw/turbulence/${EXE}_RE${RE}_Ustar${UstarRATIO}_ak${ak}_LEVEL${LEVEL}_05_refinewater_precursor"
echo $ScratchDir
#rm -rf $ScratchDir
#mkdir -p $ScratchDir
#cp /scratch/gpfs/jiarongw/executable/f$EXE/$EXE $ScratchDir 
#cp /home/jiarongw/windwave/project_spec_source/turbulence/curved/${EXE}.c $ScratchDir
cd $ScratchDir
cp ./dump43.6 ./restart
srun ./$EXE $RE $BO $LEVEL $UstarRATIO $ak $TIME $emaxRATIO 

#To move the whole directory to /tigress
#cp -r $ScratchDir /tigress/jiarongw/
#rm -rf $ScratchDir
