#!/bin/bash
#SBATCH --nodes=4
#SBATCH --ntasks=160
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=jiarongw@princeton.edu

#The executable name
EXE=curved_uniform_forcing_moving
#Parameter value
LEVEL1=9
LEVEL2=8
RE=10000 #Default 40000
UstarRATIO=0.5
ak=0.1
TIME=5
emaxRATIO=0.3
#START=0
#END=50


export ScratchDir="/scratch/gpfs/jiarongw/turbulence/${EXE}_RE${RE}_Ustar${UstarRATIO}_ak${ak}_LEVEL${LEVEL1}_03_refinewater"
echo $ScratchDir
#rm -rf $ScratchDir
#mkdir -p $ScratchDir
cp /scratch/gpfs/jiarongw/executable/f$EXE/$EXE $ScratchDir 
#cp /home/jiarongw/windwave/project_spec_source/turbulence/curved/${EXE}.c $ScratchDir
cd $ScratchDir
#cp ./dump29 ./restart
srun ./$EXE $RE $LEVEL1 $LEVEL2 $UstarRATIO $ak $TIME $emaxRATIO 

#To move the whole directory to /tigress
#cp -r $ScratchDir /tigress/jiarongw/
#rm -rf $ScratchDir
