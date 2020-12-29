#!/bin/bash
#SBATCH --nodes=2
#SBATCH --ntasks=80
#SBATCH --cpus-per-task=1
#SBATCH --time=04:00:00
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=jiarongw@princeton.edu

#The executable name
EXE=curved_uniform_forcing_moving
#Parameter value
LEVEL1=8
LEVEL2=7
RE=10000 #Default 40000
UstarRATIO=0.5
ak=0.1
TIME=50
uemaxRATIO=0.3

#START=0
#END=50


export ScratchDir="/scratch/gpfs/jiarongw/turbulence/${EXE}_RE${RE}_Ustar${UstarRATIO}_ak${ak}_LEVEL${LEVEL1}_03_notmoving"
echo $ScratchDir
rm -rf $ScratchDir
mkdir -p $ScratchDir
cp /scratch/gpfs/jiarongw/executable/f$EXE/$EXE $ScratchDir
cp /home/jiarongw/windwave/project_spec_source/turbulence/curved/${EXE}.c $ScratchDir
cd $ScratchDir
#cp ../curved_uniform_forcing_moving_RE10000_Ustar0.5_ak0.1_LEVEL9_02/dump4.6 ./restart
srun ./$EXE $RE $LEVEL1 $LEVEL2 $UstarRATIO $ak $TIME $uemaxRATIO > message 2>&1

#To move the whole directory to /tigress
#cp -r $ScratchDir /tigress/jiarongw/
#rm -rf $ScratchDir
