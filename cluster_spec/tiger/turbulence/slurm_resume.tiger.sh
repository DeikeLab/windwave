#!/bin/bash
#SBATCH --nodes=4
#SBATCH --ntasks=160
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=jiarongw@princeton.edu

#The executable name
EXE=curved_fixREtau
#Parameter value
LEVEL=10
RE_tau=1800 #Default 40000
BO=10
g=1
ak=0.1
TIME=200
emaxRATIO=0.3


export ScratchDir="/scratch/gpfs/jiarongw/turbulence/${EXE}_REtau${RE_tau}_BO${BO}_g${g}_ak${ak}_LEVEL${LEVEL}_emax${emaxRATIO}"
echo $ScratchDir
#rm -rf $ScratchDir
#mkdir -p $ScratchDir
cp /scratch/gpfs/jiarongw/executable/f$EXE/$EXE $ScratchDir 
cp /home/jiarongw/windwave/project_spec_source/turbulence/curved/${EXE}.c $ScratchDir
cd $ScratchDir
#cp ../curved_uniform_forcing_moving_RE10000_Ustar0.5_ak0.1_LEVEL9_04_refinewater_precursor/dump45 ./restart
cp ./dump203 ./restart
srun ./$EXE $RE_tau $BO $LEVEL $g $ak $TIME $emaxRATIO > message 2>&1 

#To move the whole directory to /tigress
#cp -r $ScratchDir /tigress/jiarongw/
#rm -rf $ScratchDir
