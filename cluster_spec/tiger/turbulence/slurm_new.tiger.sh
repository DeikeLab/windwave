#!/bin/bash
#SBATCH --nodes=4
#SBATCH --ntasks=160
#SBATCH --cpus-per-task=1
#SBATCH --time=4:00:00
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=jiarongw@princeton.edu

#The executable name
EXE=curved_fixREtau_nomatching
#Parameter value
LEVEL=10
RE_tau=180 #Default 40000
BO=200
g=1
ak=0.1
TIME=250
emaxRATIO=0.3
alterMU=8

export ScratchDir="/scratch/gpfs/jiarongw/turbulence/${EXE}_REtau${RE_tau}_BO${BO}_g${g}_ak${ak}_MU${alterMU}_LEVEL${LEVEL}_emax${emaxRATIO}"
echo $ScratchDir
rm -rf $ScratchDir
mkdir -p $ScratchDir
cp /scratch/gpfs/jiarongw/executable/f$EXE/$EXE $ScratchDir 
cp /home/jiarongw/windwave/project_spec_source/turbulence/curved/${EXE}.c $ScratchDir
cd $ScratchDir
cp ../curved_fixREtau_precursor_REtau180_ak0.1_LEVEL9_emax0.3/dump245 ./restart
mkdir ./field
mkdir ./eta
srun ./$EXE $RE_tau $BO $LEVEL $g $ak $TIME $emaxRATIO $alterMU > message 2>&1 

#To move the whole directory to /tigress
#cp -r $ScratchDir /tigress/jiarongw/
#rm -rf $ScratchDir
