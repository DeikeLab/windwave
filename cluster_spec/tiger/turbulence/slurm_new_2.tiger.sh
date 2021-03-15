#!/bin/bash
#SBATCH --nodes=4
#SBATCH --ntasks=160
#SBATCH --cpus-per-task=1
#SBATCH --time=12:00:00
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=jiarongw@princeton.edu

#The executable name
EXE=curved_fixREtau
#Parameter value
LEVEL=9
RE_tau=180 #Default 40000
BO=200
g=4
ak=0.2
TIME=200
emaxRATIO=0.3


export ScratchDir="/scratch/gpfs/jiarongw/turbulence/${EXE}_REtau${RE_tau}_BO${BO}_g${g}_ak${ak}_LEVEL${LEVEL}_emax${emaxRATIO}"
echo $ScratchDir
rm -rf $ScratchDir
mkdir -p $ScratchDir
cp /scratch/gpfs/jiarongw/executable/f$EXE/$EXE $ScratchDir 
cp /home/jiarongw/windwave/project_spec_source/turbulence/curved/${EXE}.c $ScratchDir
cd $ScratchDir
mkdir ./eta
mkdir ./field
cp ../curved_fixREtau_precursor_REtau1800_ak0.2_LEVEL9_emax0.3/dump25 ./restart
srun ./$EXE $RE_tau $BO $LEVEL $g $ak $TIME $emaxRATIO > message 2>&1 

#To move the whole directory to /tigress
#cp -r $ScratchDir /tigress/jiarongw/
#rm -rf $ScratchDir
