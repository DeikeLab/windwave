#!/bin/bash
#PBS -N mup5L10
#PBS -A UPRI0020
#PBS -l walltime=12:00:00
#PBS -q regular
#PBS -m abe
#PBS -M jiarongw@princeton.edu
#PBS -l select=5:ncpus=36:mpiprocs=36

#The executable name
EXE=curved_fixREtau_boundary
#Parameter value
LEVEL=10
RE_tau=720 #Default 40000
BO=200
g=1
ak=0.1
TIME=250
emaxRATIO=0.3

export ScratchDir="/glade/scratch/jiarongw/turbulence/${EXE}_REtau${RE_tau}_BO${BO}_g${g}_ak${ak}_LEVEL${LEVEL}_emax${emaxRATIO}"
echo $ScratchDir
rm -rf $ScratchDir
mkdir -p $ScratchDir
cp $SCRATCH/executable/f$EXE/$EXE $ScratchDir 
cp $SCRIPT/turbulence/curved/${EXE}.c $ScratchDir
cd $ScratchDir
cp ../restart ./
mkdir ./eta
mkdir ./field
#cp $SCRATCH/executable/fdiagnosis/diagnosis ./
#module load python
#python call.py
mpirun -np 180 ./$EXE $RE_tau $BO $LEVEL $g $ak $TIME $emaxRATIO > message 2>&1 
