#!/bin/bash
#PBS -N pre_ak01_L8
#PBS -A UPRI0020
#PBS -l walltime=12:00:00
#PBS -q regular
#PBS -m abe
#PBS -M jiarongw@princeton.edu
#PBS -l select=5:ncpus=36:mpiprocs=36

#The executable name
EXE=curved_fixREtau_boundary
#Parameter value
LEVEL=8
RE_tau=720 #Default 40000
ak=0.1
TIME=500
emaxRATIO=0.3

#Not actually relevant
BO=200
g=1
alterMU=1

export ScratchDir="/glade/scratch/jiarongw/turbulence/curved_fixREtau_precursor_REtau${RE_tau}_ak${ak}_LEVEL${LEVEL}_emax${emaxRATIO}"
echo $ScratchDir
mkdir -p $ScratchDir
cp $SCRATCH/executable/f$EXE/$EXE $ScratchDir 
cp $SCRIPT/turbulence/curved/${EXE}.c $ScratchDir
cd $ScratchDir
#cp ./dump81 ./restart
mkdir ./eta
mkdir ./field
#cp $SCRATCH/executable/fdiagnosis/diagnosis ./
#module load python
#python call.py
mpirun -np 180 ./$EXE $RE_tau $BO $LEVEL $g $ak $TIME $emaxRATIO $alterMU > message 2>&1
