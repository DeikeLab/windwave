#!/bin/bash
#PBS -N mup5L10
#PBS -A UPRI0020
#PBS -l walltime=0:30:00
#PBS -q regular
#PBS -m abe
#PBS -M jiarongw@princeton.edu
#PBS -l select=2:ncpus=36:mpiprocs=36

#The executable name
EXE=curved_fixREtau
#Parameter value
LEVEL=10
g=1
BO=200
RE_tau=720
ak=0.25
TIME=500
emaxRATIO=0.3
alterMU=8

export ScratchDir="/glade/scratch/jiarongw/turbulence/${EXE}_precursor_REtau${RE_tau}_ak${ak}_LEVEL${LEVEL}_emax${emaxRATIO}"
echo $ScratchDir
cd $ScratchDir

# Output data on uniform grid command
mkdir ./field
mkdir ./eta
cp $SCRATCH/executable/fdiagnosis_run_interface/diagnosis_run_interface ./

OUTLEVEL=9
for i in `seq 0 1`
do
    let Snapshot=29+i
    echo $Snapshot
    mpirun -np 72 ./diagnosis_run_interface $RE_tau $BO $LEVEL $g $ak $TIME $emaxRATIO $alterMU $Snapshot $OUTLEVEL > message_diagnosis 2>&1 
done

#srun ./diagnosis_run $RE_tau $BO $LEVEL $g $ak $TIME $emaxRATIO $alterMU 75 $OUTLEVEL > message_diagnosis 2>&1 

#bash ./call.sh
#module load python
#python call.py

# Plotting command
#cp $SCRATCH/executable/fplot/plot ./ 
#rm plot.ppm
#./plot 5 1 40
