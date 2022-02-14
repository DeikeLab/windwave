#!/bin/bash                                             
#PBS -N mup5L10
#PBS -A UPRI0020
#PBS -l walltime=01:00:00
#PBS -q regular
#PBS -m abe
#PBS -M jiarongw@princeton.edu
#PBS -l select=1:ncpus=36:mpiprocs=36

 
EXE=curved_fixREtau_boundary
#Parameter value                                        
LEVEL=10
RE_tau=720
ak=0
TIME=500
emaxRATIO=0.1

export ScratchDir="/glade/scratch/jiarongw/turbulence/curved_fixREtau_precursor_REtau${RE_tau}_ak${ak}_LEVEL${LEVEL}_emax${emaxRATIO}"
echo $ScratchDir
cd $ScratchDir

# Output data on uniform grid command             
mkdir ./field            
mkdir ./eta                         
cp $SCRATCH/executable/fdiagnosis/diagnosis ./

OUTLEVEL=9
for i in `seq 0 3`
do
    let Snapshot=75+i
    mpirun -np 36 ./diagnosis $Snapshot $OUTLEVEL > message_diagnosis 2>&1
done
