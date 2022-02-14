#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=40
#SBATCH --cpus-per-task=1
#SBATCH --time=00:60:00
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=jiarongw@princeton.edu

#The executable name
EXE=curved_fixREtau
#Parameter value
LEVEL=10
RE_tau=720
BO=200
g=4
ak=0.2
TIME=500
emaxRATIO=0.3
alterMU=8


export ScratchDir="/scratch/gpfs/jiarongw/turbulence/${EXE}_precursor_REtau${RE_tau}_ak${ak}_LEVEL${LEVEL}_emax${emaxRATIO}"
echo $ScratchDir
cd $ScratchDir

# Output data on uniform grid command
mkdir ./field
mkdir ./eta
cp $SCRATCH/executable/fdiagnosis_run_interface/diagnosis_run_interface ./

OUTLEVEL=9
for i in `seq 0 5`
do
    let Snapshot=50+i
    echo $Snapshot
    srun ./diagnosis_run_interface $RE_tau $BO $LEVEL $g $ak $TIME $emaxRATIO $alterMU $Snapshot $OUTLEVEL > message_diagnosis 2>&1 
done

#srun ./diagnosis_run $RE_tau $BO $LEVEL $g $ak $TIME $emaxRATIO $alterMU 75 $OUTLEVEL > message_diagnosis 2>&1 

#bash ./call.sh
#module load python
#python call.py

# Plotting command
#cp $SCRATCH/executable/fplot/plot ./ 
#rm plot.ppm
#./plot 5 1 40
