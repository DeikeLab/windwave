#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=40
#SBATCH --cpus-per-task=1
#SBATCH --time=00:60:00
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=jiarongw@princeton.edu

#The executable name
EXE=onephase
#Parameter value
LEVEL1=9
LEVEL2=8
RE=10000 #Default 40000
UstarRATIO=0.5
ak=0.1
#START=0
#END=50


export ScratchDir="/scratch/gpfs/jiarongw/turbulence/channel/${EXE}_test3"
echo $ScratchDir
cd $ScratchDir

# Output data on uniform grid command
cp $SCRATCH/executable/fdiagnosis_onephase/diagnosis_onephase ./
python call.py

# Plotting command
#cp $SCRATCH/executable/fplot/plot ./ 
#rm plot.ppm
#./plot 5 1 40
