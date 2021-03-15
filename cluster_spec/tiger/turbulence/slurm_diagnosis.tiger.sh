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
LEVEL=9
RE_tau=180 #Default 40000
BO=200
g=1
ak=0.1
emaxRATIO=0.3


export ScratchDir="/scratch/gpfs/jiarongw/turbulence/${EXE}_REtau${RE_tau}_BO${BO}_g${g}_ak${ak}_LEVEL${LEVEL}_emax${emaxRATIO}"
echo $ScratchDir
cd $ScratchDir

# Output data on uniform grid command
mkdir ./field
mkdir ./eta
cp $SCRATCH/executable/fdiagnosis/diagnosis ./
bash ./call.sh
#module load python
#python call.py

# Plotting command
#cp $SCRATCH/executable/fplot/plot ./ 
#rm plot.ppm
#./plot 5 1 40
