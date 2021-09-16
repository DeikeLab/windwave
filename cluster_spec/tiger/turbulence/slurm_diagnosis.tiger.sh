#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --time=00:30:00
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=jiarongw@princeton.edu

#The executable name
EXE=curved_fixREtau_boundary
#Parameter value
LEVEL=10
RE_tau=720
BO=200
g=4
ak=0.2
TIME=57
emaxRATIO=0.3
alterMU=8

export ScratchDir="/scratch/gpfs/jiarongw/turbulence/${EXE}_REtau${RE_tau}_BO${BO}_g${g}_ak${ak}_MU${alterMU}_LEVEL${LEVEL}_emax${emaxRATIO}"
echo $ScratchDir
cd $ScratchDir

# Output data on uniform grid command
# mkdir ./field
# mkdir ./eta
cp $SCRATCH/executable/fdiagnosis/diagnosis ./

OUTLEVEL=9
for i in `seq 0 1`
do
    let Snapshot=100+i
    srun ./diagnosis $Snapshot $OUTLEVEL > message_diagnosis 2>&1
done 

#bash ./call.sh
#module load python
#python call.py

# Plotting command
#cp $SCRATCH/executable/fplot/plot ./ 
#rm plot.ppm
#./plot 5 1 40
