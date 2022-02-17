#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=40
#SBATCH --cpus-per-task=1
#SBATCH --time=00:30:00
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=jiarongw@princeton.edu

#The executable name
EXE=onephase
#Parameter value
OUTLEVEL=8
#SNAPSHOT=180

export ScratchDir="/scratch/gpfs/jiarongw/turbulence/channel/${EXE}_180_LEVEL9"
echo $ScratchDir
cd $ScratchDir

# Output data on uniform grid command
cp $SCRATCH/executable/fdiagnosis_onephase/diagnosis_onephase ./

for i in `seq 0 5`;
do
    let t=260+i*10
    echo $t
    ./diagnosis_onephase $t $OUTLEVEL
done

# Plotting command
#cp $SCRATCH/executable/fplot/plot ./ 
#rm plot.ppm
#./plot 5 1 40
