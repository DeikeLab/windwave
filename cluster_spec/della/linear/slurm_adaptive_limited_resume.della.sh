#!/bin/bash
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=20
#SBATCH --time=16:00:00
#SBATCH --mail-type=end
#SBATCH --mail-user=jiarongw@princeton.edu

# This is a slurm script for Della
#The executable name
EXE=wavewind_linear_adaptive_limited
#Parameter value
LEVEL=11
ak=0.05
BO=0.53
RE=20000. #Default 40000
m=5
B=0
UstarRATIO=1

export ScratchDir="/scratch/gpfs/jiarongw/miscellaneous/linear_adaptive_limited_m${m}B${B}Ustar${UstarRATIO}ak${ak}Bo${BO}Re${RE}LEVEL${LEVEL}"
echo $ScratchDir
#Copy the executable from executable directory
cp /scratch/gpfs/jiarongw/executable/f$EXE/$EXE $ScratchDir
cd $ScratchDir
cp ./dump2.4375 ./restart
srun ./$EXE $LEVEL $ak $BO $RE $m $B $UstarRATIO

#To move the whole directory to /tigress
#cp -r $ScratchDir /tigress/jiarongw/
#rm -rf $ScratchDir
