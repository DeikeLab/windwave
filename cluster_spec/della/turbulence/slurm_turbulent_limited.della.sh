#!/bin/bash
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=20
#SBATCH --time=20:00:00
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=jiarongw@princeton.edu

#The executable name
EXE=turbulence_RE_limited
#Parameter value
LEVEL=5
RE=64000 #Default 40000
START=0
END=50


export ScratchDir="/scratch/gpfs/jiarongw/turbulence/adlim_RE${RE}LEVEL${LEVEL}START${START}END${END}"
echo $ScratchDir
rm -f $ScratchDir
mkdir -p $ScratchDir
cp /scratch/gpfs/jiarongw/executable/f$EXE/$EXE $ScratchDir
cd $ScratchDir
srun ./$EXE $RE $LEVEL $START $END 

#To move the whole directory to /tigress
#cp -r $ScratchDir /tigress/jiarongw/
#rm -rf $ScratchDir
