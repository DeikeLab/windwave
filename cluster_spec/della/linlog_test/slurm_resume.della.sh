#!/bin/bash
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=20
#SBATCH --time=04:00:00
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=jiarongw@princeton.edu

#The executable name
EXE=wavewind_linlog_adaptive_limited
#Parameter value
LEVEL=11
ak=0.05
BO=200
RE=200000. #Default 40000
m=5
B=0
UstarRATIO=0.4


export ScratchDir="/scratch/gpfs/jiarongw/miscellaneous/linlog_adaptive_limited_m${m}B${B}Ustar${UstarRATIO}ak${ak}Bo${BO}Re${RE}LEVEL${LEVEL}"
echo $ScratchDir

#Copy the executable from executable directory
cp /scratch/gpfs/jiarongw/executable/f$EXE/$EXE $ScratchDir
cd $ScratchDir
cp ./dump9.125 ./restart
srun ./$EXE $LEVEL $ak $BO $RE $m $B $UstarRATIO

#To move the whole directory to /tigress
#cp -r $ScratchDir /tigress/jiarongw/
#rm -rf $ScratchDir
