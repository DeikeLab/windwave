#!/bin/bash
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=40
#SBATCH --time=0:30:00
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=jiarongw@princeton.edu

LEVEL=11
ak=0.05
BO=200
RE=4000 #Default 40000
m=8
B=0
UstarRATIO=1

export ScratchDir="/scratch/gpfs/jiarongw/parameter/m${m}B${B}Ustar${UstarRATIO}ak${ak}RE${RE}"
echo $ScratchDir
#rm -rf $ScratchDir
#mkdir -p $ScratchDir
cp ./wavewind_parallel $ScratchDir

cd $ScratchDir
srun ./wavewind_parallel $LEVEL $ak $BO $RE $m $B $UstarRATIO

#cp -r $ScratchDir /tigress/jiarongw/
#rm -rf $ScratchDir
