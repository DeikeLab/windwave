#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --time=2:00:00
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=jiarongw@princeton.edu

LEVEL=11
ak=0.05
BO=200
RE=40000.
m=8
B=0
UstarRATIO=2.5

export ScratchDir="/scratch/jiarongw/parameter/m$m\_B$B\_Ustar$UstarRATIO\_ak$ak"
echo $ScratchDir
rm -rf $ScratchDir
mkdir -p $ScratchDir

cp ./wavewind $ScratchDir
cd $ScratchDir
./wavewind $LEVEL $ak $BO $RE $m $B $UstarRATIO 
