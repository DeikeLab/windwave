#!/bin/bash
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=40
#SBATCH --time=24:00:00
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=jiarongw@princeton.edu

#The executable name
EXE=wavewind_linlog_withforcing
#Parameter value
LEVEL=11
ak=0.2
BO=200
RE=100000. #Default 40000
UstarRATIO=0.6
Nwave=1
g=1
TEND=50


export ScratchDir="/scratch/gpfs/jiarongw/windwave/linlog_withforcing_${Nwave}waves_Ustar${UstarRATIO}ak${ak}Bo${BO}Re${RE}LEVEL${LEVEL}"
echo $ScratchDir
cd $ScratchDir
cp ./dump12.7 ./restart
srun ./$EXE $LEVEL $ak $BO $RE $UstarRATIO $Nwave $g $TEND

#To move the whole directory to /tigress
#cp -r $ScratchDir /tigress/jiarongw/
#rm -rf $ScratchDir
