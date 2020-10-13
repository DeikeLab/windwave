#!/bin/bash
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=40
#SBATCH --time=4:00:00
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=jiarongw@princeton.edu

#The executable name
EXE=wavewind_real_linlog_ve
#Parameter value
LEVEL=11
ak=0.1
UstarRATIO=0.6
Nwave=4
L0=0.25
TEND=50
REGION=0.2

export ScratchDir="/scratch/gpfs/jiarongw/windwave/${EXE}_coarser_waves${Nwave}L0${L0}Ustar${UstarRATIO}ak${ak}Bo${BO}Re${RE}LEVEL${LEVEL}REGION${REGION}"
echo $ScratchDir
cd $ScratchDir
cp ./dump2.6 ./restart
srun ./$EXE $LEVEL $ak $UstarRATIO $Nwave $L0 $TEND $REGION

#To move the whole directory to /tigress
#cp -r $ScratchDir /tigress/jiarongw/
#rm -rf $ScratchDir
