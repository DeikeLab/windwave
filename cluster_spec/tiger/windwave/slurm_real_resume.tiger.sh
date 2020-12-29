#!/bin/bash
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=40
#SBATCH --time=4:00:00
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=jiarongw@princeton.edu

#The executable name
EXE=wavewind_real_linlog_coarser
#Parameter value
LEVEL=11
ak=0.05
UstarRATIO=0.4
Nwave=1
L0=0.03
TEND=50
REGION=0.05

export ScratchDir="/scratch/gpfs/jiarongw/windwave/${EXE}_waves${Nwave}L0${L0}Ustar${UstarRATIO}ak${ak}Bo${BO}Re${RE}LEVEL${LEVEL}REGION${REGION}_comparetolinear"
echo $ScratchDir
cd $ScratchDir
cp ./dump0.3 ./restart
srun ./$EXE $LEVEL $ak $UstarRATIO $Nwave $L0 $TEND $REGION > message 2>&1

#To move the whole directory to /tigress
#cp -r $ScratchDir /tigress/jiarongw/
#rm -rf $ScratchDir
