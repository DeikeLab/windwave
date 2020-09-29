#!/bin/bash
#SBATCH --nodes=2
#SBATCH --ntasks=80               # total number of tasks across all nodes
#SBATCH --cpus-per-task=1        # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --time=04:00:00
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
REGION=0.05

export ScratchDir="/scratch/gpfs/jiarongw/windwave/${EXE}_coarser_waves${Nwave}L0${L0}Ustar${UstarRATIO}ak${ak}Bo${BO}Re${RE}LEVEL${LEVEL}REGION${REGION}_testing2"
echo $ScratchDir
rm -rf $ScratchDir
mkdir -p $ScratchDir
#Copy the executable from executable directory
cp /scratch/gpfs/jiarongw/executable/f$EXE/$EXE $ScratchDir
cd $ScratchDir
mkdir ./pressure
srun ./$EXE $LEVEL $ak $UstarRATIO $Nwave $L0 $TEND $REGION > message 2>&1

#To move the whole directory to /tigress
#cp -r $ScratchDir /tigress/jiarongw/
#rm -rf $ScratchDir
