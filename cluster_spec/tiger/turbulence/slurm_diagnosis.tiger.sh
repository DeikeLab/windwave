#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --cpus-per-task=1
#SBATCH --time=00:60:00
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=jiarongw@princeton.edu

#The executable name
EXE=curved_uniform_forcing_moving_region
#Parameter value
LEVEL1=8
LEVEL2=7
RE=10000 #Default 40000
UstarRATIO=0.5
ak=0
#START=0
#END=50


export ScratchDir="/scratch/gpfs/jiarongw/turbulence/${EXE}_RE${RE}_Ustar${UstarRATIO}_ak${ak}_LEVEL${LEVEL1}"
echo $ScratchDir
cd $ScratchDir
#cp ../curved_uniform_forcing_RE40000/dump46 ./restart
python call.py 
