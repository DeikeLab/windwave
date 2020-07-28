#!/bin/bash
#SBATCH --nodes=2
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
#SBATCH --time=00:10:00
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=jiarongw@princeton.edu

#The executable name
EXE=stokes_ml
#Parameter value
NLAYER=30
N=256
coeff=0.05
RE=40000
ak=0.35


export ScratchDir="/scratch/gpfs/jiarongw/multilayer/${EXE}_RE${RE}_nl${NLAYER}_N${N}_coeff${coeff}_ak${ak}"
echo $ScratchDir
rm -rf $ScratchDir
mkdir -p $ScratchDir
#Copy the executable from executable directory
cp /scratch/gpfs/jiarongw/executable/f$EXE/$EXE $ScratchDir
cp ${SCRATCH}/multilayer/F_kxky_${FIELD} $ScratchDir/F_kxky 
#cp -r ${SCRATCH}/multilayer/${EXE}_${FIELD}/dump59 $ScratchDir/restart
cd $ScratchDir
echo srun ./$EXE $RE $NLAYER $coeff $N $ak
srun ./$EXE $RE $NLAYER $coeff $N $ak

#To move the whole directory to /tigress
#cp -r $ScratchDir /tigress/jiarongw/
#rm -rf $ScratchDir
