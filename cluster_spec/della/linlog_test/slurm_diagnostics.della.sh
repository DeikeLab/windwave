#!/bin/bash
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=16
#SBATCH --time=00:10:00
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=jiarongw@princeton.edu

#The executable name
EXE1=wavewind_linlog
EXE2=norun_multigrid
#Parameter value
LEVEL=11
ak=0.05
BO=3.34
RE=31000. #Default 40000
m=5
B=0
UstarRATIO=0.8


export ScratchDir="/scratch/gpfs/jiarongw/miscellaneous/linlog_m${m}B${B}Ustar${UstarRATIO}ak${ak}Bo${BO}Re${RE}LEVEL${LEVEL}"
echo $ScratchDir
#Copy the executable from executable directory
cp /scratch/gpfs/jiarongw/executable/f$EXE2/$EXE2 $ScratchDir
cd $ScratchDir
rm -r ./diagnostics
mkdir ./diagnostics
START=0
END=200


for i in $(seq 0 $END); do
    echo ${i}
    srun ./$EXE2 $LEVEL $ak $BO $RE $m $B $UstarRATIO ${i}
    cat ./diagnostics/eta${i}_* > ./diagnostics/eta${i}.dat
    cat ./diagnostics/field_direct${i}_*.dat > ./diagnostics/field_direct${i}.dat
    rm ./diagnostics/eta${i}_*
    rm ./diagnostics/field_direct${i}_*.dat
done 

#To move the whole directory to /tigress
#cp -r $ScratchDir /tigress/jiarongw/
#rm -rf $ScratchDir
