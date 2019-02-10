LEVEL=11
ak=0.05
BO=200
RE=40000.
m=8
B=0
UstarRATIO=1

export ScratchDir="/scratch/jiarongw/parameter/m$m\_B$B\_Ustar$UstarRATIO\_ak$ak"
echo $ScratchDir
rm -rf $ScratchDir
mkdir -p $ScratchDir

cp ./wavewind_serial $ScratchDir
#cp ./wavewind_parallel $ScratchDir
cd $ScratchDir
./wavewind_serial $LEVEL $ak $BO $RE $m $B $UstarRATIO
#./wavewind_parallel $LEVEL $ak $BO $RE $m $B $UstarRATIO 
