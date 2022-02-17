#!/bin/bash
#PBS -N mup5L10
#PBS -A UPRI0020
#PBS -l walltime=00:30:00
#PBS -q regular
#PBS -m abe
#PBS -M jiarongw@princeton.edu
#PBS -l select=1:ncpus=36:mpiprocs=36

#The executable name
EXE=curved_fixREtau
#Parameter value
LEVEL=9
RE_tau=180 #Default 40000
BO=10
g=4
ak=0.1
TIME=250
emaxRATIO=0.3

export ScratchDir="/glade/scratch/jiarongw/turbulence/${EXE}_REtau${RE_tau}_BO${BO}_g${g}_ak${ak}_LEVEL${LEVEL}_emax${emaxRATIO}_precursor"
echo $ScratchDir
#rm -rf $ScratchDir
#mkdir -p $ScratchDir
#cp /glade/u/home/jiarongw/windwave/cluster_spec/cheyenne/executable/curved

#cp /glade/scratch/jiarongw/executable/f$EXE/$EXE $ScratchDir 
#cp /home/jiarongw/windwave/project_spec_source/turbulence/curved/${EXE}.c $ScratchDir
cd $ScratchDir
#cp ../curved_fixREtau_REtau_BO10_g4_ak0.1_LEVEL8_emax0.5/dump180 ./restart
#cp ../curved_uniform_forcing_moving_RE10000_Ustar0.5_ak0.1_LEVEL9_04_refinewater_precursor/dump45 ./restart
#cp ../curved_fixREtau_REtau1800_BO10_g4_ak0.1_LEVEL8_emax0.3/dump80 ./restart
#cp ../curved_fixREtau_REtau720_BO10_g4_ak0.1_LEVEL9_emax0.3/dump245 ./restart
#cp ./dump100 ./restart
mpirun -np 36 ./curved $RE_tau $BO $LEVEL $g $ak $TIME $emaxRATIO > message 2>&1 
