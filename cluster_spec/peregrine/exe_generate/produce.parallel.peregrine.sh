#!/bin/bash
# This is the compile script for Della


OUTFOLDER="/home/jiarongw/scratch/executable/f$1"
# Note that space is needed between = and []
if [ "$2" = "clean" ]
then
    rm -r $OUTFOLDER
    echo "Directory cleaned!"
else
    mkdir -p $OUTFOLDER
    make -f Makefile.parallel.peregrine OUTPUT=$1
    #rm the copy of the source code to avoid version error 
    rm ./$1.c
    rm ./_$1.c
fi
