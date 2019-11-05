#!/bin/bash
OUTFOLDER="$WORK/f$1"
# Note that space is needed between = and []
if [ "$2" = "clean" ]; then
    rm -r $OUTFOLDER
    echo "Directory cleaned!"
fi
mkdir $OUTFOLDER
make -f Makefile.parallel OUTPUT=$1
