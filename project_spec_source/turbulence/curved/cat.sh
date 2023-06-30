#!/bin/bash

for i in `seq 0 13`;
do
    let t=91+i*1
    echo $t
    cat eta_t${t}_* >> eta_t${t}
    rm eta_t${t}_*
    for j in `seq 1 9`
    do
	cat eta_t${t}.${j}_* >> eta_t$t.$j
	rm eta_t${t}.${j}_*
    done
done
