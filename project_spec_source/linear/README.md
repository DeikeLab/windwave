This folder contains source code for the 2D laminar wind paper.

Most of the analyses are done on cases that were run with wavewind_rerun_test.c
The major difference between wavewind_rerun and wavewind_rerun_test is that the 
latter has modified g accounting for capillary effect. The backward traveling wave
seen in most previous simulations (which causes the energy curve to wiggle) might be
due to the inconsistent initialization of wave speed c.

There are some other cases that are run for testing purposes. For example two waves in a box etc.
