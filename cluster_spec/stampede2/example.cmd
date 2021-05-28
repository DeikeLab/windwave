#!/bin/bash
#SBATCH -p skx-normal
## # -p normal  #if running on KNL. Otherwise skx-normal will run on SKX nodes
#SBATCH -J example
#SBATCH -o example.o%j
#SBATCH -e example.e%j
#SBATCH -N 2
#SBATCH --ntasks-per-node 48
## # --ntasks-per-node 68 # if running on KNL. Otherwise the skx-normal queue uses 48 tasks per node.
#SBATCH -t 11:00:00
## # You may also need to specify the account if you have access to more than one. In this case use: # -A BB-BB012345 # replacing BB-BB012345 with your account number, whatever it is.

ibrun ./example
