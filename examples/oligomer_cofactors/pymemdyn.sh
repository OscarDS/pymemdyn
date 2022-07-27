#!/bin/bash -l

#SBATCH -N 1
#SBATCH -n 32
#SBATCH -t 3-00:00:00

conda activate pymemdyn

pymemdyn -p 10mer_prot.pdb -i 10mer_ions -w 10mer_water --res ca
