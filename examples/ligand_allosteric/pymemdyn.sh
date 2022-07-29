#!/bin/bash -l

#SBATCH -N 1
#SBATCH -n 32
#SBATCH -t 3-00:00:00

conda activate pymemdyn

pymemdyn -p cnr1_human.pdb -l 9gf -a 9gl --res bw
