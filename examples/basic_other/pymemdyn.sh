#!/bin/bash -l

#SBATCH -N 1
#SBATCH -n 32
#SBATCH -t 3-00:00:00

conda activate pymemdyn

pymemdyn -p sc6a2_human.pdb --res ca