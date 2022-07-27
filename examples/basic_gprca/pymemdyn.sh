#!/bin/bash -l

#SBATCH -N 1
#SBATCH -n 32
#SBATCH -t 3-00:00:00

conda activate pymemdyn

pymemdyn -p aa2br_human.pdb --res bw
