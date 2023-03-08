#!/bin/bash -l

#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 00:10:00

# Software/Packages | Version | .bashrc
# Python            | 3.7     | ml anaconda
# GROMACS           | 2021    | ml gromacs/2021
# PyMemDyn          | 1.5     | export PYMEMDYN=/home/rbroek/apps/pymemdyn_1    .5
#                   |         | export PATH=$PYMEMDYN:$PATH

ml gromacs/2021
export PYMEMDYN=/home/rkupper/apps/pymemdyn # change path to your pymemdyn installation
export PATH=$PYMEMDYN:$PATH
source /home/rkupper/.bashrc 	# Change path to your .bashrc file 
conda activate py37 		# Activate conda env of your ligpargen installation

# run pymemdyn on protein called cnr1_human
# with ligand called 9gf_charged, with charge of -1 with 2 ligpargen optimizations
# with allosteric calld 9gl

pymemdyn -p cnr1_human.pdb -l 9gf_charged --llo 2 --llc -1 -a 9gl --res bw
