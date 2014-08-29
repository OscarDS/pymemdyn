#!/bin/bash
# Launcher to sbatch
source /home/apps/gromacs-4.6.5/bin/GMXRC.bash
source /home/apps/bin/apps.sh
/home/apps/bin/python2.7 /home/gpcruser/gpcr_modsim/pymemdyn/run_db.py $*
#python /home/slurm/pymoldyn/run_db.py $*
