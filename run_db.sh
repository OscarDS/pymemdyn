#!/bin/bash -l
#source /home/apps/gromacs-4.6.5/bin/GMXRC.bash
#source /home/apps/gromacs-4.6.7/bin/GMXRC.bash
source /home/apps/gromacs/4.6.7/bin/GMXRC.bash
source /home/apps/OLD/bin/apps.sh
/home/gpcradmin/gpcr_modsim/pymemdyn/run_db.py $*

