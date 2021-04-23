#!/bin/bash -l
#source /home/apps/gromacs-4.6.5/bin/GMXRC.bash
#source /home/apps/gromacs-4.6.7/bin/GMXRC.bash
ml gromacs/4.6.7
source /home/apps/apps/gromacs/4.6.7/bin/GMXRC.bash
#source /home/apps/OLD/bin/apps.sh
which python
/home/gpcradmin/gpcr_modsim/pymemdyn/run_db.py $*

