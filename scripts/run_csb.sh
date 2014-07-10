#!/bin/bash
echo "#!/bin/bash -l
source home/apps/gromacs-4.6.5/bin/GMXRC
source /home/apps/bin/apps.sh
run.py -p a2b_filter.pdb
" > temp.sh
chmod +x temp.sh
sbatch -A snic2013-26-1 -n 8 --exclusive -t 48:00:00 -J pymemdyn temp.sh
