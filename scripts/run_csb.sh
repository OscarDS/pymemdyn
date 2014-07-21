#!/bin/bash
echo "#!/bin/bash -l
source /home/apps/gromacs-4.6.5/bin/GMXRC.bash
source /home/apps/bin/apps.sh
run.py -p a2a_ag.pdb -l lig  --cho cho
" > temp.sh
chmod +x temp.sh
sbatch -n 8 -t 47:59:00 -J pymemdyn temp.sh
