#!/bin/csh -f
#
# $1: nombre de la proteina o del conjunto de proteinas que se quiera extender
# $2: numero correspondiente a cada una de las seed de cada proteina

# EJEMPLO: 
#	 ./equilibra.sh " a1 a2 a3 " " 1 2 3 "
#	 ./equilibra.sh " a3 a1 " " 1 2 3 "

#
# HUGO: Note changes in ndx file!!
#
cat > ndx << estodo
1 || r LIG
name 19 protlig
r HOH || a Na* || r SOL || a Cl*
name 20 wation
q
estodo
module load gromacs/4.0.5
#set gromacspath = /usr/local/bin
foreach x ($1)
	cd {$x}
	editconf -f Rmin/confout.gro -o min.pdb
	cat ../ndx | make_ndx -f min.pdb 
	echo 3 | genrestr -f Rmin/topol.tpr -fc 200 200 200 -o ../${x}ca200.itp
	foreach y ($2)
		mkdir S{$y}
		mkdir S{$y}/eq
		grompp -f ../eq.mdp -c min.pdb -n -debug 1
		mv topol.tpr S{$y}/eq
#		cp ../eq.mdp S{$y}/eq
		cp ../run_ib.sh S{$y}/eq/{$x}S{$y}.sh
		cp posre.itp S{$y}/eq/posre1000_1.itp
		cp posre_hoh.itp S{$y}/eq/posre_hoh1000_1.itp
		cp posre_ion.itp S{$y}/eq/posre_ion1000_1.itp
                cp posre_lig.itp S{$y}/eq/posre_lig1000_1.itp
		cp ../eq.mdp S{$y}/eq/eq1000_1.mdp
		cd S{$y}/eq
#		mdrun -v -s
		#qsub -P icts132 -l num_proc=1,s_rt=24:00:00,s_vmem=4G,h_fsize=10G,h_stack=16M -pe mpi 16 ./{$x}S{$y}.sh
		#sbatch -p batch -N2 -n16 ./{$x}S{$y}.sh
		qsub run_ib.sh
                cd ../../
	end
	cd ..
end

