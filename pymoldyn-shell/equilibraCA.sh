#!/bin/csh -f
#
# $1: nombre de la proteina o del conjunto de proteinas que se quiera extender
# $2: numero correspondiente a cada una de las seed de cada proteina

# EJEMPLO: 
#	 ./equilibraCA.sh " a1 a2 a3 " " 1 2 3 "
#	 ./equilibraCA.sh " a3 a1 " " 1 2 3 "

module load gromacs/4.0.5
foreach x ($1)
	cd {$x}
	foreach y ($2)
		cp S{$y}/eq/conf*.gro ./eq.gro
		mkdir S{$y}/eqCA
		cp ../{$1}ca200.itp ./posre.itp
		grompp -f ../eqCA.mdp -c eq.gro -n
		mv topol.tpr S{$y}/eqCA
		cp ../eqCA.mdp S{$y}/eqCA
		cp ../run_ib.sh S{$y}/eqCA/{$x}S{$y}CA.sh
		cp posre.itp S{$y}/eqCA/posreCA_1.itp
		cp ../eqCA.mdp S{$y}/eqCA/eqCA_1.mdp
		cd S{$y}/eqCA
		qsub {$x}S{$y}CA.sh
                #sbatch -p batch -N2 -n16 ./{$x}S{$y}CA.sh
#		qsub -P icts132 -l num_proc=1,s_rt=90:00:00,s_vmem=1G,h_fsize=10G,h_stack=16M -pe mpi 16 ./{$x}S{$y}CA.sh
		cd ../../
	end
	cd ..
end

