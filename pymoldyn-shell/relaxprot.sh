#!/bin/csh -f
#
# $1: nombre de la proteina o del conjunto de proteinas 
# $2: numero correspondiente a cada una de las seed de cada proteina
# $3: constante de fuerza utilizada para las restricciones de los atomos pesados en kJ/mol/nm^2
# $4: tiempo de simulacion en picosegundos


# EJEMPLO: 
#	 ./relaxprot.sh "a1 a2 a3" "1 2 3" 800 1000
#	 ./relaxprot.sh "a3 a1" "1 2 3" 400 1000

module load gromacs/4.0.5
echo 1 > uno
@ time=`awk -v tempo=$4 '{time=tempo/0.002; print time}' uno `
rm uno

foreach x ($1)
	cd {$x}
	foreach y ($2)
		cd S{$y}/eq
		ls -lhtr >> ls
		grep drwxr ls > tmp
		@ nd=`awk 'END {print $NF+1}' tmp`
		rm ls tmp
		mkdir $nd
		mv *.* $nd
		mv nodelist $nd
		cd $nd
		echo 1 > uno
		@ ndmasuno=`awk -v var=$nd '{ndmasuno=var+1; print ndmasuno}' uno `
		rm uno
		# Moving posres*itp and *mdp files
		cp ../1/posre1000_1.itp ./posre{$3}_{$ndmasuno}.itp
		cp ../1/posre_hoh1000_1.itp ./posre_hoh{$3}_{$ndmasuno}.itp
		cp ../1/posre_ion1000_1.itp ./posre_ion{$3}_{$ndmasuno}.itp
		cp ../1/posre_lig1000_1.itp ./posre_lig{$3}_{$ndmasuno}.itp
		cp ../1/eq1000_1.mdp ./eq{$4}_{$ndmasuno}.mdp
		# Reducing the constant for posre.itp (protein)
		sed "s|1  1000  1000  1000|1  $3  $3  $3|" posre{$3}_{$ndmasuno}.itp > tmp
		mv tmp posre{$3}_{$ndmasuno}.itp
		cp posre{$3}_{$ndmasuno}.itp ../../../posre.itp
		# Reducing the constant for posre_hoh.itp (crystallographic waters)
		sed "s|1       1000       1000       1000|1  $3  $3  $3|" posre_hoh{$3}_{$ndmasuno}.itp > tmp
		mv tmp posre_hoh{$3}_{$ndmasuno}.itp
		cp posre_hoh{$3}_{$ndmasuno}.itp ../../../posre_hoh.itp
		# Reducing the constant for posre_ion.itp (ions)
		sed "s|1       1000       1000       1000|1  $3  $3  $3|" posre_ion{$3}_{$ndmasuno}.itp > tmp
		mv tmp posre_ion{$3}_{$ndmasuno}.itp
		cp posre_ion{$3}_{$ndmasuno}.itp ../../../posre_ion.itp
		# Reducing the constant for posre_lig.itp (ligand)
		sed "s|1       1000       1000       1000|1  $3  $3  $3|" posre_lig{$3}_{$ndmasuno}.itp > tmp
		mv tmp posre_lig{$3}_{$ndmasuno}.itp
		cp posre_lig{$3}_{$ndmasuno}.itp ../../../posre_lig.itp
		# Changing eq.mdp, in order to specify simulation time
		sed "s|nsteps              =  250000    ; total 0.5 ns|nsteps = $time|" eq{$4}_{$ndmasuno}.mdp > tmp
		mv tmp eq{$4}_{$ndmasuno}.mdp
		sed "s|gen_vel             =  yes|gen_vel = no|" eq{$4}_{$ndmasuno}.mdp > tmp
		mv tmp eq{$4}_{$ndmasuno}.mdp
		grompp -f eq{$4}_{$ndmasuno}.mdp -c confout.gro -n ../../../index.ndx -p ../../../topol.top -o trelax.tpr
		cp trelax.tpr ../topol.tpr
		# Moving *mdp and *itp files one folder upwards
		mv eq{$4}_{$ndmasuno}.mdp ..
		mv posre{$3}_{$ndmasuno}.itp ../posre{$3}_{$ndmasuno}.itp
		mv posre_hoh{$3}_{$ndmasuno}.itp ../posre_hoh{$3}_{$ndmasuno}.itp
		mv posre_ion{$3}_{$ndmasuno}.itp ../posre_ion{$3}_{$ndmasuno}.itp
		mv posre_lig{$3}_{$ndmasuno}.itp ../posre_lig{$3}_{$ndmasuno}.itp
		sed 's/mdrun_mpi -v/mdrun_mpi -v -cpi/' ../../../../run_ib.sh > ../{$x}S{$y}.sh
		cd ..
                qsub {$x}S{$y}.sh
		#sbatch -p batch -N2 -n16 ./{$x}S{$y}.sh
		#qsub -P icts132 -l num_proc=1,s_rt=24:00:00,s_vmem=4G,h_fsize=10G,h_stack=16M -pe mpi 16 ./{$x}S{$y}.sh
		cd ../../
	end
	cd ..
end

