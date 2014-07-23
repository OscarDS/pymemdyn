#!/bin/csh -f

# $1: nombre de la proteina o del conjunto de proteinas que se quiera extender
# $2: numero correspondiente a cada una de las seed de cada proteina
# $3: numero de picosegundos que se quieren simular. La trayectoria se extendera en esta cantidad de picosegundos.
#
# Los archivos se almacenan en un directorio cuyo nombre es un numero correspondiente al numero de extensiones realizadas sobre esa simulacion utilizando este script
#
# EJEMPLO: si quiero extender durante 5 ns mas =>
#	 ./extendprod "a1 a2 a3" "1 2 3" 5000
#	 ./extendprod "a1 a3" "3" 5000
#
# no olvidar poner comillas si se quieren extender las simulaciones de varias proteinas/seeds
#
module load gromacs/4.0.5
foreach x ($1)
	cd {$x}
	foreach y ($2)
		cd S{$y}/prod
		@ nd=`ls -l |grep -c dr` + 1
		mkdir $nd
		mv *.* $nd
		#mv nodelist $nd
		cd $nd
		cp state.cpt ..
		tpbconv -s -o -extend $3
		mv tpxout.tpr ../topol.tpr
		#cp ../../../../runcolacpt.sh ../{$x}S{$y}.sh
		cd ..
		qsub ../$nd/{$x}S{$y}.sh
		#qsub -P icts132 -l num_proc=1,s_rt=300:00:00,s_vmem=4G,h_fsize=10G,h_stack=16M -pe mpi_4p 16 ./{$x}S{$y}.sh
		cd ../..
	end
	cd ..
end

