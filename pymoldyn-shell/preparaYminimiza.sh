#!/bin/csh -f
#

# $1: nombre de la proteina o del conjunto de proteinas que se quiera revisar
# $2: nombre del residuo de la proteina que se quiere centrar en la bicapa (si se pone cero se centra toda la proteina)
#
# EJEMPLOS: 
#	 ./preparaYminimiza.sh a1 0
#	 ./preparaYminimiza.sh a2a 12

foreach x ($1)
		mkdir {$x}_$3
		cd {$x}_$3
		cp ../x4bilayer.pdb .
		cp ../{$x}.pdb ../hoh* ../ions* ../$3* .
		cp ../ffopls* .
		../inspopc.sh {$x}.pdb 2.0 0.8 $2 protbi.pdb 0 $3
		rm x4bilayer.pdb
#		cd ..
		mkdir Rmin
		cd Rmin
		mv ../topol.tpr .
		cp ../../runcola.sh {$x}min.sh
		#mdrun -v -s
#		qsub -l num_proc=1,s_rt=24:00:00,s_vmem=1G,h_fsize=10G -pe mpi 8 ./{$x}min.sh
		cd ../..
end
