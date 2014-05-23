#!/bin/bash -f
#
# This script was designed to work with GROMACS 4.0.3
#
# It introduces a membrane protein dimer into a water-solvated POPC-bilayer, all inside of a hexagonal prism. It considers the total charge of the protein and add ions to neutralize it. In addition to the ions necessary to neutralize the system, it adds 0 Cl- and 0 Na+ to allow their interaction with the system (change accordingly if you want salt concentration, i.e. 28 for 0.1M). The OPLS-topology for the protein is generated and the script prepares the tpr-binary to minimize the whole system. The topol.top file specifies the number of POPC molecules in each leaflet of the bilayer.
#
# In order to see the final system BEFORE minimization open hexagon.pdb with a molecular viewer. If you want to minimize the system just execute "mdrun -s topol.tpr", preferibly in another directory.
#
#
# Please, contact to Angel Pineiro at angel.pineiro@usc.es for any question or comment
# Version: June 10, 2009
#
#
# This script NEEDS THE FOLLOWING FILES IN THE SAME DIRECTORY: x4bilayer.pdb (a bilayer of POPC equilibrated at 260 K), ffoplsaabon_mod.itp, ffoplsaa_mod.itp and ffoplsaanb_mod.itp (opls ff modified following instructions provided by Chris Neale in the GROMACS mailing list).
#
#
# INPUTS:
# 1.- PDB Filename containing the protein structure oriented with the helices parallel to the z axis
# 2.- Minimum distance from any protein atom to the xy box walls
# 3.- Minimum distance from any protein atom to the z box walls
# 4.- A reference residue which CA will be in the center of the bilayer
# 5.- Output pdb file
#
# EXAMPLE: " ./dinspopc.sh protein.pdb 2.0 0.9 241 out.pdb "
#
#
 sed -e 's/HID /HISA/g' -e 's/HIE /HISB/g' -e 's/HIP /HISH/g' -e 's/OXT/O  /g' $1 > proteinhis.pdb
 #sed -i 's/OXT/O  /' proteinhis.pdb
 #grep -v TER proteinhis.pdb | grep -v ENDMDL > tmp
 #cat hoh_exp.pdb >>tmp
 #echo TER >> tmp
 #echo ENDMDL >> tmp
 #mv tmp proteinhis.pdb

pdb2gmx_d -f proteinhis.pdb -o proteinopls.pdb -p protein.top -ignh << EOF
14
4
EOF

 #rm proteinhis.pdb tmp
# awk '{if ($3=="Position") fin=5}; {if (($1!="#include") && (fin!=5)) print}' protein.top > protein.itp
 #rm protein.top

# The size of the box (in the xy plane and in the z dimension) is measured using cubic boxes.
# The same (a, b, c) box vectors will be employed for the hexagonal prism shaped box
# The volume of this box is 13.4% smaller than that of a cubic box with the same distance between periodic images
#
 editconf_d -f proteinopls.pdb -o proteinopls.pdb -d $2
 xyprotein=`awk '{if (($1=="CRYST1") && ($2>$3)) print $2}; {if (($1=="CRYST1") && ($3>$2)) print $3}' proteinopls.pdb`
 editconf_d -f proteinopls.pdb -o proteinopls.pdb -d $3
 hprotein=`awk '{if ($1=="CRYST1") print $4}' proteinopls.pdb`

# The dimensions of the bilayer are measured
 xbilayer=`awk '{if ($1=="CRYST1") print $2}' x4bilayer.pdb`
 ybilayer=`awk '{if ($1=="CRYST1") print $3}' x4bilayer.pdb`
 hbilayerw=`awk '{if ($1=="CRYST1") print $4}' x4bilayer.pdb`

# A box containing only the POPC molecules is created and its height is measured
 grep POP x4bilayer.pdb > popc.pdb
 editconf_d -f popc.pdb -o popc.pdb -d 0 
 zbilayer=`awk '{if ($1=="CRYST1") print $4}' popc.pdb`

# Conversion from Anstromgs to GROMACS units (nanometers)
 echo 1 > uno
 gmxxbilayer=`awk -v x=$xbilayer '{xbilayer=x/10; print xbilayer}' uno `
 gmxybilayer=`awk -v y=$ybilayer '{ybilayer=y/10; print ybilayer}' uno `
 gmxhbilayerw=`awk -v h=$hbilayerw '{hbilayerw=h/10; print hbilayerw}' uno `
#******************************************************************************
 gmxzbilayer=`awk -v z=$zbilayer '{zbilayer=z/9; print zbilayer}' uno `
#******************************************************************************
 gmxxyprotein=`awk -v xy=$xyprotein '{xyprotein=xy/10; print xyprotein}' uno `
 gmxhprotein=`awk -v h=$hprotein '{hprotein=h/10; print hprotein}' uno `

# Introducing the protein in a hexagonal prism shaped box and the dried bilayer in a cubic box
 editconf_d -f proteinopls.pdb -o proteinopls.pdb -box $gmxxyprotein $gmxxyprotein $gmxzbilayer -angles 90 90 120 -bt tric
 editconf_d -f popc.pdb -o popc.pdb -box $gmxxbilayer $gmxybilayer $gmxzbilayer
 
#creating posre_hoh.itp, for a gradual reduction of force constants
cat > posre_hoh.itp << EOF
; position restraints for crystallographic waters (resn HOH)

[ position_restraints ]
;  i funct       fcx        fcy        fcz
   1    1       1000       1000       1000
EOF

# Creating a topol.top file 
cat > topol.top << estodo
#include "ffoplsaa_mod.itp"
#include "protein_Protein_chain_A.itp"
; Include Position restraint file
#ifdef POSRES
#include "posre_Protein_chain_A.itp"
#endif
#include "protein_Protein_chain_B.itp"
; Include Position restraint file
#ifdef POSRES
#include "posre_Protein_chain_B.itp"
#endif
#include "../popc.itp"
#include "ions.itp"
#include "spc.itp"
[ system ]
; Name
Membrane protein
[ molecules ]
; Compound  #mols
Protein_chain_A 1
Protein_chain_B 1
estodo

# Taking the z coordinate of the reference-residue CA and 
# Moving the protein to put the reference residue in the center of the bilayer

if [ "$4" -gt "0" ]; then
 grep "   $4    " proteinopls.pdb | grep CA  > tmp
 zCA=`awk -v zcenter=$gmxzbilayer '{zCA=-(($8/10)-(zcenter/2)); print zCA}' tmp `
 editconf_d -f proteinopls.pdb -o proteinopls.pdb -translate 0 0 $zCA
 rm tmp
fi

if [ "$4" -eq "0" ]; then
 zCA="0"
fi

# Solvating the protein with POPC
 genbox_d -cp proteinopls.pdb -cs popc.pdb -o protpopc.pdb

# Solvating the protein-POPC system with water
 ini=`awk -v hw=$hbilayerw -v hp=$hprotein '{delta_2=(hw-hp)/2; print delta_2}' uno `
 fin=`awk -v princ=$ini -v hp=$hprotein '{final=princ+hp; print final}' uno `
 awk -v start=$ini -v end=$fin '{if ($3=="OW" && ($8>end || $8<start)) res=$5; if (($5!=res) && ($4=="SOL")) print;}' x4bilayer.pdb > water.pdb
 editconf_d -f water.pdb -o water.pdb -box $gmxxbilayer $gmxybilayer $gmxhprotein
 editconf_d -f protpopc.pdb -o protpopc.pdb -box $gmxxyprotein $gmxxyprotein $gmxhprotein -angles 90 90 120 -bt tric
 genbox_d -cp protpopc.pdb -cs water.pdb -o tmp.pdb

# Determining the number of lipids in each leaflet
 hprotein_2=`awk -v h=$hprotein '{half=h/2; print half}' uno `
 rm popc.pdb
 grep N4 tmp.pdb > popc.pdb
 awk -v half=$hprotein_2 '{if ((NF==11) && ($9>half)) print}; {if ((NF==10) && ($8>half)) print}' popc.pdb > up.pdb
 awk -v half=$hprotein_2 '{if ((NF==11) && ($9<half)) print}; {if ((NF==10) && ($8<half)) print}' popc.pdb > down.pdb
 nup=`awk 'END {print NR}' up.pdb`
 ndown=`awk 'END {print NR}' down.pdb`

# Writing the number of lipid and solvent molecules to the topol.top file
 echo "; Number of POPC molecules in the layer with higher z-coordinate value:" >> topol.top
 echo POPC $nup >> topol.top
 echo "; Number of POPC molecules in the layer with lower z-coordinate value:" >> topol.top
 echo POPC $ndown >> topol.top

# grep 'OW&&SOL' tmp.pdb > ow.pdb
 `awk '{if (($3=="OW") && ($4=="SOL")) print}' tmp.pdb > ow.pdb`
 nwat=`awk 'END {print NR}' ow.pdb`
 echo "; Total number of water molecules:" >> topol.top
 echo SOL $nwat >> topol.top
 rm water.pdb up.pdb down.pdb protpopc.pdb popc.pdb ow.pdb

# Generating a simple mdp file for structure minimization
 cat > steep.mdp << estodo
	constraints = none
        integrator      = steep
        nsteps          = 500
        nstlist         = 1
        ns-type         = Grid
        pbc             = xyz
        rlist           = 0.8
        coulombtype     = cut-off
        rcoulomb        = 1.4
        vdw-type        = cut-off
        rvdw            = 1.4
        nstenergy       = 20
estodo

# # Generating the tpr binary that will be used to add ions
grompp_d -f steep.mdp -c tmp.pdb 
#
trjconv_d -f tmp.pdb -center -o tmp.pdb -pbc mol << EOF
1
0
EOF

#
grompp_d -f steep.mdp -c tmp.pdb &> tmp

# Neutralize the total charge by adding ions

# grep "total charge" tmp > tmp1
# mv tmp1 tmp
#  charge=`awk 'END {charge=$NF*10/10 ; print int(charge)}' tmp`

# echo "Protein charge = " $charge


# if [ "$charge" -lt "0" ]; then
# charge=`awk -v x=$charge -v nions=$6 '{carga=-x; print carga+nions}' uno `
# echo 14 | genion_d -s -o $5 -np $charge -nn $6 -pname NA+ -nname CL- -p &> tmp1
# fi

# if [ "$charge" -gt "0" ]; then
# charge=`awk -v x=$charge -v nions=$6 '{carga=x; print carga+nions}' uno `
# echo 14 | genion_d -s -o $5 -nn $charge -np $6 -nname CL- -pname NA+ -p &> tmp1
# fi

# if [ "$charge" -eq "0" ]; then
# echo 14 | genion_d -s -o $5 -nn $6 -np $6 -nname CL- -pname NA+ -p &> tmp1
# fi

# As we already know the charge of both proteins (10), we will introduce it manually
# The 'grep' command only takes the charge of one protein (not both at the same time)
echo 15 | genion_d -s -o $5 -nn 19 -np 0 -nname CL- -pname NA+ -p &> tmp1

# # Centering the bilayer in the box
 grompp_d -f steep.mdp -c $5
 echo 0 | trjconv_d -f $5 -pbc mol -o $5 -trans 0 0 $zCA

# # Generating the tpr binary for minimization (the input for the mdrun file)
 grompp_d -f steep.mdp -c $5

# # Generating a pdb file of the structure prior to minimization, just for visualization. Don't use this pdb to run simulations!

echo 1 0 | trjconv_d -f $5 -ur compact -s -pbc mol -center -o hexagon.pdb 

#rm tmp tmp1 uno tmp.pdb genion.log

echo $zCA

