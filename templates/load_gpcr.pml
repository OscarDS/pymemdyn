#!/usr/bin/env python
set defer_builds_mode, 3
load hexagon.pdb, ini-state
load hexagon.pdb, equi
color grey70, ini-state
hide lines, resn pop
hide lines, (element h and neighbor element c)
select protein, chain A and equi
select solvent, resn SOL
select membrane, resn POP
select memblimi, equi and name n4+p8+na*
hide everything, membrane
hide everything, solvent
show cartoon, protein
show spheres, solvent and elem O
show spheres, memblimi
set sphere_scale, 0.2
cmd.spectrum("count",selection="(protein)&e. c")

### If you have a ligand ###
select ligand, equi and resn L01
show sticks, ligand and not (element h and neighbor element c)
select bsite, br. (equi and not solvent within 5 of ligand)
show lines, bsite and not (element h and neighbor element C)
util.cba(154,"ligand",_self=cmd)
cmd.disable('ligand')
### End of ligand settings ###


### Load trajectory and align using c-alphas ###
load traj_pymol.xtc, equi
intra_fit equi and name ca

### Center on ligand ###
center ligand
zoom ligand

### If the ligand and the protein are too close in the initial 
### pymol draws non-existing bonds based on distance.
unbond protein, ligand

set auto_zoom, 0
deselect
