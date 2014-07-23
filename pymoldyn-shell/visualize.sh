#!/bin/bash
# This command is often needed to get a file in the correct format for visualization 
# in pymol.
echo 1 ; echo 0 | trjconv -f traj.xtc -pbc mol -ur compact -s topol.tpr -o trajpbc.xtc -center
