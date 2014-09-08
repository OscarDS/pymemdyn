CONGRATULATIONS!!

You have performed a MD equilibration of your receptor, including lipids, water
molecules and counterions.
If you want more technical details see reference [1].

The performed equilibration includes the following stages:

STAGE           RESTRAINED ATOMS       FORCE CONSTANT          TIME
                    kJ·mol^(-1)·nm^(-2) ns
Minimization    -                       -                       (Max. 500 steps)
Eq1             Protein Heavy Atoms     1000                    0.5
Eq2             Protein Heavy Atoms     800                     0.5
Eq3             Protein Heavy Atoms     600                     0.5
Eq4             Protein Heavy Atoms     400                     0.5
Eq5             Protein Heavy Atoms     200                     0.5
Eq6             Protein C-alfa Atoms    200                     2.5

In this folder you will find several files related to this simulation:
INPUTS:
- popc.itp              # Topology of the lipids
- ffoplsaa_mod.itp      # Modified OPLSAA-FF, to account for lipid modifications
- ffoplsaabon_mod.itp   # Modified OPLSAA-FF(bonded), to account for lipid modifications
- ffoplsaanb_mod.itp    # Modified OPLSAA-FF(non-bonded), to account for lipid modifications
- topol.tpr             # Input for the first equilibration stage
- topol.top             # Topology of the system
- protein.itp           # Topology of the protein
- index.ndx             # Index file with appropriate groups for GROMACS
- prod_example.mdp      # Example of a parameter file to configure a production phase (see TIPS)

STRUCTURES:
- hexagon.pdb           # Initial structure of the system, with the receptor centered in the box 
- confout.gro           # Final structure of the system (see TIPS)

TRAJECTORY FILES
- traj_EQ.xtc           # Trajectory of the whole system in .xtc format: 1 snapshot/50 ps     
- ener_EQ.edr           # Energy file of the trajectory

REPORTS:
In the "reports" subfolder, you will find the following files:
- tot_ener.xvg, tot_ener.log    # System total energy plot and log
- temp.xvg, temp.log            # System temperature plot and log
- pressure.xvg, pressure.log    # System pressure plot and log
- volume.xvg, volume.log        # System volume plot and log

LOGS:
In the "logs" subfolder, you will find the log files of mdrun:
- eq_{force_constant}.log       # log of stages with restrained heavy atoms of the receptor
- eqCA.log                      # log of the stage with restrained C-alfa atoms of the receptor

* TIPS *
- If you want to configure a .tpr input file for production phase, you can use the template
'prod.mdp' file by introducing the number steps (nsteps), and thus the simulation time,
you want to run.
After that, you just have to type:

grompp -f prod.mdp -c confout.gro -p topol.top -n index.ndx -o topol_prod.tpr

- If you want to create a PDB file of your system after the equilibration, with the
receptor centered in the box, type:

echo 1 0 | trjconv -pbc mol -center -ur compact -f confout.gro -o confout.pdb

NOTE: these tips work for GROMACS version 4.0.5.

[1] Rodríguez D., Piñeiro Á. and Gutiérrez-de-Terán H.: Molecular Dynamics Simulations Reveal Insights into Key Structural
Elements of Adenosine Receptors (2011) Biochemistry. 50(19):4194-208

