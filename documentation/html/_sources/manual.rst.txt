Manual
======

The fully automated pipeline available by using **PyMemDyn** allows any
researcher, without prior experience in computational chemistry, to
perform an otherwise tedious and complex process of membrane insertion
and thorough MD equilibration, as outlined in Figure 1.

.. image:: pymemdyn.png

In the simplest scenario, only the receptor structure is considered. In
such case the GPCR is automatically surrounded by a pre-equilibrated
POPC (Palmitoyl-Oleoyl-Phosphatidyl-Choline) membrane model in a way
that the TM (Trans-Membrane) bundle is parallel to the vertical axis of
the membrane. The system is then soaked with bulk water and inserted
into an hexagonal prism-shaped box, which is energy-minimized and
carefully equilibrated in the framework of periodic boundary conditions
(PBC). A thorough MD equilibration protocol lasting 2.5 ns follows.

But the simulation of an isolated receptor can only account for one part
of the problem, and the influence of different non-protein elements in
receptor dynamics such as the orthosteric (primary) ligand, allosteric
modulator, or even specific cholesterol, lipid, water or ion molecules
are key for a more comprehensive characterization of GPCRs. **PyMemDyn**
can explicitly handle these elements allowing a broader audience in the
field of GPCRs to use molecular dynamics simulations. These molecules
should be uploaded in the same way they're present in the original PDB
file of the receptor, so they are properly integrated into the membrane
insertion protocol described above, together with force-field associated
files (which can be either generated with external software like
Macromodel, or by manual parameterization). In addition, it is also
possible to perform MD simulations of receptor dimers, provided that a
proper dimerization model exists (i.e., coming from X-ray
crystallography or from a protein-protein docking protocol). The ease of
use, flexibility and public availability of the **PyMemDyn** library
makes it a unique tool for researchers in the GPCR field interested in
exploring dynamic processes of these receptors.


Running with queues
-------------------

$\approx$ 90% of the time you will want to use some queueing system. We
deal with queue systems tweaks as we stumble into them and it's out of
our scope to cover them all. If you take a look at the source code dir,
you'll found some files called "run_pbs.sh", "run_svgd.sh" and so on.
Also there are specific queue objects in the source file queue.py we
have to tweak for every and each queue. In you want to run your
simulation in a supported queue, copy the "run_queuename.sh" file to
your working directory, and edit it. E.g. the workdir to run an A2a.pdb
simulation in svgd.cesga.es looks like: . .. A2a.pdb run_svgd.sh And
run_svgd.sh looks like:

::

   \$!/bin/bash   
   module   load  python/3.7
   module   load  gromacs/2021
   python ~/bin/pymoldyn/pymemdyn  -p  a2a.pdb 

Now we just launch this script with:

::

   qsub -l arch=amd,num_proc=1,s_rt=50:00:00,s_vmem=1G,h_fsize=1G -pe mpi 8 run_svgd.sh 

and wait for the results. Note that we launch 1 process, but flag the
run as mpi with reservation of 8 cores in SVGD queue.

Debugging 
---------

If you are to set up a new system, it is a good idea to just run a few
steps of each stage in the equilibration protocol just to test that the
pdb file is read correctly and the membrane-insertion protocol works
fine and the system can be minimized and does not "explode" during the
equilibration protocol (i.e., detect if atom clashes and so on exist on
your system).

To do this, use the --debug option, like:

::

   pymemdyn -p gpcr.pdb -l lig --waters hoh --debug 

If everythings works fine, you will see the list of output directories
and files just as in a regular equilibration protocol, but with much
smaller files (since we only use here 1000 steps of MD in each stage).
NOTE that sometimes, due to the need of a smooth equilibration procedure
(i.e. when a new ligand is introduced in the binding site without
further refinement of the complex, or with slight clashes of existing
water molecules) this kind of debugging equilibration procedure might
crash during the first stages due to hot atoms or LINCS failure. This is
normal, and you have two options: i) trust that the full equilibration
procedure will fix the steric clashes in your starting system, and then
directly run the pymoldyn without the debugging option, or ii) identify
the hot atoms (check the mdrun.log file in the last subdirectory that
was written in your output (generally eq/mdrun.log and look for "LINCS
WARNING"). What if you want to check partial functions of pymoldyn? In
order to do this you must edit the file pymemdyn and change:

1. Line 260 comment with "#" this line [that states: run.clean()], which
   is the one that deletes all the output files present in the working
   directory.

2. In the last two lines of this file, comment (add a "#") the line:
   run.moldyn()

3. And uncomment (remove the "#") the line: run.light_moldyn()

4. In the line 153 and within that block (ligh_moldyn) change the lines
   stating steps = ["xxxx"] and include only those steps that you want
   to test, which should be within a list of strings.

For the sake of clarity, these have been subdivided in two lines:

::

   line 1-  steps = ["Init", "Minimization", "Equilibration", "Relax", "CARelax"] 

Here you remove those strings that you do not want to be executed, i.e.
if only membrane insertion and minimization is wished, remove
"Equilibration", "Relax", "CARelax" so the line states:

::

   steps   =   ["Init",  "Minimization"]  

   line  2-   steps  =  ["CollectResults"]  

This only accounts for the preparation of the output files for analysis,
so if you only are interested on this stage, comment the previous line.
The last assignment is the one that runs. NOTE that you must know what
you do, otherwise you might have crashes in the code if needed files to
run intermediate stages are missing!

Output 
------

The performed equilibration includes the following stages:


.. list-table:: Output
   :widths: 20 35 25 20
   :header-rows: 1

   * - STAGE
     - RESTRAINED ATOMS
     - FORCE CONSTANT
     - TIME
   * - 
     -
     - kJ/(mol A ·nm^2) 
     - ns
   * - Minimization 
     -                       
     -              
     - (Max. 500 steps)
   * - Equil. 1     
     - Protein Heavy Atoms     
     - 1000           
     - 0.5
   * - Equil. 2     
     - Protein Heavy Atoms     
     - 800            
     - 0.5
   * - Equil. 3     
     - Protein Heavy Atoms     
     - 600            
     - 0.5
   * - Equil. 4     
     - Protein Heavy Atoms     
     - 400            
     - 0.5
   * - Equil. 5     
     - Protein Heavy Atoms     
     - 200            
     - 0.5 
   * - Equil. 6 
     - Venkatakrishnan Pairs / 
     - 200 /          
     - 2.5 /
   * - Equil. 6 
     - C-alpha Atoms           
     - 200            
     - 2.5


In this folder you will find several files related to this simulation:

INPUT:
~~~~~~

::

   - popc.itp              # Topology of the lipids  
   - ffoplsaa_mod.itp      # Modified OPLSAA-FF, to account for lipid modifications  
   - ffoplsaabon_mod.itp   # Modified OPLSAA-FF(bonded), to account for lipid modifications   
   - ffoplsaanb_mod.itp    # Modified OPLSAA-FF(non-bonded), to account for lipid modifications
   - topol.tpr             # Input for the first equilibration stage
   - topol.top             # Topology of the system
   - protein.itp           # Topology of the protein
   - index.ndx             # Index file with appropriate groups for GROMACS
   - prod.mdp              # Example of a parameter file to configure a production run (see TIPS)

STRUCTURES:
~~~~~~~~~~~

::

   - hexagon.pdb           # Initial structure of the system, with the receptor centered in the box 
   - confout.gro           # Final structure of the system (see TIPS)
   - load_gpcr.pml         # Loads the initial structure and the trajectory in pymol

TRAJECTORY FILES:
~~~~~~~~~~~~~~~~~

::

   - traj_pymol.xtc        # Trajectory of the whole system for visualization in pymol. 1 snapshot/100 ps
   - traj_EQ.xtc           # Trajectory of the whole system in .xtc format: 1 snapshot/50 ps 
   - ener_EQ.edr           # Energy file of the trajectory
   - load_gpcr.pml         # Script to load the equilibration trajectory in pymol.

REPORTS:
~~~~~~~~

| In the "reports" subfolder, you will find the following files:
| - tot_ener.xvg, tot_ener.log # System total energy plot and log -
  temp.xvg, temp.log # System temperature plot and log - pressure.xvg,
  pressure.log # System pressure plot and log - volume.xvg, volume.log #
  System volume plot and log - rmsd-all-atom-vs-start # All atoms RMSD
  plot - rmsd-backbone-vs-start.xvg # Backbone RMSD plot -
  rmsd-calpha-vs-start.xvg # C-Alpha RMSD plot - rmsf-per-residue.xvg #
  Residue RMSF plot

LOGS:
-----

| In the "logs" subfolder, you will find the log files of mdrun:
| - eq_{force_constant}.log # log of stages with restrained heavy atoms
  of the receptor - eqCA.log # log of the stage with restrained C-alfa
  atoms of the receptor

**NOTE ON GROMACS METHODS** To integrate the equations of motion we have
selected the leap-frog integrator with a 2 femtosecond timestep.
Longe-range electrostatic interactions in periodic boundary conditions
are treated with the particle mesh Ewald method. We use a Nose-Hoover
thermostat with a tau_t of 0.5 picoseconds and a Parinello-Rahman
barostat with a tau_p of 2.0. The pressure coupling is semiisotropic,
meaning that it's isotropic in the x and y directions but different in
the z direction. Since we are using pressure coupling we are working
with an NPT ensemble. This is done both in the all-atom restrained steps
and in the alpha-carbon atom restrained part. All of these details are
more explicitly stated in the Rodriguez et al. [1] publication.

**TIPS**

NOTE: these tips work for GROMACS version >= 4.5 and < 5.0. For later
versions, adjustments are required, but the principle remains the same.

-  If you want to configure a .tpr input file for a **production** run,
   you can use the template 'prod.mdp' file by introducing the number of
   steps (nsteps), and thus the simulation time, you want to run.

After that, you just have to type:

::

   grompp -f prod.mdp -c confout.gro -p topol.top -n index.ndx -o topol_prod.tpr  
   mdrun -s topol_prod.tpr -o traj.trr -e ener.edr -c confout.gro -g production.log -x traj_prod.xtc

-  If you want to create a PDB file of your system after the
   equilibration, with the receptor centered in the box, type:

   echo 1 0 \| trjconv -pbc mol -center -ur compact -f confout.gro -o
   confout.pdb

-  If you want to create an xmgrace graph of the root mean square
   deviation for c-alpha atoms in the 5.0 ns of simulation you can use:

   echo 3 3 \| g_rms -f traj_EQ.xtc -s topol.tpr -o
   rmsd-calpha-vs-start.xvg

-  You may want to get a pdb file of your last frame. You can first
   check the total time of your trajectory and then use this time to
   request the last frame with:

   gmxcheck -f traj_pymol.xtc echo 1 \| trjconv -b 5000 -e 5000 -f
   traj_pymol.xtc -o last51.pdb

References
----------

| [1] Rodríguez D., Piñeiro A. and Gutiérrez-de-Terán H.
| Molecular Dynamics Simulations Reveal Insights into Key Structural
  Elements of Adenosine Receptors
| Biochemistry (2011), 50, 4194-208.
