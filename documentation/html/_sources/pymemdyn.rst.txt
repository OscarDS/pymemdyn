Pymemdyn 
========

This is the main script for the pymemdyn commandline tool.
In this script the following things are accomplished:

        1. Command line arguments are parsed.
        2. (If necessary) a working directory is created.
        3. (If necessary) previous Run files are removed.
        4. A run is done. 


Usage
-----

.. code-block:: bash

   usage: pymemdyn [-h] [-v] [-b OWN_DIR] [-r REPO_DIR] -p PDB [-l LIGAND]
                [-a ALLOSTERIC] [-w WATERS] [-i IONS] [-c CHO]
                [--res RESTRAINT] [--llc LIGPARGEN_LIGAND_CHARGE]
                [--llo LIGPARGEN_LIGAND_NROFOPTIMIZATIONS]
                [--lac LIGPARGEN_ALLOSTERIC_CHARGE]
                [--lao LIGPARGEN_ALLOSTERIC_NROFOPTIMIZATIONS] [-q QUEUE] [-d]
   
   == setup molecular dynamics for membrane proteins given a pdb. ==
   
   optional arguments:
     -h, --help            show this help message and exit
     -v, --version         show program's version number and exit
     -b OWN_DIR            Working dir if different from actual dir
     -r REPO_DIR           Path to templates of fixed files. If not provided,
                           take the value from settings.TEMPLATES_DIR.
     -p PDB                Name of the pdb to insert into membrane for MD
                           (mandatory). Use the pdb extension. (e.g. -p
                           myprot.pdb)
     -l LIGAND, --lig LIGAND
                           Name of the ligand, without extension. See
                           input_guide.txt for details on how to generate the
                           required pdb and forcefield files.
     -a ALLOSTERIC, --alo ALLOSTERIC
                           Name of the allosteric, without extension. See
                           input_guide.txt for details on how to generate the
                           required pdb and forcefield files.
     -w WATERS, --waters WATERS
                           Crystalized water molecules. File name without
                           extension.
     -i IONS, --ions IONS  Crystalized ions file name without extension.
     -c CHO, --cho CHO     Crystalized cholesterol molecules file name without
                           extension.
     --res RESTRAINT       Position restraints during MD production run. Options:
                           bw (Ballesteros-Weinstein Restrained Relaxation -
                           default), ca (C-Alpha Restrained Relaxation)
     --llc LIGPARGEN_LIGAND_CHARGE
                           Charge of ligand for ligpargen (when itp file should
                           be generated)
     --llo LIGPARGEN_LIGAND_NROFOPTIMIZATIONS
                           Number of optimizations that ligpargen should use to
                           generate itp file for ligand (only needed when itp is
                           not provided)
     --lac LIGPARGEN_ALLOSTERIC_CHARGE
                           Charge of allosteric for ligpargen (when itp file
                           should be generated)
     --lao LIGPARGEN_ALLOSTERIC_NROFOPTIMIZATIONS
                           Number of optimizations that ligpargen should use to
                           generate itp file for allosteric (only needed when itp
                           is not provided)
     -q QUEUE, --queue QUEUE
                           Queueing system to use (slurm, pbs, pbs_ib and svgd
                           supported)
     -d, --debug

