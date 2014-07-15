
pyMEMdyn Version 1.1
================================================================================

pyMEMdyn  is a  standalone  python package  to  setup calculations  of
membrane molecular dynamics using the gromacs set of programs. The package
can be used either in a desktop environment, or in a cluster with popular 
queuing systems such as Torque or Slurm.


### Modeling Modules 

The following modules define the objects to be modelled.

- **protein.py**.  This module defines the ProteinComplex, Protein, Monomer,
Dimer, Compound, Ligand, CrystalWaters, Ions, Cholesterol, Lipids, 
and Alosteric objects. These  objects are  started with  the required files, 
and can then be passed  to other objects.   
- **membrane.py**. Defines the cellular membrane.  
- **complex.py**.  Defines the full complex, protein + membrane.   
  It can  include any  of the previous objects.

### Auxiliary Modules

- **queue.py**.   Queue  manager.  That  is,  it  receives  objects to  be
  executed.   
- **recipes.py**.   Applies  step by  step instructions  for  carrying a 
  modeling  step.  
- **utils.py**.  Puts the  functions done by the previous objects on demand.
  For example, manipulate files, copy  folders, etc.

### Execution Modules

- **gromacs.py**. Defines the Gromacs and Wrapper objects.  * Gromacs will
  load the  objects to be modeled,  the modeling recipe, and  run it.  *
  Wrapper is a  proxy for gromacs commands. When a  recipe entry is sent
  to it this returns the command to be run.

### Examples

- **example.py**  An example  showing how  to use  the  previously defined
  libraries.


Changelog
--------------------------------------------------------------------------------

### Changes from version 1.0 to 1.1

- Wednesday, July 9, 2014

Among many changes to get pyMEMdyn up an running with gromacs 4.6.5 instead of
4.0.5 a new substitution for HIE, HID, and HIP is done. Previously the 
substitution was HIE:HISB, HID:HISA, HIP:HISH, now it's HIE:HISE, HID:HISD,
HIP:HISH



