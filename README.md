
pyMEMdyn Version 1.1
================================================================================

pyMEMdyn  is a  standalone  python library  to  setup calculations  of
membrane molecular dynamics with the gromacs set of programs either in
a desktop environment or in a cluster using popular queuing systems 
such as Torque or Slurm.


### Modeling Objects

The following three files define the objects to be modelled.

- protein.py. Includes  the definition of the  Monomer, Ligand, Dimer,
  and  ProteinComplex objects.   These  objects are  started with  the
  required  files,   and  then  are   passed  to  other   objects.   
- membrane.py.  Defines the cellular membrane.  
- complex.py.  Defines   the  full object, protein  + membrane.   
  It can  include any  of the previous objects.

### Auxiliary Objects

- queue.py.   Queue  manager.  That  is,  it  receives  objects to  be
  executed.   
- recipes.py.   Applies  step by  step instructions  for  carrying a 
  modeling  step.  
- utils.py.  Puts the  functions done by the previous objects on demand.
  For example, manipulate files, copy  folders, etc.

### Execution Objects

- gromacs.py. Defines the Gromacs and Wrapper objects.  * Gromacs will
  load the  objects to be modeled,  the modeling recipe, and  run it.  *
  Wrapper is a  proxy for gromacs commands. When a  recipe entry is sent
  to it returns the command to be run.

### Examples

- example.py  An example  showing how  to use  the  previously defined
  libraries.


Changelog
--------------------------------------------------------------------------------

- Wednesday, July 9, 2014

Among many changes to get pyMEMdyn up an running with gromacs 4.6.5 instead of
4.0.5 a new substitution for HIE, HID, and HIP is done. Previously the 
substitution was HIE:HISB, HID:HISA, HIP:HISH, now it's HIE:HISE, HID:HISD,
HIP:HISH



