import complex
import gromacs
import protein
import recipes
import membrane

import sys

#First we define all parts to be used
monomer = protein.Monomer(pdb = "Y1_min.pdb")
ligand = protein.Ligand(pdb = "lig.pdb", itp = "lig.itp")
membrane = membrane.Membrane()
g = gromacs.Gromacs()

#Now we create a complex membrane + protein(s) + ligand
prot_complex = protein.ProteinComplex(monomer = monomer, ligand = ligand)
full_complex = complex.MembraneComplex()
full_complex.complex = prot_complex
full_complex.membrane = membrane

#Now we call gromacs to make all the operations
g = gromacs.Gromacs(membrane_complex = full_complex)

#Now we can peek inside any object to look for its properties:
#print g.membrane_complex.box_height
#print g.membrane_complex.complex.monomer.pdb
#print g.membrane_complex.complex.ligand.pdb
#print g.membrane_complex.membrane.pdb
# print g.repo_dir
# ... and so on

print g.run_recipe()

sys.exit()

#Execute function by its name
#getattr(utils, step["command"])(step["options"])

