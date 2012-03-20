import inspect
import os
import shutil
from string import Template

def _arrange_dir(src_dir, new_dir, useful_files=[], useful_fixed=[]):
    '''Copy the files in useful files from src_dir and
    fixed files from self.own_dir to new dir, which is created if needed'''

    if not os.path.isdir(new_dir):
        os.makedirs(new_dir)

    for u_f in [x for x in os.listdir(src_dir) if x in useful_files]:
        src = os.path.join(src_dir, u_f)
        tgt = os.path.join(new_dir, u_f)
        shutil.copy(src, tgt)

    for u_f in [x for x in os.listdir(self.own_dir) if x in useful_fixed]:
        src = os.path.join(self.own_dir, u_f)
        tgt = os.path.join(new_dir, u_f)

        shutil.copy(src, tgt)

    return True

def concat(**kwargs):
    '''Add (concatenate) a tgt pdb file to another src pdb file'''
    src = open(kwargs["src"], "r")
    f2a = open(kwargs["tgt"], "r")
    tgt_name = kwargs["src"].replace(".pdb", "-ligand.pdb")
    tgt = open(tgt_name, "w")

    for line in src:
        if ("TER" or "ENDMDL") not in line:
            tgt.write(line)
        else:
            for line_ligand in f2a:
                tgt.write(line_ligand)
            break
    tgt.write("TER\nENDMDL\n")
    tgt.close()
    f2a.close()
    src.close()

    shutil.copy(tgt_name,
                kwargs["src"])

    return True

def make_cat(dir1, dir2, name):
    '''Very tight function to make a list of files to inject 
    in some GROMACS suite programs
    '''
    traj_src = [os.path.join(dir1, name)]
    traj_src.extend([os.path.join(dir1, "{0}", name).format(x)
                     for x in range(800, 0, -200)])
    #traj_src.extend([os.path.join(dir2, name)])

    return traj_src

def make_topol(template_dir = "templates",
    target_dir = "", #Dir where topol.top should land
    working_dir = "", #Dir where script is working
    protein = "",
    ligand = "",
    crystal_waters = "",
    ions = ""):
    '''Make the topol starting from our topol.top template'''

    comps = {"protein": {"itp_name": "protein.itp",
                 "ifdef_name": "POSRES",
                 "posre_name": "posre.itp"},
             "ligand": {"itp_name": ligand.replace(".pdb", ".itp"),
                 "ifdef_name": "POSRESLIG",
                 "posre_name": "posre_lig.itp"},
             "crystal_waters": {"itp_name": "hoh.itp",
                 "ifdef_name": "POSRESHOH",
                 "posre_name": "posre_hoh.itp"},
             "ions": {"itp_name": "ions_local.itp",
                 "ifdef_name": "POSRESION",
                 "posre_name": "posre_ion.itp"}}

    src = open(os.path.join(template_dir, "topol.top"), "r")
    tgt = open(os.path.join(target_dir, "topol.top"), "w")

    t = Template("".join(src.readlines()))
    src.close()

    #for comp in comps.keys():
    #    if locals()[comp]:  
    #        return comp

    if(ligand):
        ligand_name = ligand.replace(".pdb", ".itp")
        lig_itp = '#include "{0}"\n'.format(ligand_name) + \
                  '; Include Position restraint file\n' + \
                  '#ifdef POSRESLIG\n' + \
                  '#include "posre_lig.itp"\n' + \
                  '#endif'
        lig_name = (os.path.join(working_dir, ligand_name + ".itp")) + " 1"
    else:
        lig_itp = ";"
        lig_name = ";"

    tgt.write(t.substitute(working_dir = working_dir,
                           protein_name = "protein 1",
                           lig_itp = lig_itp,
                           lig_name = "lig 1"))
    tgt.close()

    return True

def make_topol_lines(itp_name = "",
    ifdef_name = "",
    posre_name = ""):
    '''Make the topol lines to be included'''

    return "\n".join(['#include "{it}"',
        '; Include Position restraint file',
        '#ifdef {id}',
        '#include "{po}"',
        '#endif']).format(it = itp_name,
                          id = ifdef_name,
                          po = posre_name)

def tune_mdp(groups):
    '''Adjust the tc-groups of eq.mdp to be in line with our system'''
    shutil.move("Rmin/eq.mdp", "Rmin/eq.mdp~")
    eq = open("Rmin/eq.mdp~", "r")
    eq_out = open("Rmin/eq.mdp", "w")
    
    for line in eq:
        new_line = line
        if line.startswith("tc-grps"):
            new_line = line.replace("POP", groups["lipids"])
            new_line = line.replace("wation", groups["solvent"])
            new_line = line.replace("Protein", groups["complex"])
        eq_out.write(new_line)
    eq.close()
    eq_out.close()

    return True
    
