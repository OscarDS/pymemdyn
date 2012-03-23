import inspect
import os
import shutil
from string import Template

import sys

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
    protein = "", # Number of proteins
    lig = "", # Number of ligands
    hoh = "", # Number of crystal waters
    na = ""): # Number of ions
    '''Make the topol starting from our topol.top template'''

    order = ("protein", "hoh", "lig", "na")
    comps = {"protein": {"itp_name": "protein.itp",
                 "ifdef_name": "POSRES",
                 "posre_name": "posre.itp"},
             "lig": {"itp_name": "ligand.itp",
                 "ifdef_name": "POSRESLIG",
                 "posre_name": "posre_lig.itp"},
             "hoh": {"itp_name": "hoh.itp",
                 "ifdef_name": "POSRESHOH",
                 "posre_name": "posre_hoh.itp"},
             "na": {"itp_name": "ions_local.itp",
                 "ifdef_name": "POSRESION",
                 "posre_name": "posre_ion.itp"}}

    src = open(os.path.join(template_dir, "topol.top"), "r")
    tgt = open(os.path.join(target_dir, "topol.top"), "w")

    t = Template("".join(src.readlines()))
    src.close()

    itp_include = []
    for c in order:
        if locals()[c]:
            itp_name = comps[c]["itp_name"]
            itp_include.extend(['#include "{0}"'.format(comps[c]["itp_name"]),
                '; Include Position restraint file',
                '#ifdef {0}'.format(comps[c]["ifdef_name"]),
                '#include "{0}/{1}"'.format(working_dir,
                                           comps[c]["posre_name"]),
                '#endif'])

            comps[c]["line"] = "{0} {1}".format(c, locals()[c])
        else:
            comps[c]["line"] = ";"

    if working_dir: working_dir += "/" #Root dir doesn't need to be slashed

    tgt.write(t.substitute(working_dir = working_dir,
                           protein = comps["protein"]["line"],
                           lig = comps["lig"]["line"],
                           hoh = comps["hoh"]["line"],
                           na = comps["na"]["line"],
                           itp_includes = "\n".join(itp_include)))
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
    
