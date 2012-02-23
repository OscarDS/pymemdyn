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

def concat(src, file2add):
    '''Add (concatenate) a pdb file to another pdb'''
    src_f = open(src, "r")
    f2a = open(file2add, "r")
    tgt_name = src.replace(".pdb", "-ligand.pdb")
    tgt = open(tgt_name, "w")

    for line in src_f:
        if ("TER" or "ENDMDL") not in line:
            tgt.write(line)
        else:
            for line_ligand in f2a:
                tgt.write(line_ligand)
            break
    tgt.write("TER\nENDMDL\n")
    tgt.close()
    f2a.close()
    src_f.close()

    shutil.copy(tgt_name,
                src)

    return True

def count_lipids(src, tgt, prot_z = 0):
    '''Count the lipids in src and write a target with N4 tags'''
    src = open(src, "r")
    tgt = open(tgt, "w")

    half = prot_z / 2
    lipids_up = 0
    lipids_down = 0
    n_wats = 0

    for line in src:
        if len(line.split()) > 2:
            if line.split()[2] == "N4": #Lipid marker
                tgt.write(line)
                if(float(line.split()[7]) >= half):
                    lipids_up += 1
                elif(float(line.split()[7]) < half):
                    lipids_down += 1
            elif line.split()[2] == "OW": #Water marker
                n_wats += 1

    src.close()
    tgt.close()

    return {"lipids_up": lipids_up,
            "lipids_down": lipids_down,
            "n_wats": n_wats}

def make_topol(template_dir = "templates",
    target_dir = "", #Dir where topol.top should land
    working_dir = "", #Dir where script is working
    protein_name="",
    ligand_name = ""):
    '''Make the topol starting from our topol.top template'''

    src = open(os.path.join(template_dir, "topol.top"), "r")
    tgt = open(os.path.join(target_dir, "topol.top"), "w")

    t = Template("".join(src.readlines()))
    src.close()

    if(ligand_name):
        lig_itp = '#include "lig.itp"\n' + \
                  '; Include Position restraint file\n' + \
                  '#ifdef POSRESLIG\n' + \
                  '#include "posre_lig.itp"\n' + \
                  '#endif'
        lig_name = (os.path.join(working_dir, ligand_name + ".itp")) + " 1"
    else:
        lig_itp = ";"
        lig_name = ";"


    tgt.write(t.substitute(template_dir = template_dir,
                           working_dir = working_dir,
                           protein_name = protein_name,
                           lig_itp = lig_itp,
                           lig_name = lig_name))
    tgt.close()

    return True

def make_topol_lipids(topol_dir = "", lipids_up = 0, lipids_down = 0, n_wats = 0):
    '''Actualize topol.top with lipids and waters'''

    topol = open(os.path.join(topol_dir, "topol.top"), "a")
    topol.write("; Number of POPC molecules with higher z-coord value:\n")
    topol.write("POPC " + lipids_up + "\n")
    topol.write("; Number of POPC molecules with lower z-coord value:\n")
    topol.write("POPC " + lipids_down + "\n")
    topol.write("; Total number of water molecules:\n")
    topol.write("SOL " + n_wats + "\n")
    topol.close()

    return True

def _make_steep(template_dir = "", tgt_dir = ""):
    '''Copy the template steep.mdp to the target dir'''
    shutil.copy(os.path.join(template_dir, "steep.mdp"),
                os.path.join(tgt_dir, "steep.mdp"))

    return True

def _make_water(self,
    bilayer_zw = 0,
    prot_z = 0,
    x4pdb = "",
    tgt = ""):
    '''Create water from x4bilayer.pdb'''
    start = (bilayer_zw - prot_z) / 2
    end = start + prot_z

    src = open(x4pdb, "r")
    tgt = open(tgt, "w")

    res = "NULL"
    for line in src:
        if len(line.split()) > 7:
            if ((line.split()[2] == "OW") and
                ((float(line.split()[7]) > end) or
                 (float(line.split()[7]) < start))):
                res = line.split()[4]
            if ((line.split()[4] != res) and
                (line.split()[3] == "SOL")):
                tgt.write(line)

    tgt.close()
    src.close()

    return True
    return True

def set_itp(kwargs):
    '''Cut a top file to be usable later'''
    src = open(kwargs["src"], "r")
    tgt = open(kwargs["tgt"], "w")

    get_name = False

    for line in src:
        if line.startswith("#include"):
            pass
        elif line.startswith("; Include Position restraint file"):
            break
        else:
            tgt.write(line)

        if get_name and not line.startswith(";"):
            protein_name = line.split()[0]
            get_name = False

        if line.startswith("[ moleculetype ]"):
            get_name = True

    tgt.close()
    src.close()

    return protein_name
