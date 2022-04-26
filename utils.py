import os
import shutil
from string import Template

import bw4posres


def clean_topol(src = [], tgt = []):
    """
    Clean the src topol of path specifics, and paste results in target
    """
    source = open(src, "r")
    target = open(tgt, "w")

    for line in source:
        newline = line
        if line.startswith("#include"):
            newline = line.split()[0] + ' "'
            newline += os.path.split(line.split()[1][1:-1])[1]
            newline += '"\n'
        target.write(newline)

    target.close()
    source.close()

    return True

def concat(**kwargs):
    """
    Make a whole pdb file with all the pdb provided
    """
    for compound_class in ["ligand", "ions", "cho", "alosteric", "waters"]:
        # Does the complex carry the group?
        if hasattr(kwargs["tgt"], compound_class):
            if getattr(kwargs["tgt"], compound_class):
                _file_append(kwargs["src"],
                             getattr(kwargs["tgt"], compound_class).pdb)

def getbw(**kwargs):
    """
    Call the Ballesteros-Weistein based pair-distance restraint
    module.
    """
    bw4posres.Run(kwargs["src"]).pdb2fas()
    bw4posres.Run(kwargs["src"]).clustalalign()
    bw4posres.Run(kwargs["src"]).getcalphas()
    bw4posres.Run(kwargs["src"]).makedisre()

def _file_append(f_src, f2a):
    """
    Add (concatenate) a f2a pdb file to another src pdb file
    """
    src = open(f_src, "r")
    f2a = open(f2a, "r")
    tgt = open("tmp_" + f_src, "w")

    for line in src:
        if ("TER" or "ENDMDL") not in line:
            tgt.write(line)
        else:
            for line_2_add in f2a:
                tgt.write(line_2_add)
            break
    tgt.write("TER\nENDMDL\n")
    tgt.close()
    f2a.close()
    src.close()

    shutil.copy(tgt.name, f_src)

    return True

def make_cat(dir1, dir2, name):
    """
    Very tight function to make a list of files to inject
    in some GROMACS suite programs
    """
    traj_src = [os.path.join(dir1, name)]
    traj_src.extend([os.path.join(dir1, "{0}", name).format(x)
                     for x in range(800, 0, -200)])
    if dir2 != "":
        traj_src.extend([os.path.join(dir2, name)])

    return traj_src

def make_ffoplsaanb(complex = None):
    """
    Join all OPLS force fields needed to run the simulation
    """
    ff = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                      "templates", "ffoplsaanb_")

    base = "{0}base.itp".format(ff) # This is the ff for proteins and other
    lip = "{0}lip.itp".format(ff)   # This for the lipids
    cho = "{0}cho.itp".format(ff)   # This for cholesterol

    to_concat = []
    if hasattr(complex, "ligand"):
        if hasattr(complex.ligand, "force_field"):
            to_concat.append(complex.ligand.force_field)
    if hasattr(complex, "alosteric"):
        if hasattr(complex.alosteric, "force_field"):
            to_concat.append(complex.alosteric.force_field)
    if hasattr(complex, "cho"):
        to_concat.append(cho)

    to_concat.extend([lip, base])

    output = "[ atomtypes ]\n"
    for ff_i in to_concat:
        output += open(ff_i).read()
        if not output.endswith("\n"): output += "\n"

    open("ffoplsaanb_mod.itp", "w").write(output)

    return True

def make_topol(template_dir = \
    os.path.join(os.path.dirname(os.path.realpath(__file__)), "templates"),
    target_dir = "",  # Dir where topol.top should land
    working_dir = "", # Dir where script is working
    complex = None):  # The MembraneComplex object to deal
    """
    Make the topol starting from our topol.top template
    """
    protein = bw =  dimer = lig = na = cho = alo = hoh = 0
    lig_name = ions_name = cho_name = alosteric_name = hoh_name = ""
    if hasattr(complex, "monomer"):
        protein = 1
        bw = 1
        if getattr(complex, "monomer").__class__.__name__ == "Dimer":
            dimer = 1
    if hasattr(complex, "ligand"):
        if complex.ligand:
            lig = 1
            lig_name = complex.ligand.itp
    if hasattr(complex, "ions"):
        if hasattr(complex.ions, "number"):
            na = complex.ions.number
            ions_name = complex.ions.itp
    if hasattr(complex, "cho"):
        if hasattr(complex.cho, "number"):
            cho = complex.cho.number
            cho_name = complex.cho.itp
    if hasattr(complex, "alosteric"):
        if complex.alosteric:
            alo = 1
            alosteric_name = complex.alosteric.itp
    if hasattr(complex, "waters"):
        if hasattr(complex.waters, "number"):
            hoh = complex.waters.number
            hoh_name = complex.waters.itp

    order = ("protein", "bw", "dimer", "lig", "na", "cho", "alo", "hoh")
    comps = {"protein": {"itp_name": "protein.itp",
                         "ifdef_name": "POSRES",
                         "posre_name": "posre.itp"},

             "bw": {"ifdef_name": "DISRE",
                    "posre_name": "disre.itp"},

             "dimer": {"name": "Protein_chain_B",
                 "itp_name": "protein_Protein_chain_B.itp",
                 "ifdef_name": "POSRES",
                 "posre_name": "posre_Protein_chain_B.itp"},

             "lig": {"itp_name": lig_name,
                 "ifdef_name": "POSRESLIG",
                 "posre_name": "posre_lig.itp"},

             "na": {"itp_name": ions_name,
                 "ifdef_name": "POSRESION",
                 "posre_name": "posre_ion.itp"},

             "cho": {"itp_name": cho_name,
                 "ifdef_name": "POSRESCHO",
                 "posre_name": "posre_cho.itp"},

             "alo": {"itp_name": alosteric_name,
                 "ifdef_name": "POSRESALO",
                 "posre_name": "posre_alo.itp"},

             "hoh": {"itp_name": hoh_name,
                 "ifdef_name": "POSRESHOH",
                 "posre_name": "posre_hoh.itp"},
             }

    if dimer:
        comps["protein"] = {"name": "Protein_chain_A",
            "itp_name": "protein_Protein_chain_A.itp",
            "ifdef_name": "POSRES",
            "posre_name": "posre_Protein_chain_A.itp"}

    src = open(os.path.join(template_dir, "topol.top"), "r")
    tgt = open(os.path.join(target_dir, "topol.top"), "w")

    t = Template("".join(src.readlines()))
    src.close()

    itp_include = []
    for c in order:
        if locals()[c]:
            if "itp_name" in comps[c].keys():
                posre_name = comps[c]["posre_name"]
                itp_include.append('#include "{0}"'.format(comps[c]["itp_name"]))
            if "posre_name" in comps[c].keys():
                itp_include.extend(['; Include Position restraint file',
                '#ifdef {0}'.format(comps[c]["ifdef_name"]),
                '#include "{0}"'.format(os.path.join(target_dir,
                    comps[c]["posre_name"])),
                '#endif'])

            if ("name") in comps[c]:
                comps[c]["line"] = "{0} {1}".format(
                    comps[c]["name"], locals()[c])
            else:
                comps[c]["line"] = "{0} {1}".format(c, locals()[c])
        else:
            comps[c]["line"] = ";"

    if working_dir: working_dir += "/" # Root dir doesn't need to be slashed

    tgt.write(t.substitute(working_dir = working_dir,
                           protein = comps["protein"]["line"],
                           dimer = comps["dimer"]["line"],
                           lig = comps["lig"]["line"],
                           na = comps["na"]["line"],
                           cho = comps["cho"]["line"],
                           alosteric = comps["alo"]["line"],
                           hoh = comps["hoh"]["line"],
                           itp_includes = "\n".join(itp_include)))
    tgt.close()

    return True

def tar_out(src_dir = [], tgt = []):
    """
    Tar everything in a src_dir to the tar_file
    """
    import tarfile

    t_f = tarfile.open(tgt, mode="w:gz")
    base_dir = os.getcwd()
    os.chdir(src_dir) # To avoid the include of all parent dirs
    for to_tar in os.listdir(os.path.join(base_dir, src_dir)):
        t_f.add(to_tar)
    t_f.close()
    os.chdir(base_dir)