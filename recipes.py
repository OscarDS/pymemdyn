import os

class MonomerRecipe(object):
    def __init__(self):
        self.recipe = \
        [{"gromacs": "pdb2gmx", #0
          "options": {"src": "",
                      "tgt": "proteinopls.pdb",
                      "top": "protein.top"}},
         {"command": "set_itp", #1
          "options": {"src": "protein.top",
                      "tgt": "protein.itp"}},
         {"gromacs": "editconf", #2
          "options": {"src": "proteinopls.pdb",
                      "tgt": "proteinopls.pdb",
                      "dist": ""}},
         {"command": "set_protein_size", #3
          "options": {"src": "proteinopls.pdb",
                      "dir": "xy"}},
         {"gromacs": "editconf", #4
          "options": {"src": "proteinopls.pdb",
                      "tgt": "proteinopls.pdb",
                      "dist": ""}},
         {"command": "set_protein_size", #5
          "options": {"src": "proteinopls.pdb",
                      "dir": "z"}},
         {"command": "set_popc", #6
          "options": {"tgt": "popc.pdb"}},
         {"gromacs": "editconf", #7
          "options": {"src": "proteinopls.pdb",
                      "tgt": "proteinopls.pdb",
                      "box": "",
                      "angles": ["90", "90", "120"],
                      "bt": "tric"}},
         {"gromacs": "editconf", #8
          "options": {"src": "popc.pdb",
                      "tgt": "popc.pdb",
                      "box": ""}},
         {"command": "make_topol"}, #9
         {"gromacs": "editconf", #10
          "options": {"src": "proteinopls.pdb",
                      "tgt": "proteinopls.pdb",
                      "translate": ["0", "0", "0"]}},
         {"gromacs": "genbox", #11
          "options": {"cp": "proteinopls.pdb",
                      "cs": "popc.pdb",
                      "tgt": "protpopc.pdb",
                      "top": "topol.top"}},
          {"command": "set_water", #12
           "options": {"tgt": "water.pdb"}},
          {"gromacs": "editconf", #13
           "options": {"src": "water.pdb",
                       "tgt": "water.pdb",
                       "box": ""}},
          {"gromacs": "editconf", #14
           "options": {"src": "protpopc.pdb",
                       "tgt": "protpopc.pdb",
                       "box": "",
                       "angles": ["90", "90", "120"],
                       "bt": "tric"}},
          {"gromacs": "genbox", #15
           "options": {"cp": "protpopc.pdb",
                       "cs": "water.pdb",
                       "tgt": "tmp.pdb",
                       "top": "topol.top"}},
          {"command": "count_lipids", #16
           "options": {"src": "tmp.pdb",
                       "tgt": "popc.pdb"}},
          {"command": "make_topol"}, #17
          {"command": "make_topol_lipids"}, #18
          {"command": "set_grompp"}, #19
          {"gromacs": "grompp", #20
           "options": {"src": "steep.mdp",
                       "src2": "tmp.pdb",
                       "tgt": "topol.tpr",
                       "top": "topol.top"}},
          {"gromacs": "trjconv", #21
           "options": {"src": "tmp.pdb",
                       "src2": "topol.tpr",
                       "tgt": "tmp.pdb",
                       "pbc": "mol"},
           "input": "1\n0\n"},
          {"command": "get_charge",
           "options": {"src": "steep.mdp",
                       "src2": "tmp.pdb",
                       "tgt": "topol.tpr",
                       "top": "topol.top"}}, #22
          {"gromacs": "genion", #23
           "options": {"src": "topol.tpr",
                       "src2": "topol.top",
                       "tgt": "output.pdb",
                       "np": "",
                       "nn": ""},
           "input": "SOL\n"},
          {"gromacs": "grompp", #24
           "options": {"src": "steep.mdp",
                       "src2": "output.pdb",
                       "tgt": "topol.tpr",
                       "top": "topol.top"}},
          {"gromacs": "trjconv", #25
           "options": {"src": "output.pdb",
                       "src2": "topol.tpr",
                       "tgt": "output.pdb",
                       "trans": [],
                       "pbc": "mol"},
           "input": "0\n"},
          {"gromacs": "grompp", #26
           "options": {"src": "steep.mdp",
                       "src2": "output.pdb",
                       "tgt": "topol.tpr",
                       "top": "topol.top"}},
          {"gromacs": "trjconv", #27
           "options": {"src": "output.pdb",
                       "src2": "topol.tpr",
                       "tgt": "hexagon.pdb",
                       "ur": "compact",
                       "pbc": "mol"},
           "input": "1\n0\n"}
        ]

        self.breaks = {0: {"src": "membrane_complex.complex.monomer.pdb"},
                       2: {"dist": "membrane_complex.box_height"},
                       4: {"dist": "membrane_complex.box_width"},
                       7: {"box": "membrane_complex.trans_box_size"},
                       8: {"box": "membrane_complex.bilayer_box_size"},
                       13: {"box": "membrane_complex.embeded_box_size"},
                       14: {"box": "membrane_complex.protein_box_size"},
                      }
                       

# This recipe modifies the previous one taking a ligand into account
class MonomerLigandRecipe(MonomerRecipe):
    def __init__(self):
        super(MonomerLigandRecipe, self).__init__()
        self.recipe[17] = \
            {"command": "make_topol",
             "options": {"ligand": ""}}
        self.recipe[9] = \
            {"command": "make_topol",
             "options": {"ligand": ""}}

        self.recipe.insert(9,
            {"gromacs": "genrestr",
             "options": {"src": "",
                         "tgt": "posre_lig.itp",
                         "index": "ligand_ha.ndx",
                         "forces": ["1000", "1000", "1000"]},
             "input": "2\n"})

        self.recipe.insert(9,
            {"gromacs": "make_ndx",
             "options": {"src": "",
                         "tgt": "ligand_ha.ndx",
                         "ligand": True},
             "input": "! a H*\nq\n"})

        self.recipe.insert(2,
            {"command": "concat",
             "options": {"src": "proteinopls.pdb",
                         "tgt": ""}})

        self.breaks = {0: {"src": "membrane_complex.complex.monomer.pdb"},
                       2: {"tgt": "membrane_complex.complex.ligand.pdb"},
                       3: {"dist": "membrane_complex.box_height"},
                       5: {"dist": "membrane_complex.box_width"},
                       8: {"box": "membrane_complex.trans_box_size"},
                       9: {"box": "membrane_complex.bilayer_box_size"},
                       10: {"src": "membrane_complex.complex.ligand.pdb"},
                       11: {"src": "membrane_complex.complex.ligand.pdb"},
                       12: {"ligand": "membrane_complex.complex.ligand.pdb",
                            "protein": "membrane_complex.complex.monomer.pdb"},
                       #13: {"box": "membrane_complex.embeded_box_size"},
                       14: {"box": "membrane_complex.protein_box_size"},
                       16: {"box": "membrane_complex.embeded_box_size"},
                       17: {"box": "membrane_complex.protein_box_size"},
                       20: {"ligand": "membrane_complex.complex.ligand.pdb",
                            "protein": "membrane_complex.complex.monomer.pdb"},
                       26: {"np": "membrane_complex.complex.positive_charge",
                            "nn": "membrane_complex.complex.negative_charge"},
                       28: {"trans": "membrane_complex.complex.trans"},
                      }

class BasicMinimization(object):
    def __init__(self):
        self.recipe = \
        [{"command": "set_minimization_init", #0
          "options": {"src_tpr": "topol.tpr",
                      "min_dir": "Rmin",
                      "mdp": "eq.mdp"}},
         {"gromacs": "mdrun", #1
          "options": {"src": "Rmin/topol.tpr",
                      "tgt": "Rmin/traj.trj",
                      "energy": "Rmin/ener.edr",
                      "conf": "Rmin/confout.gro",
                      "log": "Rmin/md.log"}},
        ]
        self.breaks = {}

class BasicEquilibration(object):
    def __init__(self):
        self.recipe = \
        [{"gromacs": "editconf", #0
          "options": {"src": "Rmin/confout.gro",
                      "tgt": "min.pdb"}},
         {"command": "make_ndx", #1
          "options": {"src": "min.pdb",
                      "tgt": "index.ndx"}},
         {"gromacs": "grompp", #2
          "options": {"src": "Rmin/eq.mdp",
                      "src2": "min.pdb",
                      "top": "topol.top",
                      "tgt": "topol.tpr",
                      "index":"index.ndx"}},
         {"command": "set_equilibration_init", #3
          "options": {"src_tpr": "Rmin/topol.tpr",
                      "src_itp": "posre.itp",
                      "mdp": "Rmin/eq.mdp",
                      "eq_dir": "eq"}},
         {"gromacs": "mdrun", #4
          "options": {"src": "eq/topol.tpr",
                      "tgt": "eq/traj.trj",
                      "energy": "eq/ener.edr",
                      "conf": "eq/confout.gro",
                      "traj": "eq/traj.xtc",
                      "log": "eq/md_eq1000.log"}},
        ]
        self.breaks = {}

class LigandEquilibration(BasicEquilibration):
    def __init__(self):
        super(LigandEquilibration, self).__init__()
        self.recipe[0]["options"]["src_lig_itp"] = "posre_lig.itp"
        self.recipe.insert(2,
            {"gromacs": "genrestr",
             "options": {"src": "Rmin/topol.tpr",
                         "tgt": "protein_ca200.itp",
                         "index": "index.ndx",
                         "forces": ["200", "200", "200"]},
             "input": "3\n"})

class BasicRelax(object):
    def __init__(self):
        self.recipe = []
        for const in range(800, 0, -200):
            tgt_dir = "eq/{0}".format(const)
            src_dir = "eq"
            self.recipe += [
            {"command": "relax", #0, 3, 6, 9
             "options": {"const": const,
                         "src_dir": src_dir,
                         "tgt_dir": tgt_dir}},
            {"gromacs": "grompp", #1, 4, 7, 10
             "options": {"src": os.path.join(tgt_dir, "eq.mdp"),
                         "src2": os.path.join(src_dir, "confout.gro"),
                         "top": os.path.join(tgt_dir, "topol.top"),
                         "tgt": os.path.join(tgt_dir, "topol.tpr"),
                         "index": "index.ndx"}},
            {"gromacs": "mdrun", #2, 5, 8, 11
             "options": {"src": os.path.join(tgt_dir, "topol.tpr"),
                         "tgt": os.path.join(tgt_dir, "traj.trr"),
                         "energy": os.path.join(tgt_dir, "ener.edr"),
                         "conf": os.path.join(src_dir, "confout.gro"),
                         "traj": os.path.join(tgt_dir, "traj.xtc"),
                         "log": os.path.join(src_dir,
                             "md_eq{0}.log".format(const))}},
            ]
        self.breaks = {}

class CAEquilibrate(object):
    def __init__(self):
        self.recipe = \
            [{"command": "set_caequil_init", #0
              "options": {"src_dir": "eq",
                          "tgt_dir": "eqCA",
                          "src_files": ["confout.gro", "eq.mdp"],
                          "tmp_files": ["eqCA.mdp"]}},
             {"gromacs": "genrestr", #1
              "options": {"src": "Rmin/topol.tpr",
                          "tgt": "posre.itp",
                          "index": "index.ndx",
                          "forces": ["200"] * 3},
              "input": "3\n"},
             {"gromacs": "grompp", #2
              "options": {"src": "eqCA/eq.mdp",
                          "src2": "eqCA/confout.gro",
                          "top": "topol.top",
                          "tgt": "eqCA/topol.tpr",
                          "index": "index.ndx"}},
             {"gromacs": "mdrun", #3
              "options": {"src": "eqCA/topol.tpr",
                          "tgt": "eqCA/traj.trr",
                          "energy": "eqCA/ener.edr",
                          "conf": "eqCA/confout.gro",
                          "traj": "eqCA/traj.xtc",
                          "log": "eqCA/md_eqCA.log"}},
             {"gromacs": "trjcat", #4
              "options": {"dir1": "eq",
                          "dir2": "",
                          "name": "traj.xtc",
                          "tgt": "traj_EQ.xtc"},
              "input": "c\n" * 6},
             {"gromacs": "eneconv", #5
              "options": {"dir1": "eq",
                          "dir2": "eqCA",
                          "name": "ener.edr",
                          "tgt": "ener_EQ.edr"},
              "input": "c\n" * 6}
             ]

        self.breaks = {}
