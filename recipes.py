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
                       "tgt": "output.pdb",#<------------
                       "np": "",
                       "nn": ""},
           "input": "SOL\n"},
          {"gromacs": "grompp", #24
           "options": {"src": "steep.mdp",
                       "src2": "output.pdb",#<------------
                       "tgt": "topol.tpr",
                       "top": "topol.top"}},
          {"gromacs": "trjconv", #25
           "options": {"src": "output.pdb",#<------------
                       "src2": "topol.tpr",
                       "tgt": "output.pdb",#<------------
                       "trans": [],
                       "pbc": "mol"},
           "input": "0\n"},
          {"gromacs": "grompp", #26
           "options": {"src": "steep.mdp",
                       "src2": "output.pdb",
                       "tgt": "output.tpr", #<------------
                       "top": "topol.top"}},
          {"gromacs": "trjconv", #27
           "options": {"src": "output.pdb",
                       "src2": "output.tpr",
                       "tgt": "hexagon.pdb",
                       "ur": "compact",
                       "pbc": "mol"},
           "input": "1\n0\n"}
        ]

        self.minimize = [{"command": "set_init", #1
                          "options": {"src": "topol.tpr",
                                      "mdp": "eq.mdp"}},
                         {"gromacs": "md_run", #0
                          "options": {"src": "topol.tpr",
                                      "tgt": "traj.trj",
                                      "energy": "ener.edr",
                                      "conf": "confout.gro",
				      "log": "md.log"}},

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
                         "index": "lig_ha.ndx",
                         "forces": ["1000, 1000, 1000"]},
             "input": "1\n"})

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
                       28: {"trans": "membrane_complex.complex.trans"}
                      }

class BasicMinimization(object):
    def __init__(self):
        self.recipe = \
        [{"command": "set_minimization_init", #1
          "options": {"src": "topol.tpr",
                      "mdp": "eq.mdp"}},
         {"gromacs": "md_run", #0
          "options": {"src": "topol.tpr",
                      "tgt": "traj.trj",
                      "energy": "ener.edr",
                      "conf": "confout.gro",
                      "log": "md.log"}},
        ]

        self.breaks = {}

