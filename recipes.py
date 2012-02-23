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
         {"command": "set_prot_size", #3
          "options": {"src": "proteinopls.pdb",
                      "dir": "xy"}},
         {"gromacs": "editconf", #4
          "options": {"src": "proteinopls.pdb",
                      "tgt": "proteinopls.pdb",
                      "dist": ""}},
         {"command": "set_prot_size", #5
          "options": {"src": "proteinopls.pdb",
                      "dir": "z"}},
         {"command": "set_popc", #6
          "options": {"tgt": "popc.pdb"}},
         {"command": "set_nanom"}, #7
         {"gromacs": "editconf", #8
          "options": {"src": "proteinopls.pdb",
                      "tgt": "proteinopls.pdb",
                      "box": "",
                      "angles": ["90", "90", "120"],
                      "bt": "tric"}},
         {"gromacs": "editconf", #9
          "options": {"src": "popc.pdb",
                      "tgt": "popc.pdb",
                      "box": ""}},
         {"command": "make_topol"}, #10
         {"gromacs": "editconf", #11
          "options": {"src": "proteinopls.pdb",
                      "tgt": "proteinopls.pdb",
                      "translate": ["0", "0", "0"]}},
         {"gromacs": "genbox", #12
          "options": {"cp": "proteinopls.pdb",
                      "cs": "popc.pdb",
                      "tgt": "protpopc.pdb",
                      "top": "topol.top"}},
          {"command": "make_water", #13
           "options": {"tgt": "water.pdb"}},
          {"gromacs": "editconf", #14
           "options": {"src": "water.pdb",
                       "tgt": "water.pdb",
                       "box": ""}},
          {"gromacs": "editconf", #15
           "options": {"src": "protpopc.pdb",
                       "tgt": "protpopc.pdb",
                       "box": "",
                       "angles": ["90", "90", "120"],
                       "bt": "tric"}},
          {"gromacs": "genbox", #16
           "options": {"cp": "protpopc.pdb",
                       "cs": "water.pdb",
                       "tgt": "tmp.pdb",
                       "top": "topol.top"}},
          {"command": "count_lipids"}, #17
          {"command": "make_topol"}, #18
          {"command": "make_topol_lipids"}, #19
          {"command": "make_steep"}, #20
          {"gromacs": "grompp", #21
           "options": {"src": "steep.mdp",
                       "src2": "tmp.pdb",
                       "tgt": "topol.tpr",
                       "top": "topol.top"}},
          {"gromacs": "trjconv", #22
           "options": {"src": "tmp.mdp",
                       "src2": "topol.tpr",
                       "tgt": "tmp.pdb",
                       "pbc": "mol"}},
          {"gromacs": "genion", #23
           "options": {"src": "topol.tpr",
                       "src2": "topol.top",
                       "tgt": "output.pdb",#<------------
                       "np": "",
                       "nn": ""}},
          {"gromacs": "grompp", #24
           "options": {"src": "steep.mdp",
                       "src2": "output.pdb",#<------------
                       "tgt": "topol.tpr",
                       "top": "topol.top"}},
          {"gromacs": "trjconv", #25
           "options": {"src": "output.pdb",#<------------
                       "src2": "topol.tpr",
                       "tgt": "output.pdb",#<------------
                       "trans": ["0", "0", "0"],
                       "pbc": "mol"}},
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
                       "pbc": "mol"}}
        ]
          
        self.breaks = [(0, "src", "membrane_complex.complex.monomer.pdb"),
                       (2, "dist", "membrane_complex.box_height"),
                       (4, "dist", "membrane_complex.box_width"),
                       

# This recipe modifies the previous one taking a ligand into account
class MonomerLigandRecipe(object):
    def __init__(self):
        super(MonomerLigandRecipe, self).__init__()
        self.recipe[18] = \
            {"command": "make_topol",
             "options": {"ligand": True}}
        self.recipe[10] = \
            {"command": "make_topol",
             "options": {"ligand": True}}

        self.recipe.insert(10,
            {"gromacs": "genrestr",
             "options": {"src": "",
                         "tgt": "posre_lig.itp",
                         "index": "lig_ha.ndx",
                         "forces": ["1000, 1000, 1000"],
                         "ligand": True}})

        self.recipe.insert(10,
            {"gromacs": "make_ndx",
             "options": {"src": "",
                         "tgt": "ligand_ha.ndx",
                         "ligand": True}})

        self.recipe.insert(2,
            {"command": "concat",
             "options": {"src": "proteinopls.pdb",
                         "tgt": "ligand.pdb"}})

        self.breaks = []#
