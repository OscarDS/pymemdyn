import os

class MonomerRecipe(object):
    def __init__(self, **kwargs):
        # First we have to make a list of ordered steps
        self.steps = ["pdb2gmx", "set_itp", "concat", "editconf",
                      "set_protein_size", "editconf2", "set_protein_size2",
                      "set_popc", "editconf3", "editconf4", "make_topol",
                      "editconf5", "genbox", "set_water", "editconf6",
                      "editconf7", "genbox2", "count_lipids", "make_topol2",
                      "make_topol_lipids", "make_ffoplsaanb", "set_grompp",
                      "grompp", "trjconv", "get_charge", "genion", "grompp2",
                      "trjconv2", "grompp3", "trjconv3"]

        # And then we have to define each step
        self.recipe = \
        {"pdb2gmx": {"gromacs": "pdb2gmx", #0
          "options": {"src": "",
                      "tgt": "proteinopls.pdb",
                      "top": "protein.top"}},
         "set_itp": {"command": "set_itp", #1
          "options": {"src": "protein.top",
                      "tgt": "protein.itp"}},
         "concat": {"command": "concat",
          "options": {"src": "proteinopls.pdb",
                      "tgt": ""}},
         "editconf": {"gromacs": "editconf", #2
          "options": {"src": "proteinopls.pdb",
                      "tgt": "proteinopls.pdb",
                      "dist": ""}},
         "set_protein_size": {"command": "set_protein_size", #3
          "options": {"src": "proteinopls.pdb",
                      "dir": "xy"}},
         "editconf2": {"gromacs": "editconf", #4
          "options": {"src": "proteinopls.pdb",
                      "tgt": "proteinopls.pdb",
                      "dist": ""}},
         "set_protein_size2": {"command": "set_protein_size", #5
          "options": {"src": "proteinopls.pdb",
                      "dir": "z"}},
         "set_popc": {"command": "set_popc", #6
          "options": {"tgt": "popc.pdb"}},
         "editconf3": {"gromacs": "editconf", #7
          "options": {"src": "proteinopls.pdb",
                      "tgt": "proteinopls.pdb",
                      "box": "",
                      "angles": ["90", "90", "120"],
                      "bt": "tric"}},
         "editconf4": {"gromacs": "editconf", #8
          "options": {"src": "popc.pdb",
                      "tgt": "popc.pdb",
                      "box": ""}},
         "make_topol": {"command": "make_topol",
          "options": {}}, #9
         "editconf5": {"gromacs": "editconf", #10
          "options": {"src": "proteinopls.pdb",
                      "tgt": "proteinopls.pdb",
                      "translate": ["0", "0", "0"]}},
         "genbox": {"gromacs": "genbox", #11
          "options": {"cp": "proteinopls.pdb",
                      "cs": "popc.pdb",
                      "tgt": "protpopc.pdb",
                      "top": "topol.top"}},
         "set_water": {"command": "set_water", #12
          "options": {"tgt": "water.pdb"}},
         "editconf6": {"gromacs": "editconf", #13
          "options": {"src": "water.pdb",
                      "tgt": "water.pdb",
                       "box": ""}},
         "editconf7":{"gromacs": "editconf", #14
          "options": {"src": "protpopc.pdb",
                       "tgt": "protpopc.pdb",
                       "box": "",
                       "angles": ["90", "90", "120"],
                       "bt": "tric"}},
         "genbox2": {"gromacs": "genbox", #15
          "options": {"cp": "protpopc.pdb",
                       "cs": "water.pdb",
                       "tgt": "tmp.pdb",
                       "top": "topol.top"}},
         "count_lipids": {"command": "count_lipids", #16
          "options": {"src": "tmp.pdb",
                       "tgt": "popc.pdb"}},
         "make_topol2": {"command": "make_topol",#17
          "options": {}},
         "make_topol_lipids": {"command": "make_topol_lipids"}, #18
         "make_ffoplsaanb": {"command": "make_ffoplsaanb",
          "options": {}},
         "set_grompp": {"command": "set_grompp", #19
          "options": {"steep.mdp": "steep.mdp",
                       "popc.itp": "popc.itp",
                       #"ffoplsaanb_mod.itp": "ffoplsaanb_mod.itp",
                       "ffoplsaabon_mod.itp": "ffoplsaabon_mod.itp",
                       "ffoplsaa_mod.itp": "ffoplsaa_mod.itp"}},
         "grompp": {"gromacs": "grompp", #20
          "options": {"src": "steep.mdp",
                       "src2": "tmp.pdb",
                       "tgt": "topol.tpr",
                       "top": "topol.top"}},
         "trjconv": {"gromacs": "trjconv", #21
          "options": {"src": "tmp.pdb",
                       "src2": "topol.tpr",
                       "tgt": "tmp.pdb",
                       "pbc": "mol"},
          "input": "1\n0\n"},
         "get_charge": {"command": "get_charge",
          "options": {"src": "steep.mdp",
                       "src2": "tmp.pdb",
                       "tgt": "topol.tpr",
                       "top": "topol.top"}}, #22
         "genion": {"gromacs": "genion", #23
          "options": {"src": "topol.tpr",
                       "src2": "topol.top",
                       "tgt": "output.pdb",
                       "np": "",
                       "nn": ""},
          "input": "SOL\n"},
         "grompp2": {"gromacs": "grompp", #24
          "options": {"src": "steep.mdp",
                       "src2": "output.pdb",
                       "tgt": "topol.tpr",
                       "top": "topol.top"}},
         "trjconv2": {"gromacs": "trjconv", #25
          "options": {"src": "output.pdb",
                       "src2": "topol.tpr",
                       "tgt": "output.pdb",
                       "trans": [],
                       "pbc": "mol"},
          "input": "0\n"},
         "grompp3": {"gromacs": "grompp", #26
          "options": {"src": "steep.mdp",
                       "src2": "output.pdb",
                       "tgt": "topol.tpr",
                       "top": "topol.top"}},
         "trjconv3": {"gromacs": "trjconv", #27
          "options": {"src": "output.pdb",
                       "src2": "topol.tpr",
                       "tgt": "hexagon.pdb",
                       "ur": "compact",
                       "pbc": "mol"},
          "input": "1\n0\n"}
           }

        self.breaks = \
            {"pdb2gmx": {"src": "membrane_complex.complex.monomer.pdb_his"},
             "concat": {"tgt": "membrane_complex.complex"},
             "editconf": {"dist": "membrane_complex.box_height"},
             "editconf2": {"dist": "membrane_complex.box_width"},
             "editconf3": {"box": "membrane_complex.trans_box_size"},
             "editconf4": {"box": "membrane_complex.bilayer_box_size"},
             "editconf6": {"box": "membrane_complex.embeded_box_size"},
             "editconf7": {"box": "membrane_complex.protein_box_size"},
             "genion": {"np": "membrane_complex.complex.positive_charge",
                        "nn": "membrane_complex.complex.negative_charge"},
             "make_topol": {"complex": "membrane_complex.complex"},
             "make_topol2": {"complex": "membrane_complex.complex"},
             "make_ffoplsaanb": {"complex": "membrane_complex.complex"},
             "trjconv2": {"trans": "membrane_complex.complex.trans"}
            }

        if kwargs["debug"] or False:
            self.recipe["set_grompp"]["options"]["steep.mdp"] = "steepDEBUG.mdp"

# This recipe modifies the previous one taking a ligand into account
class MonomerLigandRecipe(MonomerRecipe):
    def __init__(self, **kwargs):
        super(MonomerLigandRecipe, self).__init__(**kwargs)
        #self.recipe["set_grompp"]["options"]["ffoplsaanb_mod.itp"] =\
        #    "ffoplsaanb_mod_lig.itp"
        #self.recipe["set_grompp"]["options"]["ffoplsaabon_mod.itp"] =\
        #    "ffoplsaabon_mod_lig.itp"
        #self.recipe["set_grompp"]["options"]["ffoplsaa_mod.itp"] =\
        #    "ffoplsaa_mod_lig.itp"

        self.steps.insert(9, "genrestr")
        self.recipe["genrestr"] =\
            {"gromacs": "genrestr",
             "options": {"src": "",
                         "tgt": "posre_lig.itp",
                         "index": "ligand_ha.ndx",
                         "forces": ["1000", "1000", "1000"]},
             "input": "2\n"}

        self.steps.insert(9, "make_ndx")
        self.recipe["make_ndx"] =\
            {"gromacs": "make_ndx",
             "options": {"src": "",
                         "tgt": "ligand_ha.ndx",
                         "ligand": True},
             "input": "! a H*\nq\n"}

        #self.steps.insert(2, "concat")
        #self.recipe["concat"] =\
        #    {"command": "concat",
        #     "options": {"src": "proteinopls.pdb",
        #                 "tgt": ""}}

        #self.breaks["concat"] =\
        #    {"tgt": "membrane_complex.complex"}
        self.breaks["make_ndx"] =\
            {"src": "membrane_complex.complex.ligand.pdb"}
        self.breaks["genrestr"] =\
            {"src": "membrane_complex.complex.ligand.pdb"}

##########################################################################
#                 Minimization                                           #
##########################################################################

class BasicMinimization(object):
    def __init__(self, **kwargs):
        self.steps = ["set_stage_init", "mdrun"]
        self.recipe = {
         "set_stage_init": {"command": "set_stage_init", #0
          "options": {"src_dir": "",
                      "src_files": ["topol.tpr"],
                      "tgt_dir": "Rmin",
                      "repo_files": ["eq.mdp"]}},
         "mdrun": {"gromacs": "mdrun", #1
          "options": {"dir": "Rmin",
                      "src": "topol.tpr",
                      "tgt": "traj.trj",
                      "energy": "ener.edr",
                      "conf": "confout.gro",
                      "log": "md.log"}},
        }
        self.breaks = {}

        if kwargs["debug"] or False:
            self.recipe["set_stage_init"]["options"]["repo_files"] =\
                ["eqDEBUG.mdp"]

class LigandMinimization(BasicMinimization):
    def __init__(self, **kwargs):
        super(LigandMinimization, self).__init__(**kwargs)

##########################################################################
#                 Equilibration                                          #
##########################################################################

class BasicEquilibration(object):
    def __init__(self, **kwargs):
        self.steps = ["editconf", "make_ndx", "grompp", "set_stage_init",
                      "set_stage_init2", "mdrun"]
        self.recipe = {
         "editconf": {"gromacs": "editconf", #0
          "options": {"src": "Rmin/confout.gro",
                      "tgt": "min.pdb"}},
         "make_ndx": {"command": "make_ndx", #1
          "options": {"src": "min.pdb",
                      "tgt": "index.ndx"}},
         "grompp": {"gromacs": "grompp", #2
          "options": {"src": "Rmin/eq.mdp",
                      "src2": "min.pdb",
                      "top": "topol.top",
                      "tgt": "topol.tpr",
                      "index":"index.ndx"}},
         "set_stage_init": {"command": "set_stage_init", #3
          "options": {"src_dir": "Rmin",
                      "src_files": ["eq.mdp"],
                      "tgt_dir": "eq"}},
         "set_stage_init2": {"command": "set_stage_init", #4
          "options": {"src_dir": "",
                      "src_files": ["topol.tpr", "posre.itp", "posre_hoh.itp",
                                   "posre_ion.itp", "posre_lig.itp"],
                      "tgt_dir": "eq"}},
         "mdrun": {"gromacs": "mdrun", #5
          "options": {"dir": "eq",
                      "src": "topol.tpr",
                      "tgt": "traj.trr",
                      "energy": "ener.edr",
                      "conf": "confout.gro",
                      "traj": "traj.xtc",
                      "log": "md_eq1000.log"}},
        }
        self.breaks = {}

        if kwargs["debug"] or False:
            self.recipe["grompp"]["options"]["src"] = "Rmin/eqDEBUG.mdp"
            self.recipe["set_stage_init"]["options"]["src_files"] =\
                ["eqDEBUG.mdp"]

class LigandEquilibration(BasicEquilibration):
    def __init__(self, **kwargs):
        super(LigandEquilibration, self).__init__(**kwargs)
        self.steps.insert(2, "genrestr")
        self.recipe["genrestr"] = \
            {"gromacs": "genrestr",
             "options": {"src": "Rmin/topol.tpr",
                         "tgt": "protein_ca200.itp",
                         "index": "index.ndx",
                         "forces": ["200", "200", "200"]},
             "input": "3\n"}

##########################################################################
#                    Relaxation                                          #
##########################################################################

class BasicRelax(object):
    def __init__(self, **kwargs):
        self.steps = []
        self.recipe = {}
        for const in range(800, 0, -200):
            self.steps.extend(["relax%d" % const,
                               "grompp%d" % const,
                               "mdrun%d" % const])
            tgt_dir = "eq/{0}".format(const)
            src_dir = "eq"
            self.recipe["relax%d" % const] =\
             {"command": "relax", #0, 3, 6, 9
              "options": {"const": const,
                          "src_dir": src_dir,
                          "tgt_dir": tgt_dir,
                          "mdp": "eq.mdp",
                          "posres": ["posre.itp"]}}
            self.recipe["grompp%d" % const] =\
             {"gromacs": "grompp", #1, 4, 7, 10
              "options": {"src": os.path.join(tgt_dir, "eq.mdp"),
                          "src2": os.path.join(src_dir, "confout.gro"),
                          #top": os.path.join(tgt_dir, "topol.top"),
                          "top": "topol.top",
                          "tgt": os.path.join(tgt_dir, "topol.tpr"),
                          "index": "index.ndx"}}
            #TODO ese confout.gro de abaixo hai que copialo, non vale asi
            self.recipe["mdrun%d" % const] =\
             {"gromacs": "mdrun", #2, 5, 8, 11
              "options": {"dir": tgt_dir,
                          "src": "topol.tpr",
                          "tgt": "traj.trr",
                          "energy": "ener.edr",
                          "conf": "../confout.gro",
                          "traj": "traj.xtc",
                          "log": "md_eq{0}.log".format(const)}}
        self.breaks = {}

        if kwargs["debug"] or False:
            for i in [x for x in self.recipe.keys() if x.startswith("relax")]:
                self.recipe[i]["options"]["mdp"] = "eqDEBUG.mdp"

class LigandRelax(BasicRelax):
    def __init__(self, **kwargs):
        super(LigandRelax, self).__init__(**kwargs)
        self.recipe["relax800"]["options"]["posres"].append("posre_lig.itp")

##########################################################################
#                    Alpha Chains Relaxation                             #
##########################################################################

class CAEquilibrate(object):
    def __init__(self, **kwargs):
        self.steps = ["set_stage_init", "genrestr", "grompp", "mdrun",
                      "trjcat", "eneconv"]
        self.recipe = {
             "set_stage_init": {"command": "set_stage_init", #0
              "options": {"src_dir": "eq",
                          "tgt_dir": "eqCA",
                          "src_files": ["confout.gro"],
                          "repo_files": ["eqCA.mdp"]}},
             "genrestr": {"gromacs": "genrestr", #1
              "options": {"src": "Rmin/topol.tpr",
                          "tgt": "posre.itp",
                          "index": "index.ndx",
                          "forces": ["200"] * 3},
              "input": "3\n"},
             "grompp": {"gromacs": "grompp", #2
              "options": {"src": "eqCA/eqCA.mdp",
                          "src2": "eqCA/confout.gro",
                          "top": "topol.top",
                          "tgt": "eqCA/topol.tpr",
                          "index": "index.ndx"}},
             "mdrun": {"gromacs": "mdrun", #3
              "options": {"dir": "eqCA",
                          "src": "topol.tpr",
                          "tgt": "traj.trr",
                          "energy": "ener.edr",
                          "conf": "confout.gro",
                          "traj": "traj.xtc",
                          "log": "md_eqCA.log"}},
             "trjcat": {"gromacs": "trjcat", #4
              "options": {"dir1": "eq",
                          "dir2": "eqCA",
                          "name": "traj.xtc",
                          "tgt": "traj_EQ.xtc"},
              "input": "c\n" * 6},
             "eneconv": {"gromacs": "eneconv", #5
              "options": {"dir1": "eq",
                          "dir2": "eqCA",
                          "name": "ener.edr",
                          "tgt": "ener_EQ.edr"},
              "input": "c\n" * 6}
             }

        self.breaks = {}

        if kwargs["debug"] or False:
            self.recipe["set_stage_init"]["options"]["src_files"] =\
                ["confout.gro", "eqDEBUG.mdp"]
            self.recipe["grompp"]["options"]["src"] = "eqCA/eqDEBUG.mdp"
