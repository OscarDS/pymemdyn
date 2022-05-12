import os


class BasicInit(object):
    def __init__(self, **kwargs):
        # First we make a list of ordered steps
        self.steps = ["pdb2gmx", "set_itp", "concat", "editconf",
                      "set_protein_size", "editconf2", "set_protein_size2",
                      "set_popc", "editconf3", "editconf4", "make_topol",
                      "editconf5", "solvate", "set_water", "editconf6",
                      "editconf7", "solvate2", "count_lipids", "make_topol2",
                      "make_topol_lipids", "make_ffoplsaanb", "set_grompp",
                      "set_chains", "make_ndx", "grompp", "trjconv",
                      "get_charge", "genion", "grompp2", "trjconv2",
                      "grompp3", "trjconv3", "clean_pdb"]

        # And then we define each step
        self.recipe = \
            {"pdb2gmx": {"gromacs": "pdb2gmx",  # 1
                         "options": {"src": "",
                                     "tgt": "proteinopls.pdb",
                                     "top": "protein.top"}},

             "set_itp": {"command": "set_itp",  # 2
                         "options": {"src": "protein.top",
                                     "tgt": "protein.itp"}},

             "concat": {"command": "concat",  # 3
                        "options": {"src": "proteinopls.pdb",
                                    "tgt": ""}},

             "editconf": {"gromacs": "editconf",  # 4
                          "options": {"src": "proteinopls.pdb",
                                      "tgt": "proteinopls.pdb",
                                      "dist": ""}},

             "set_protein_size": {"command": "set_protein_size",  # 5
                                  "options": {"src": "proteinopls.pdb",
                                              "dir": "xy"}},

             "editconf2": {"gromacs": "editconf",  # 6
                           "options": {"src": "proteinopls.pdb",
                                       "tgt": "proteinopls.pdb",
                                       "dist": ""}},

             "set_protein_size2": {"command": "set_protein_size",  # 7
                                   "options": {"src": "proteinopls.pdb",
                                               "dir": "z"}},

             "set_popc": {"command": "set_popc",  # 8
                          "options": {"tgt": "popc.pdb"}},

             "editconf3": {"gromacs": "editconf",  # 9
                           "options": {"src": "proteinopls.pdb",
                                       "tgt": "proteinopls.pdb",
                                       "box": "",
                                       "angles": ["90", "90", "120"],
                                       "bt": "tric"}},

             "editconf4": {"gromacs": "editconf",  # 10
                           "options": {"src": "popc.pdb",
                                       "tgt": "popc.gro",
                                       "box": ""}},

             "make_topol": {"command": "make_topol",  # 11
                            "options": {}},

             "editconf5": {"gromacs": "editconf",  # 12
                           "options": {"src": "proteinopls.pdb",
                                       "tgt": "proteinopls.pdb",
                                       "translate": ["0", "0", "0"]}},

             "solvate": {"gromacs": "solvate",  # 13
                        "options": {"cp": "proteinopls.pdb",
                                    "cs": "popc.gro",
                                    "tgt": "protpopc.pdb",
                                    "top": "topol.top"}},

             "set_water": {"command": "set_water",  # 14
                           "options": {"tgt": "water.pdb"}},

             "editconf6": {"gromacs": "editconf",  # 15
                           "options": {"src": "water.pdb",
                                       "tgt": "water.gro",
                                       "box": ""}},

             "editconf7": {"gromacs": "editconf",  # 16
                           "options": {"src": "protpopc.pdb",
                                       "tgt": "protpopc.pdb",
                                       "box": "",
                                       "angles": ["90", "90", "120"],
                                       "bt": "tric"}},

             "solvate2": {"gromacs": "solvate",  # 17
                         "options": {"cp": "protpopc.pdb",
                                     "cs": "water.gro",
                                     "tgt": "tmp.pdb",
                                     "top": "topol.top"}},

             "count_lipids": {"command": "count_lipids",  # 18
                              "options": {"src": "tmp.pdb",
                                          "tgt": "popc.pdb"}},

             "make_topol2": {"command": "make_topol",  # 19
                             "options": {}},

             "make_topol_lipids": {"command": "make_topol_lipids"},  # 20

             "make_ffoplsaanb": {"command": "make_ffoplsaanb",  # 21
                                 "options": {}},

             "set_grompp": {"command": "set_grompp",  # 22
                            "options": {"steep.mdp": "steep.mdp",
                                        "popc.itp": "popc.itp",
                                        "spc.itp": "spc.itp",
                                        "ions.itp": "ions.itp",
                                        "ffoplsaabon_mod.itp": "ffoplsaabon_mod.itp",
                                        "ffoplsaa_mod.itp": "ffoplsaa_mod.itp"}},

             "set_chains": {"command": "set_chains",  # 23
                            "options": {"src": "proteinopls.pdb"}},

             "make_ndx": {"command": "make_ndx",  # 24
                          "options": {"src": "tmp.pdb",
                                      "tgt": "index.ndx"}},

             "grompp": {"gromacs": "grompp",  # 25
                        "options": {"src": "steep.mdp",
                                    "src2": "tmp.pdb",
                                    "tgt": "topol.tpr",
                                    "top": "topol.top",
                                    "index": "index.ndx"}},

             "trjconv": {"gromacs": "trjconv",  # 26
                         "options": {"src": "tmp.pdb",
                                     "src2": "topol.tpr",
                                     "tgt": "tmp.pdb",
                                     "pbc": "mol",
                                     "index": "index.ndx"},
                         "input": "1\n0\n"},

             "get_charge": {"command": "get_charge",  # 27
                            "options": {"src": "steep.mdp",
                                        "src2": "tmp.pdb",
                                        "tgt": "topol.tpr",
                                        "top": "topol.top",
                                        "index": "index.ndx"}},

             "genion": {"gromacs": "genion",  # 28
                        "options": {"src": "topol.tpr",
                                    "tgt": "output.pdb",
                                    "src2": "topol.top",
                                    "index": "index.ndx",
                                    "np": "",
                                    "nn": ""},
                        "input": " SOL \n"},

             "grompp2": {"gromacs": "grompp",  # 29
                         "options": {"src": "steep.mdp",
                                     "src2": "output.pdb",
                                     "tgt": "topol.tpr",
                                     "top": "topol.top"}},

             "trjconv2": {"gromacs": "trjconv",  # 30
                          "options": {"src": "output.pdb",
                                      "src2": "topol.tpr",
                                      "tgt": "output.pdb",
                                      "trans": [],
                                      "pbc": "mol"},
                          "input": "0\n"},

             "grompp3": {"gromacs": "grompp",  # 31
                         "options": {"src": "steep.mdp",
                                     "src2": "output.pdb",
                                     "tgt": "topol.tpr",
                                     "top": "topol.top"}},

             "trjconv3": {"gromacs": "trjconv",  # 32
                          "options": {"src": "output.pdb",
                                      "src2": "topol.tpr",
                                      "tgt": "hexagon.pdb",
                                      "ur": "compact",
                                      "pbc": "mol"},
                          "input": "1\n0\n"},
             
             "clean_pdb": {"command": "clean_pdb", # 33
                           "options": {"src": "hexagon.pdb",
                                       "tgt": "hexagon.pdb"}}
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
             "genion": {"nn": "membrane_complex.complex.positive_charge",
                        "np": "membrane_complex.complex.negative_charge"},
             "make_topol": {"complex": "membrane_complex.complex"},
             "make_topol2": {"complex": "membrane_complex.complex"},
             "make_ffoplsaanb": {"complex": "membrane_complex.complex"},
             "trjconv2": {"trans": "membrane_complex.complex.trans"}
             }

        if kwargs["debug"] or False:
            self.recipe["set_grompp"]["options"]["steep.mdp"] = "steepDEBUG.mdp"


# This recipe modifies the previous one taking a ligand into account
class LigandInit(BasicInit):
    def __init__(self, **kwargs):
        super(LigandInit, self).__init__(**kwargs)

        self.steps.insert(9, "genrestr_lig")
        self.recipe["genrestr_lig"] = \
            {"gromacs": "genrestr", 
             "options": {"src": "",
                         "tgt": "posre_lig.itp",
                         "index": "ligand_ha.ndx",
                         "forces": ["1000", "1000", "1000"]},
             "input": "2\n"}

        self.steps.insert(9, "make_ndx_lig")
        self.recipe["make_ndx_lig"] = \
            {"gromacs": "make_ndx",
             "options": {"src": "",
                         "tgt": "ligand_ha.ndx",
                         "ligand": True},
             "input": "! a H*\nq\n"}

        self.breaks["make_ndx_lig"] = \
            {"src": "membrane_complex.complex.ligand.pdb"}
        self.breaks["genrestr_lig"] = \
            {"src": "membrane_complex.complex.ligand.pdb"}


# This recipe modifies the previous one taking an alosteric ligand into account
class LigandAlostericInit(LigandInit):
    def __init__(self, **kwargs):
        super(LigandAlostericInit, self).__init__(**kwargs)

        self.steps.insert(9, "genrestr_alo")
        self.recipe["genrestr_alo"] = \
            {"gromacs": "genrestr",
             "options": {"src": "",
                         "tgt": "posre_alo.itp",
                         "index": "alosteric_ha.ndx",
                         "forces": ["1000", "1000", "1000"]},
             "input": "2\n"}

        self.steps.insert(9, "make_ndx_alo")
        self.recipe["make_ndx_alo"] = \
            {"gromacs": "make_ndx",
             "options": {"src": "",
                         "tgt": "alosteric_ha.ndx",
                         "alosteric": True},
             "input": "! a H*\nq\n"}

        self.breaks["make_ndx_alo"] = \
            {"src": "membrane_complex.complex.alosteric.pdb"}
        self.breaks["genrestr_alo"] = \
            {"src": "membrane_complex.complex.alosteric.pdb"}


##########################################################################
#                 Minimization                                           #
##########################################################################

class BasicMinimization(object):
    def __init__(self, **kwargs):
        self.steps = ["set_stage_init", "mdrun"]
        self.recipe = {
            "set_stage_init": {"command": "set_stage_init",  # 1
                               "options": {"src_dir": "",
                                           "src_files": ["topol.tpr"],
                                           "tgt_dir": "Rmin",
                                           "repo_files": ["eq.mdp"]}},

            "mdrun": {"gromacs": "mdrun",  # 2
                      "options": {"dir": "Rmin",
                                  "src": "topol.tpr",
                                  "tgt": "traj.trj",
                                  "energy": "ener.edr",
                                  "conf": "confout.gro",
                                  "log": "md.log"}},
        }
        self.breaks = {}

        if kwargs["debug"] or False:
            self.recipe["set_stage_init"]["options"]["repo_files"] = \
                ["eqDEBUG.mdp"]


class LigandMinimization(BasicMinimization):
    def __init__(self, **kwargs):
        super(LigandMinimization, self).__init__(**kwargs)


class LigandAlostericMinimization(BasicMinimization):
    def __init__(self, **kwargs):
        super(LigandAlostericMinimization, self).__init__(**kwargs)


##########################################################################
#                 Equilibration                                          #
##########################################################################

class BasicEquilibration(object):
    def __init__(self, **kwargs):
        self.steps = ["editconf", "make_ndx", "grompp", "set_stage_init",
                      "set_stage_init2", "mdrun"]
        self.recipe = {
            "editconf": {"gromacs": "editconf",  # 1
                         "options": {"src": "Rmin/confout.gro",
                                     "tgt": "min.pdb"}},

            "make_ndx": {"command": "make_ndx",  # 2
                         "options": {"src": "min.pdb",
                                     "tgt": "index.ndx"}},

            "grompp": {"gromacs": "grompp",  # 3
                       "options": {"src": "Rmin/eq.mdp",
                                   "src2": "min.pdb",
                                   "top": "topol.top",
                                   "tgt": "topol.tpr",
                                   "index": "index.ndx"}},

            "set_stage_init": {"command": "set_stage_init",  # 4
                               "options": {"src_dir": "Rmin",
                                           "src_files": ["eq.mdp"],
                                           "tgt_dir": "eq"}},

            "set_stage_init2": {"command": "set_stage_init",  # 5
                                "options": {"src_dir": "",
                                            "src_files": ["topol.tpr",
                                                          "posre.itp",
                                                          "posre_Protein_chain_A.itp",
                                                          "posre_Protein_chain_B.itp",
                                                          "posre_hoh.itp",
                                                          "posre_ion.itp",
                                                          "posre_lig.itp",
                                                          "posre_alo.itp",
                                                          "posre_cho.itp"],
                                            "tgt_dir": "eq"}},

            "mdrun": {"gromacs": "mdrun",  # 6
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
            self.recipe["set_stage_init"]["options"]["src_files"] = \
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


class LigandAlostericEquilibration(LigandEquilibration):
    def __init__(self, **kwargs):
        super(LigandAlostericEquilibration, self).__init__(**kwargs)
        self.steps.insert(2, "genrestr")
        self.recipe["genrestr"] = \
            {"gromacs": "genrestr",
             "options": {"src": "Rmin/topol.tpr",
                         "tgt": "protein_ca200.itp",
                         "index": "index.ndx",
                         "forces": ["200", "200", "200"]},
             "input": "3\n"}


##########################################################################
#                 Relaxation                                             #
##########################################################################

class BasicRelax(object):
    def __init__(self, **kwargs):
        self.steps = []
        self.recipe = {}
        for const in range(800, 0, -200):
            self.steps.extend(["relax{0}".format(const),
                               "grompp{0}".format(const),
                               "mdrun{0}".format(const)])
            tgt_dir = "eq/{0}".format(const)
            src_dir = "eq"
            self.recipe["relax%d" % const] = \
                {"command": "relax",  # 1, 4, 7, 10
                 "options": {"const": const,
                             "src_dir": src_dir,
                             "tgt_dir": tgt_dir,
                             "posres": [],
                             "mdp": "eq.mdp"}}

            self.recipe["grompp%d" % const] = \
                {"gromacs": "grompp",  # 2, 5, 8, 11
                 "options": {"src": os.path.join(tgt_dir, "eq.mdp"),
                             "src2": os.path.join(src_dir, "confout.gro"),
                             "top": "topol.top",
                             "tgt": os.path.join(tgt_dir, "topol.tpr"),
                             "index": "index.ndx"}}
            # TODO the confout.gro coming next needs to be copied, it feels a
            # bit like cheating

            self.recipe["mdrun%d" % const] = \
                {"gromacs": "mdrun",  # 3, 6, 9, 12
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


class LigandAlostericRelax(LigandRelax):
    def __init__(self, **kwargs):
        super(LigandAlostericRelax, self).__init__(**kwargs)
        self.recipe["relax800"]["options"]["posres"].append("posre_alo.itp")


##########################################################################
#                 C-Alpha Atoms Relaxation                               #
##########################################################################

class BasicCARelax(object):
    def __init__(self, **kwargs):
        self.steps = ["set_stage_init", "genrestr", "grompp", "mdrun"]
        self.recipe = {
            "set_stage_init": {"command": "set_stage_init",  # 1
                               "options": {"src_dir": "eq",
                                           "tgt_dir": "eqProd",
                                           "src_files": ["confout.gro"],
                                           "repo_files": ["eqCA.mdp"]}},

            "genrestr": {"gromacs": "genrestr",  # 2
                         "options": {"src": "Rmin/topol.tpr",
                                     "tgt": "posre.itp",
                                     "index": "index.ndx",
                                     "forces": ["200"] * 3},
                         "input": "3\n"},

            "grompp": {"gromacs": "grompp",  # 3
                       "options": {"src": "eqProd/eqCA.mdp",
                                   "src2": "eqProd/confout.gro",
                                   "top": "topol.top",
                                   "tgt": "eqProd/topol.tpr",
                                   "index": "index.ndx"}},

            "mdrun": {"gromacs": "mdrun",  # 4
                      "options": {"dir": "eqProd",
                                  "src": "topol.tpr",
                                  "tgt": "traj.trr",
                                  "energy": "ener.edr",
                                  "conf": "confout.gro",
                                  "traj": "traj.xtc",
                                  "log": "md_eqProd.log"}},
        }

        self.breaks = {}

        if kwargs["debug"] or False:
            self.recipe["set_stage_init"]["options"]["src_files"] = \
                ["confout.gro", "eqDEBUG.mdp"]
            self.recipe["grompp"]["options"]["src"] = "eqProd/eqDEBUG.mdp"


##########################################################################
#                 Ballesteros-Weinstein Restrained Relaxation            #
##########################################################################

class BasicBWRelax(object):
    def __init__(self, **kwargs):
        self.steps = ["getbw","set_stage_init", "grompp", "mdrun"]
        self.recipe = {
            "getbw": {"command": "getbw",  # 1
                      "options": {"src": "proteinopls.pdb"}},

            "set_stage_init": {"command": "set_stage_init",  # 2
                               "options": {"src_dir": "eq",
                                           "tgt_dir": "eqProd",
                                           "src_files": ["confout.gro",
                                                         "../disre.itp"],
                                           "repo_files": ["dres.mdp"]}},

            "grompp": {"gromacs": "grompp",  # 3
                       "options": {"src": "eqProd/dres.mdp",
                                   "src2": "eqProd/confout.gro",
                                   "top": "topol.top",
                                   "tgt": "eqProd/topol.tpr",
                                   "index": "index.ndx"}},

            "mdrun": {"gromacs": "mdrun",  # 4
                      "options": {"dir": "eqProd",
                                  "src": "topol.tpr",
                                  "tgt": "traj.trr",
                                  "energy": "ener.edr",
                                  "conf": "confout.gro",
                                  "traj": "traj.xtc",
                                  "log": "md_eqBW.log"}},
        }

        self.breaks = {}

        if kwargs["debug"] or False:
            self.recipe["set_stage_init"]["options"]["src_files"] = \
                ["confout.gro", "eqDEBUG.mdp"]
            self.recipe["grompp"]["options"]["src"] = "eqProd/eqDEBUG.mdp"


##########################################################################
#                 Collect All Results & Output                           #
##########################################################################
class BasicCollectResults(object):
    def __init__(self, **kwargs):
        """
        This recipe navigates through the output of *GROMACS*, generating and
        collecting various files and logs, and finally putting it all together
        in a tar file.

        Provides a list of *steps* as the order to execute the recipe.
        dict *recipe* as the variables to pass to functions.
        dict *breaks* as points where object calling can put their vars.
        """
        self.breaks = {}
        self.steps = ["trjcat", "trjconv", "rms1", "rms2",
                      "rms3", "rmsf", "tot_ener", "temp", "pressure",
                      "volume", "set_end", "clean_topol", "set_end_2",
                      "set_end_3", "set_end_4", "set_end_5", "set_end_6",
                      "tar_it"]

        self.recipe = {"trjcat":
                           {"gromacs": "trjcat",  # 1
                            "options": {"dir1": "eq",
                                        "dir2": "eqProd",
                                        "name": "traj.xtc",
                                        "tgt": "traj_EQ.xtc"},
                            "input": "c\n" * 6},

                       "trjconv":
                           {"gromacs": "trjconv",  # 2
                                   "options": {"src": "traj_EQ.xtc",
                                               "src2": "topol.tpr",
                                               "tgt": "traj_pymol.xtc",
                                               "ur": "compact",
                                               "skip": "2",
                                               "pbc": "mol"},
                                   "input": "1\n0\n"},

                       "rms1":
                           {"gromacs": "rms",  # 4
                                 "options": {"src": "eq/topol.tpr",
                                             "src2": "traj_EQ.xtc",
                                             "tgt": "rmsd-all-atom-vs-start.xvg"},
                                 "input": "1\n1\n"},

                       "rms2":
                           {"gromacs": "rms",  # 5
                                 "options": {"src": "eq/topol.tpr",
                                             "src2": "traj_EQ.xtc",
                                             "tgt": "rmsd-backbone-vs-start.xvg"},
                                 "input": "4\n4\n"},

                       "rms3":
                           {"gromacs": "rms",  # 6
                                 "options": {"src": "eq/topol.tpr",
                                             "src2": "traj_EQ.xtc",
                                             "tgt": "rmsd-calpha-vs-start.xvg"},
                                 "input": "3\n3\n"},

                       "rmsf":
                           {"gromacs": "rmsf",  # 7
                                 "options": {"src": "eq/topol.tpr",
                                             "src2": "traj_EQ.xtc",
                                             "tgt": "rmsf-per-residue.xvg"},
                                 "input": "1\n"},

                       "set_end":
                           {"command": "set_stage_init",  # 12
                                   "options": {"src_dir": "eqProd",
                                               "src_files": ["confout.gro",
                                                             "topol.tpr"],
                                               "repo_files": ["popc.itp",
                                                              "README.md",
                                                              "prod.mdp",
                                                              "load_gpcr.pml"],
                                               "tgt_dir": "finalOutput"}},

                       "clean_topol":
                           {"command": "clean_topol",  # 13
                                       "options": {"src": "topol.top",
                                                   "tgt": "finalOutput/topol.top"}},

                       "set_end_2":
                           {"command": "set_stage_init",  # 14
                                     "options": {"src_dir": "",
                                                 "src_files": [
                                                     "ffoplsaa_mod.itp",
                                                     "ffoplsaabon_mod.itp",
                                                     "ffoplsaanb_mod.itp",
                                                     "hexagon.pdb",
                                                     "protein.itp",
                                                     "lig.itp",
                                                     "ions_local.itp",
                                                     "cho.itp",
                                                     "hoh.itp",
                                                     "index.ndx", 
                                                     "traj_EQ.xtc",
                                                     "ener_EQ.edr",
                                                     "traj_pymol.xtc"],
                                                 "tgt_dir": "finalOutput"}},

                       "set_end_3":
                           {"command": "set_stage_init",  # 15
                                     "options": {"src_dir": "",
                                                 "src_files": ["tot_ener.xvg",
                                                               "tot_ener.log",
                                                               "temp.xvg",
                                                               "temp.log",
                                                               "pressure.xvg",
                                                               "pressure.log",
                                                               "volume.xvg",
                                                               "volume.log",                                                     
                                                               "rmsd-all-atom-vs-start.xvg",
                                                               "rmsd-calpha-vs-start.xvg",
                                                               "rmsd-backbone-vs-start.xvg",
                                                               "rmsf-per-residue.xvg"],
                                                 "tgt_dir": "finalOutput/reports"}},

                       "set_end_4":
                           {"command": "set_stage_init",  # 16
                                     "options": {"src_dir": "eq",
                                                 "src_files": ["md_eq1000.log"],
                                                 "tgt_dir": "finalOutput/logs"}},

                       "set_end_5":
                           {"command": "set_stage_init",  # 17
                                     "options": {"src_dir": "eq",
                                                 "src_files": [
                                                     "{0}/md_eq{0}.log".format(
                                                         const)
                                                     for const in
                                                     range(800, 0, -200)],
                                                 "tgt_dir": "finalOutput/logs"}},

                       "set_end_6":
                           {"command": "set_stage_init",  # 18
                                     "options": {"src_dir": "eqProd",
                                                 "src_files": ["md_eqProd.log"],
                                                 "tgt_dir": "finalOutput/logs"}},

                       "tar_it":
                           {"command": "tar_out", # 19
                                  "options": {"src_dir": "finalOutput",
                                              "tgt": "MD_output.tgz"}}, 
                       }

        options = {"tot_ener": "Total-Energy\n",
                   "temp": "Temperature\n",
                   "pressure": "Pressure\n",
                   "volume": "Volume\n"} # 8 to 11

        for option, gro_key in options.items():
            self.recipe[option] = \
                {"gromacs": "energy",
                 "options": {"src": "ener_EQ.edr",
                             "tgt": "{0}.xvg".format(option),
                             "log": "{0}.log".format(option)},
                 "input": gro_key}


class BasicCACollectResults(BasicCollectResults):
    def __init__(self, **kwargs):
        super(BasicCACollectResults, self).__init__(**kwargs)
        self.steps.insert(2, "eneconv")
        self.recipe["eneconv"] = \
                {"gromacs": "eneconv", 
                 "options": {"dir1": "eq",
                         "dir2": "eqProd",
                         "name": "ener.edr",
                         "tgt": "ener_EQ.edr"},
                 "input": "y\nc\nc\nc\nc\nc\nc\n"}


class BasicBWCollectResults(BasicCollectResults):
    def __init__(self, **kwargs):
        super(BasicBWCollectResults, self).__init__(**kwargs)
        self.steps.insert(14, "set_end_BW")
        self.recipe["set_end_BW"] = \
                {"command": "set_stage_init",
                          "options": {"src_dir": "",
                                      "src_files": ["tot_ener2.xvg",
                                                    "tot_ener2.log",
                                                    "temp2.xvg",
                                                    "temp2.log",
                                                    "pressure2.xvg",
                                                    "pressure2.log",
                                                    "volume2.xvg",
                                                    "volume2.log"],                                                     
                                      "tgt_dir": "finalOutput/reports"}}

        self.steps.insert(10, "volume2")
        self.steps.insert(10, "pressure2")
        self.steps.insert(10, "temp2")
        self.steps.insert(10, "tot_ener2")
            
        options2 = {"tot_ener2": "Total-Energy\n",
                        "temp2": "Temperature\n",
                        "pressure2": "Pressure\n",
                        "volume2": "Volume\n"}

        for option, gro_key in options2.items():
            self.recipe[option] = \
                {"gromacs": "energy",
                 "options": {"src": "eqProd/ener.edr",
                             "tgt": "{0}.xvg".format(option),
                             "log": "{0}.log".format(option)},
                 "input": gro_key}

        self.steps.insert(2, "eneconv")
        self.recipe["eneconv"] = \
                {"gromacs": "eneconv", 
                 "options": {"dir1": "eq",
                             "dir2": "",
                             "name": "ener.edr",
                             "tgt": "ener_EQ.edr"},
                 "input": "y\nc\nc\nc\nc\nc\nc\n"}