import groerrors
import recipes
import settings
import utils

import logging
import os
import shutil
import subprocess
import sys

class Gromacs(object):
    def __init__(self, *args, **kwargs):
        self._membrane_complex = None
        self.wrapper = Wrapper()
        logging.basicConfig(filename="GROMACS.log",
                            level=logging.DEBUG,
                            format='%(message)s')
        if "membrane_complex" in kwargs.keys():
            self.set_membrane_complex(kwargs["membrane_complex"])
            self.tpr = \
                self.membrane_complex.complex.monomer.pdb.replace(".pdb",
                    ".tpr")

    def set_membrane_complex(self, value):
        '''Sets the monomer object'''
        self._membrane_complex = value
    def get_membrane_complex(self):
        return self._membrane_complex
    membrane_complex = property(get_membrane_complex, set_membrane_complex)

    def count_lipids(self, **kwargs):
        '''Count the lipids in source and write a target with N4 tags'''
        src = open(kwargs["src"], "r")
        tgt = open(kwargs["tgt"], "w")

        half = self.membrane_complex.complex.prot_z / 2
        self.membrane_complex.membrane.lipids_up = 0
        self.membrane_complex.membrane.lipids_down = 0
        self.membrane_complex.membrane.n_wats = 0

        if getattr(self.membrane_complex.complex, "waters"):
            #Careful, some waters belong to crystal, not solvent
            self.membrane_complex.membrane.n_wats -=\
            self.membrane_complex.complex.waters.number

        for line in src:
            if len(line.split()) > 2:
                if line.split()[2] == "N4": #Lipid marker
                    tgt.write(line)
                    if(float(line.split()[7]) >= half):
                        self.membrane_complex.membrane.lipids_up += 1
                    elif(float(line.split()[7]) < half):
                        self.membrane_complex.membrane.lipids_down += 1
                elif line.split()[2] == "OW": #Water marker
                    self.membrane_complex.membrane.n_wats += 1

        src.close()
        tgt.close()
    
        return True

    def get_charge(self, **kwargs):
        '''Get the total charge of a system using gromacs grompp command'''
        #wrapper = Wrapper()

        out, err = self.wrapper.run_command({"gromacs": "grompp",
                                             "options": kwargs})

        # Now we are looking for this line:
        # System has non-zero total charge: 6.000002e+00
        charge = 0
        for line in err.split("\n"):
            if "total charge" in line:
                charge = abs(int(float(line.split()[-1])))
                break

        self.membrane_complex.complex.negative_charge = 0
        self.membrane_complex.complex.positive_charge = 0
        if charge > 0:
            self.membrane_complex.complex.positive_charge = charge
            self.membrane_complex.complex.negative_charge = 0
        elif charge < 0:
            self.membrane_complex.complex.negative_charge = charge
            self.membrane_complex.complex.positive_charge = 0
        else:
            self.membrane_complex.complex.negative_charge = 0
            self.membrane_complex.complex.positive_charge = 0

        return True

    def get_ndx_groups(self, **kwargs):
        '''Run make_ndx and sets the group Other'''
        #wrapper = Wrapper()

        out, err = self.wrapper.run_command({"gromacs": "make_ndx",
                                             "options": kwargs,
                                             "input": "q\n"})

        for line in out.split("\n"):
            if "Other" in line and "atoms" in line:
                self.n_groups = int(line.split()[0])
                return True
        return False

    def make_ndx(self, **kwargs):
        '''Wraps the make_ndx command tweaking the input to reflect the
        complex characteristics.'''

        if not (self.get_ndx_groups(**kwargs)): return False

        n_group = self.n_groups + 1

        #This makes the (always present) group wation
        input =  "r SOL || r HOH || atom Cl* || atom Na*\n"
        input += "name {0} wation\n".format(n_group)

        #This makes the group protlig
        n_group += 1
        input += "1 || r LIG || r ALO\n"
        input += "name {0} protlig\n".format(n_group)

        #And this makes the membrane group as membr
        n_group += 1
        input += "r POP || r CHO || r LIP\n"
        input += "name {0} membr\n".format(n_group)

        if hasattr(self.membrane_complex.membrane, "ions"):
            #This makes the group ions TODO
            #n_group += 1
            #input += "1 || r LIG\nname {0} protlig\n".format(n_group)
            pass
        input += "q\n"

        #Now the wrap itself
        out, err = self.wrapper.run_command({"gromacs": "make_ndx",
                                             "options": kwargs,
                                             "input": input})

        logging.debug("make_ndx command")
        logging.debug(err)
        logging.debug(out)
        #We need to set the eq.mdp with these new groups
        #NO LONGER NEEDED AS eq.mdp IS NOW GENERIC
        #utils.tune_mdp(groups)

        return True

    def make_topol_lipids(self, **kwargs):
        '''Add the lipid positions to our topol.top'''
        topol = open("topol.top", "a")
        topol.write("; Number of POPC molecules with higher z-coord value:\n")
        topol.write("POPC " + str(self.membrane_complex.membrane.lipids_up))
        topol.write("\n; Number of POPC molecules with lower z-coord value:\n")
        topol.write("POPC " + str(self.membrane_complex.membrane.lipids_down))
        topol.write("\n; Total number of water molecules:\n")
        topol.write("SOL " + str(self.membrane_complex.membrane.n_wats) + "\n")
        topol.close()

        return True

    def passing(self):
        '''Do nothing, to respect some orders'''
        return True

    def relax(self, **kwargs):
        '''Relax a protein'''

        if not os.path.isdir(kwargs["tgt_dir"]): os.makedirs(kwargs["tgt_dir"])

        if hasattr(self.membrane_complex.complex, "waters") and\
            self.membrane_complex.complex.waters:
            kwargs["posres"].append("posre_hoh.itp")
        if hasattr(self.membrane_complex.complex, "ions") and\
            self.membrane_complex.complex.ions:
            kwargs["posres"].append("posre_ion.itp")
        if hasattr(self.membrane_complex.complex, "cho") and\
            self.membrane_complex.complex.cho:
            kwargs["posres"].append("posre_cho.itp")

        for posre in kwargs["posres"]:
            new_posre = open(os.path.join(kwargs["tgt_dir"], posre), "w")

            for line in open(os.path.join(kwargs["src_dir"], posre), "r"):
                if line.split()[-3:] == ["1000", "1000", "1000"]:
                    new_posre.write(" ".join(line.split()[:2] +\
                                             [str(kwargs["const"])] * 3))
                    new_posre.write("\n")
                else:
                    new_posre.write(line)
            new_posre.close()

        new_mdp = open(os.path.join(kwargs["tgt_dir"], "eq.mdp"), "w")
        src_mdp = open(os.path.join(kwargs["src_dir"], kwargs["mdp"]), "r")

        for line in src_mdp:
            if(line.startswith("gen_vel")):
                new_mdp.write("gen_vel = no\n")
            else:
                new_mdp.write(line)
        src_mdp.close()
        new_mdp.close()

        utils.make_topol(target_dir = kwargs["tgt_dir"],
            working_dir = os.getcwd(),
            complex = self.membrane_complex.complex)

    def run_recipe(self, debug = False):
        '''Run the recipe for the complex'''
        if not hasattr(self, "recipe"):
            self.select_recipe(debug = debug)

        #wrapper = Wrapper()
        self.repo_dir = self.wrapper.repo_dir

        for n, command_name in enumerate(self.recipe.steps):
            command = self.recipe.recipe[command_name]

            if command_name in self.recipe.breaks.keys():
                command["options"] = self.set_options(command["options"],
                                     self.recipe.breaks[command_name])

            # NOW RUN IT !
            print self.recipe.__class__.__name__, 
            print "Step {0}: ".format(command_name),
            if command.has_key("gromacs"):
                # Either run a Gromacs pure command...
                print command["gromacs"]
                if hasattr(self, "queue"): command["queue"] = self.queue
                out, err = self.wrapper.run_command(command)
                logging.debug(" ".join(self.wrapper.generate_command(command)))
                logging.debug(err)
                logging.debug(out)
                #This test the Gromacs output checking for known errors
                groerror.GromacsMessages(gro_err=err, command=command[0])
                
            else:
                # ...or run a local function
                logging.debug(command)
                try:
                    f = getattr(self, command["command"])
                except AttributeError:
                    #Fallback to the utils module
                    f = getattr(utils, command["command"])

                if command.has_key("options"): f(**command["options"])
                else: f()
                logging.debug("This function does: " + str(f.__doc__))
                print command["command"]

        return True

    def select_recipe(self, stage = "", debug = False):
        '''Select the appropiate recipe for the complex'''
        recipe = ""
        stage = stage or "Init"

        if self.membrane_complex:
            if not hasattr(self.membrane_complex.complex, "ligand"):
                recipe += "Basic"
            if hasattr(self.membrane_complex.complex, "ligand"):
                if self.membrane_complex.complex.ligand:
                    recipe += "Ligand"
            if hasattr(self.membrane_complex.complex, "alosteric"):
                if self.membrane_complex.complex.alosteric:
                    recipe += "Alosteric"

        recipe += stage #This kwarg carries the proper recipe:
                        #Init, Minimization, Equilibration...

        self.recipe = getattr(recipes, recipe)(debug = debug)

        return True

    def set_box_sizes(self):
        '''Set some meassurements of the different boxes'''

        self.membrane_complex.complex.set_nanom()
        self.membrane_complex.trans_box_size = \
            [str(self.membrane_complex.complex.gmx_prot_xy),
             str(self.membrane_complex.complex.gmx_prot_xy),
             str(self.membrane_complex.membrane.gmx_bilayer_z)]
        self.membrane_complex.bilayer_box_size = \
            [str(self.membrane_complex.membrane.gmx_bilayer_x),
             str(self.membrane_complex.membrane.gmx_bilayer_y),
             str(self.membrane_complex.membrane.gmx_bilayer_z)]
        self.membrane_complex.embeded_box_size = \
            [str(self.membrane_complex.membrane.gmx_bilayer_x),
             str(self.membrane_complex.membrane.gmx_bilayer_y),
             str(self.membrane_complex.complex.gmx_prot_z)]
        self.membrane_complex.protein_box_size = \
            [str(self.membrane_complex.complex.gmx_prot_xy),
             str(self.membrane_complex.complex.gmx_prot_xy),
             str(self.membrane_complex.complex.gmx_prot_z)]

        return True

    def set_grompp(self, **kwargs):
        '''Copy the template files to the working dir'''
        for repo_src in kwargs.keys():
            shutil.copy(os.path.join(
                self.repo_dir, kwargs[repo_src]),
                repo_src)

        return True
 
    def set_itp(self, **kwargs):
        '''Cut a top file to be usable later as itp'''
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
                self.membrane_complex.complex.monomer.name = line.split()[0]
                get_name = False

            if line.startswith("[ moleculetype ]"):
                get_name = True

        tgt.close()
        src.close()

        return True

    def set_options(self, options, breaks):
        '''Set the break options from a recipe'''
        for option, value in breaks.iteritems():
            # This is a hack to get the attribute recursively,
            # feeding getattr with dot-splitted string thanks to reduce
            # Here we charge some commands with options calculated
            new_option = reduce(getattr,
                         value.split("."),
                         self)
            options[option] = new_option

        return options

    def set_popc(self, tgt = ""):
        '''Create a pdb file only with the lipid bilayer (POP), no waters.
        Set some measures on the fly (height of the bilayer)'''
        tgt = open(tgt, "w")
        pops = []

        for line in open(self.membrane_complex.membrane.pdb, "r"):
            if len(line.split()) > 7:
                if line.split()[3] == "POP":
                    tgt.write(line)
                    pops.append(float(line.split()[7]))

        tgt.close()
        self.membrane_complex.membrane.bilayer_z = max(pops) - min(pops)

        return True

    def set_protein_size(self, **kwargs):
        '''Get the protein max base wide from a pdb file'''
        for line in open(kwargs["src"], "r"):
            if line.startswith("CRYST1"):
                if kwargs["dir"] == "xy":
                    self.membrane_complex.complex.prot_xy = \
                        max(float(line.split()[1]), #xyprotein
                            float(line.split()[2]))
                elif kwargs["dir"] == "z":
                    self.membrane_complex.complex.prot_z = \
                        float(line.split()[3]) # hprotein
                    self.set_box_sizes()

                break

        return True

    def set_stage_init(self, **kwargs):
        '''Copy a set of files from the source to the target dir'''
        if not os.path.isdir(kwargs["tgt_dir"]): os.mkdir(kwargs["tgt_dir"])

        for src_file in kwargs["src_files"]:
            if(os.path.isfile(os.path.join(kwargs["src_dir"], src_file))):
                shutil.copy(os.path.join(kwargs["src_dir"], src_file),
                            os.path.join(kwargs["tgt_dir"], src_file))

        if "repo_files" in kwargs.keys():
            for repo_file in kwargs["repo_files"]:
                shutil.copy(os.path.join(self.repo_dir, repo_file),
                            os.path.join(kwargs["tgt_dir"], repo_file))

        return True

    def set_steep(self, **kwargs):
        '''Copy the template steep.mdp to the target dir'''
        shutil.copy(os.path.join(self.repo_dir, "steep.mdp"),
                    "steep.mdp")
        return True

    def set_water(self, **kwargs):
        '''Creates a water layer for a box'''
        start = (self.membrane_complex.membrane.bilayer_zw - \
                 self.membrane_complex.complex.prot_z) / 2
        end = start + self.membrane_complex.complex.prot_z

        src = open(self.membrane_complex.membrane.pdb, "r")
        tgt = open(kwargs["tgt"], "w")

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

class Wrapper(object):
    def __init__(self, *args, **kwargs):
        self.work_dir = os.getcwd()
        #The gromacs to be used
        self.gromacs_dir = settings.GROMACS_PATH
        #And this is our directory to refer the "fixed" files
        self.own_dir = os.path.dirname(os.path.abspath(__file__))
        #Repo dir is under pymoldyn file directory
        self.repo_dir = os.path.join(self.own_dir, "templates")

    def _common_io(self, src, tgt):
        '''Autoexpand many Gromacs commands that uses -f for the input
        and -o for the output file'''
        return ["-f", src, "-o", tgt]

    def generate_command(self, kwargs):
        '''Received some variables in kwargs, generate the appropiate command 
        to be run. Return a set in the form of a string "command -with flags"'''
        try:
            mode = kwargs["gromacs"]
        except KeyError:
            raise

        if "src" in kwargs["options"].keys():
            src = self._setDir(kwargs["options"]["src"])
        if "tgt" in kwargs["options"].keys():
            tgt = self._setDir(kwargs["options"]["tgt"])
        options = kwargs["options"]

        command = [os.path.join(self.gromacs_dir, mode)]
        if "queue" in kwargs.keys():
            if hasattr(kwargs["queue"], mode):
               # If we got a queue enabled for this command, use it
               command = list(kwargs["queue"].command) # Already a list
               kwargs["queue"].make_script(
                   workdir = kwargs["options"]["dir"],
                   options = self._mode_mdrun(options))
            
        # Standard -f input -o output
        if mode in ["pdb2gmx", "editconf", "grompp", "trjconv",
                    "make_ndx", "genrestr", "g_energy"]:
            command.extend(self._common_io(src, tgt))

            if (mode == "pdb2gmx"): #PDB2GMX
                command.extend(self._mode_pdb2gmx(options))
            if (mode == "editconf"): #EDITCONF
                command.extend(self._mode_editconf(options))
            if (mode == "grompp"): #GROMPP
                command.extend(self._mode_grompp(options))
            if (mode == "trjconv"): #TRJCONV
                command.extend(self._mode_trjconv(options))
            if (mode == "make_ndx"): #MAKE_NDX
                pass
            if (mode == "genrestr"): #GENRSTR
                command.extend(self._mode_genrest(options))
            if (mode == "g_energy"): #G_ENERGY
                pass

        else:
            if (mode == "genbox"): #GENBOX
                command.extend(self._mode_genbox(options))
            if (mode == "genion"): #GENION
                command.extend(self._mode_genion(options))
            if (mode == "tpbconv"): #TPBCONV
                command.extend(self._mode_tpbconv(options))
            if (mode == "trjcat"): #TRJCAT
                command.extend(self._mode_trjcat(options)) 
            if (mode == "eneconv"): #ENECONV
                command.extend(self._mode_eneconv(options))
            if (mode == "mdrun"): #MDRUN_SLURM
                pass
                #command.extend(self._mode_mdrun(options))

        return command

    def _mode_editconf(self, kwargs):
        '''Wraps the editconf command options'''
        command = []

        if "dist" in kwargs.keys():
            command.extend(["-d", str(kwargs["dist"])])
        if "box" in kwargs.keys():
            command.extend(["-box"])
            command.extend(kwargs["box"])
        if "angles" in kwargs.keys():
            command.extend(["-angles"])
            command.extend(kwargs["angles"])
            command.extend(["-bt", kwargs["bt"]])
        if "translate" in kwargs.keys():
            command.extend(["-translate"])
            command.extend(kwargs["translate"])

        return command

    def _mode_eneconv(self, kwargs):
        '''Wraps the eneconv command options'''
        src_files = utils.make_cat(kwargs["dir1"],
                                   kwargs["dir2"],
                                   kwargs["name"])
        command = ["-f"]
        command.extend(src_files)
        command.extend(["-o", self._setDir(kwargs["tgt"]),
                        "-settime"])

        return command

    def _mode_genbox(self, kwargs):
        '''Wraps the genbox command options'''
        return ["-cp", self._setDir(kwargs["cp"]),
                "-cs", self._setDir(kwargs["cs"]),
                "-p", self._setDir(kwargs["top"]),
                "-o", self._setDir(kwargs["tgt"])]

    def _mode_genion(self, kwargs):
        '''Wraps the genion command options'''
        command = ["-s", kwargs["src"],
                   "-o", kwargs["tgt"],
                   "-p", self._setDir(kwargs["src2"]),
                   "-g", self._setDir("genion.log"),
                   "-np", str(kwargs["np"]),
                   "-nn", str(kwargs["nn"]),
                   "-pname", "NA+",
                   "-nname", "CL-"]

        return command

    def _mode_genrest(self, kwargs):
        '''Wraps the genrest command options'''
        command = ["-fc"] + kwargs["forces"]

        if "index" in kwargs.keys():
            command.extend(["-n", kwargs["index"]])        
        
        return command

    def _mode_grompp(self, kwargs):
        '''Wraps the grompp command options'''
        command = ["-c", self._setDir(kwargs["src2"]),
                   "-p", self._setDir(kwargs["top"]),
                   "-po", self._setDir("mdout.mdp")]
        if "index" in kwargs.keys():
            command.extend(["-n", self._setDir(kwargs["index"])])

        return command

    def _mode_mdrun(self, kwargs):
        '''Wraps the mdrun command options'''

        command = ["-s", kwargs["src"],
                   "-o", kwargs["tgt"],
                   "-e", kwargs["energy"],
                   "-c", kwargs["conf"],
                   "-g", kwargs["log"]]

        if("traj" in kwargs.keys()):
            command.extend(["-x", kwargs["traj"]])

        if("cpi" in kwargs.keys()):
            command.extend(["-cpi", kwargs["cpi"]])

        return command

    def _mode_pdb2gmx(self, kwargs):
        '''Wraps the pdb2gmx command options'''
        return ["-p", self._setDir(kwargs["top"]),
                "-i", self._setDir("posre.itp"),
                "-ignh", "-ff", "oplsaa"]

    def _mode_tpbconv(self, kwargs):
        '''Wraps the tpbconv command options'''
        return ["-s", self._setDir(kwargs["src"]),
                "-o", self._setDir(kwargs["tgt"]),
                "-extend", kwargs["extend"]]

    def _mode_trjcat(self, kwargs):
        '''Wraps the trjcat command options'''
        src_files = utils.make_cat(kwargs["dir1"],
                                   kwargs["dir2"],
                                   kwargs["name"])
        command = ["-f"]
        command.extend(src_files)
        command.extend(["-o", self._setDir(kwargs["tgt"]),
                        "-settime"])

        return command

    def _mode_trjconv(self, kwargs):
        '''Wraps the trjconf command options'''
        command = ["-s", self._setDir(kwargs["src2"]),
                   "-pbc", kwargs["pbc"]]
        if "ur" in kwargs.keys():
            command.extend(["-ur", kwargs["ur"]])
        if "trans" in kwargs.keys():
            command.extend(["-trans"])
            command.extend([str(x) for x in kwargs["trans"]])
        else:
            command.extend(["-center"])

        return command

    def run_command(self, kwargs):
        '''Run a command that comes in kwargs in a subprocess, and returns
        the output as (output, errors)'''

        command = self.generate_command(kwargs)

        #my_dir = os.getcwd()

        #if kwargs["gromacs"] == "mdrun":
            #MDRUN depends on the local .mdp, no option to set it.
        #    os.chdir(kwargs["options"]["dir"])

        if("input" in kwargs.keys()):
            p = subprocess.Popen(command,
                stdin = subprocess.PIPE,
                stdout = subprocess.PIPE,
                stderr = subprocess.PIPE)
            gro_out, gro_errs = p.communicate(kwargs["input"])
        else:
            p = subprocess.Popen(command,
                stdout = subprocess.PIPE,
                stderr = subprocess.PIPE)

            gro_out, gro_errs = p.communicate()

        #os.chdir(my_dir)

        return gro_out, gro_errs

    def _setDir(self, filename):
        '''Expand a filename with the work dir, just to save code space'''
        return os.path.join(self.work_dir, filename)
