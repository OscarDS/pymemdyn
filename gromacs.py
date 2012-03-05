import recipes
import utils

import logging
import os
import shutil
import subprocess
import sys

class Gromacs(object):
    def __init__(self, *args, **kwargs):
        self.membrane_complex = None
        if "membrane_complex" in kwargs.keys():
            self.set_membrane_complex(kwargs["membrane_complex"])

    def set_membrane_complex(self, value):
        '''Sets the monomer object'''
        self.membrane_complex = value
    def get_membrane_complex(self):
        return self.membrane_complex
    property(get_membrane_complex, set_membrane_complex)

    def count_lipids(self, **kwargs):
        '''Count the lipids in source and write a target with N4 tags'''
        src = open(kwargs["src"], "r")
        tgt = open(kwargs["tgt"], "w")

        half = self.membrane_complex.complex.prot_z / 2
        self.membrane_complex.membrane.lipids_up = 0
        self.membrane_complex.membrane.lipids_down = 0
        self.membrane_complex.membrane.n_wats = 0

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
        wrapper = Wrapper()

        out, err = wrapper.run_command({"gromacs": "grompp",
                                        "options": kwargs})

        # Now we are looking for this line:
        # System has non-zero total charge: 6.000002e+00
        for line in err.split("\n"):
            if "total charge" in line:
                charge = int(float(line.split()[-1]))
                break

        self.membrane_complex.complex.negative_charge = 0
        self.membrane_complex.complex.positive_charge = 0
        if charge > 0:
            self.membrane_complex.complex.positive_charge = charge
        elif charge < 0:
            self.membrane_complex.complex.positive_charge = 0

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

    def run_recipe(self):
        '''Run the recipe for the complex'''
        if not hasattr(self, "recipe"):
            self.select_recipe()

        wrapper = Wrapper()
        self.repo_dir = wrapper.repo_dir

        for n, command in enumerate(self.recipe.recipe):
            if n in self.recipe.breaks.keys():
                command["options"] = self.set_options(command["options"],
                                     self.recipe.breaks[n])

            # NOW RUN IT !
            print n, command
            if command.has_key("gromacs"):
                # Either run a Gromacs pure command...
                out, err = wrapper.run_command(command)
            else:
                # ...or run a local function
                try:
                    if command.has_key("options"):
                        getattr(self, command["command"])(**command["options"])
                    else:
                        getattr(self, command["command"])()
                except AttributeError:
                    #Fallback to the utils module
                    getattr(utils, command["command"])(**command["options"])

            if n == 30:
                print " ".join(wrapper.generate_command(command))
                print out, err
                #print self.membrane_complex.complex.gmx_prot_z
                sys.exit()

    def select_recipe(self):
        '''Select the appropiate recipe for the complex'''
        if self.membrane_complex:
            if hasattr(self.membrane_complex.complex, "ligand"):
                self.recipe = recipes.MonomerLigandRecipe()
            else:
                self.recipe = recipes.MonomerRecipe()

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
        itp = ""
        if hasattr(self.membrane_complex.complex, "ligand"): itp = "_lig"

        shutil.copy(os.path.join(
            self.repo_dir, "steep.mdp"),
            "steep.mdp")
        shutil.copy(os.path.join(
            self.repo_dir, "popc.itp"),
            "popc.itp")
        shutil.copy(os.path.join(
            self.repo_dir, "ffoplsaanb_mod{0}.itp".format(itp)),
            "ffoplsaanb_mod.itp")
        shutil.copy(os.path.join(
            self.repo_dir, "ffoplsaabon_mod{0}.itp".format(itp)),
            "ffoplsaabon_mod.itp")
        shutil.copy(os.path.join(
            self.repo_dir, "ffoplsaa_mod{0}.itp".format(itp)),
            "ffoplsaa_mod.itp")
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
        self.gromacs_dir = "/opt/gromacs405/bin/"
        #And this is our directory to refer the "fixed" files
        self.own_dir = os.path.dirname(os.path.abspath(__file__))
        #Repo dir is under gromacs.py file directory
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
            if (mode == "mdrun_slurm"): #MDRUN_SLURM
                command.extend(self._mode_mdrun_slurm(options))

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
        command = ["-f"]
        command.extend(kwargs["src_files"])
        command.extend(["-o", self._setDir(kwargs["tgt"]),
                        "-settime"])

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
        return ["-fc"] + kwargs["forces"]

    def _mode_grompp(self, kwargs):
        '''Wraps the grompp command options'''
        command = ["-c", self._setDir(kwargs["src2"]),
                   "-p", self._setDir(kwargs["top"]),
                   "-po", self._setDir("mdout.mdp")]
        if "index" in kwargs.keys():
            command.extend(["-n", self._setDir(kwargs["index"])])

        return command

    def _mode_mdrun_slurm(self, kwargs):
        '''Wraps the mdrun_slurm command options'''
        command = ["-s", self._setDir(kwargs["src"]),
                   "-o", self._setDir(kwargs["tgt"]),
                   "-e", self._setDir(kwargs["energy"]),
                   "-c", self._setDir(kwargs["conf"]),
                   "-g", self._setDir(kwargs["log"])]

        if("traj" in kwargs.keys()):
            command.extend(["-x", self._setDir(kwargs["traj"])])

        if("cpi" in kwargs.keys()):
            command.extend(["-cpi", self._setDir(kwargs["cpi"])])

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

        if("multiproc" in kwargs.keys()):
            #Command is (probably) a mdrun, so it requires parallelization
            #TODO: This initial command should be generated inside "Queue"
            # class
            command = ["srun", "-n", kwargs["multiproc"], "-t", "50:00:00"]
        else:
            command = []

        command.extend(self.generate_command(kwargs))

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

        return gro_out, gro_errs

    def _setDir(self, filename):
        '''Expand a filename with the work dir, just to save code space'''
        return os.path.join(self.work_dir, filename)
