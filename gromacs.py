import utils
import os

class Gromacs(object):
    def __init__(self, *args, **kwargs):
        self.work_dir = os.getcwd() # This is where the library is called
        #And this is our directory to refer the "fixed" files
        self.own_dir = os.path.dirname(os.path.abspath(__file__))
        #Repo dir is under gromacs.py file directory
        self.repo_dir = os.path.join(self.own_dir, "templates")
        #The gromacs to be used
        self.gromacs_dir = "/opt/gromacs405/bin/"

        self.membrane_complex = None
        if "membrane_complex" in kwargs.keys():
            self.set_membrane_complex(kwargs["membrane_complex"])

    def set_membrane_complex(self, value):
        '''Sets the monomer object'''
        self.membrane_complex = value
    def get_membrane_complex(self):
        return self.membrane_complex
    property(get_membrane_complex, set_membrane_complex)

    def make_hexagon(self):
        
        if hasattr(self.membrane_complex.complex, "ligand"):
            #We got a ligand included
            from recipes import MonomerRecipe as recipe
        else:
            from recipes import MonomerLigandRecipe as recipe
            #We got no ligand

        r = recipe()
        return r.recipe

    def _set_popc(self, tgt = ""):
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

    def set_protein_size(self, src):
        '''Get the protein max base wide from a pdb file'''
        for line in open(src, "r"):
            if line.startswith("CRYST1"):
                if dir == "xy":
                    self.membrane_complex.complex.prot_xy = \
                        max(float(line.split()[1]), #xyprotein
                            float(line.split()[2]))
                elif dir == "z":
                    self.membrane_complex.complex.prot_z = \
                        float(line.split()[3]) # hprotein
                break

        return True

    def test_wrapper(self):
        w = Wrapper()
        return w.own_dir

class Wrapper(Gromacs):
    def __init__(self, *args, **kwargs):
        super(Wrapper, self).__init__(*args, **kwargs)

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

        # STD -f input -o output
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
                   "-np", kwargs["np"],
                   "-nn", kwargs["nn"],
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
            command.extend(kwargs["trans"])
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
