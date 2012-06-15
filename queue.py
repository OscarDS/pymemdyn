import os
import settings

class Queue(object):
    def __init__(self, *args, **kwargs):
        self.num_proc = 8 #Default number of processors to be used
        self.max_time = "50:00:00"
        self.sh = "./mdrun.sh"

    def set_mdrun(self, value):
        '''Sets the md_run command'''
        self._mdrun = value
    def get_mdrun(self):
        return self._mdrun
    mdrun = property(get_mdrun, set_mdrun)

class NoQueue(Queue):
    '''Dummy queue when no queue is selected'''
    def __init__(self, *args, **kwargs):
        super(NoQueue, self).__init__(self, *args, **kwargs)
        self.command = [self.sh]

        self._mdrun = os.path.join(settings.GROMACS_PATH, "mdrun")

    def make_script(self, workdir, options):
        '''binary is the executable
        options is a list with all the options'''
        sh = open(self.sh, "w")
        sh.write("#!/bin/bash\n")
        sh.write("cd %s\n" % workdir)
        sh.write("%s %s\n" % (self.mdrun, " ".join(options)))
        sh.close()
        os.chmod(self.sh, 0755)

        return True

class Slurm(Queue):
    def __init__(self, *args, **kwargs):
        super(Slurm, self).__init__(self, *args, **kwargs)
        self.command = ["srun",
            "-n", str(self.num_proc),
            "-t", self.max_time,
            self.sh]

        self._mdrun = os.path.join(settings.GROMACS_PATH, "mdrun_slurm")

    def make_script(self, workdir, options):
        '''binary is the executable
        options is a list with all the options'''
        sh = open(self.sh, "w")
        sh.write("#!/bin/bash\n")
        sh.write("cd %s\n" % workdir)
        sh.write("%s %s -v&>mdrun.log\n" % (self.mdrun, " ".join(options)))
        sh.close()
        os.chmod(self.sh, 0755)

        return True

class PBS(Queue):
    '''Queue for the PBS system'''
    def __init__(self, *args, **kwargs):
        super(PBS, self).__init__(self, *args, **kwargs)
        '''Setting the command to run mdrun in pbs queue with mpi'''
        # These values are here for reference, doesn't do NOTHING      #
        # Calling file run.sh should resemble this lines               #
        self.num_nodes = 5                                             #
        self.proc_per_node = 8                                         #
        self.max_time = "36:00:00"                                     #
        self.max_cpu_time = "1440:00:00"                               #
        self.max_mem = "12gb"                                          #
        self.command = ["qsub",                                        #
            "-nodes=%d:ppn=%d" % (self.num_nodes, self.proc_per_node), #
            "-walltime=%s" % self.max_time,                            #
            "-cput=%s" % self.max_cpu_time,                            #
            "-mem=%s" % self.max_mem,                                  #
            self.sh]                                                   #
                                                                       #
        ################################################################

        self._mdrun=os.path.join(settings.GROMACS_PATH, "mdrun_mpi")
        self.command = [self.sh]

    def make_script(self, workdir, options):
        '''PBS must load some modules in each node by shell scripts
        options is a list with all the options'''
        sh = open(self.sh, "w")
        sh.write("#!/bin/bash\n")
        sh.write("cd %s\n" % os.path.join(os.getcwd(), workdir))
        sh.write("module load gromacs/4.0.5-gige\n")
        sh.write("mpirun %s %s -v &>mdrun.log\n" \
            % (self.mdrun, " ".join(options)))
        sh.close()
        os.chmod(self.sh, 0755)

        return True

class PBS_IB(Queue):
    def __init__(self, *args, **kwargs):
        super(PBS, self).__init__(self, *args, **kwargs)
        # USELESS, see class PBS for explanation                            #
        self.num_nodes = 10                                                 #
        self.proc_per_node = 4                                              #
        self.max_time = "08:00:00"                                          #
        self.max_cpu_time = "320:00:00"                                     #
        self.max_mem = "12gb"                                               #
                                                                            #
        self.command = [                                                    #
            "qsub",                                                         #
            "-l nodes=%d:ppn=%d:ib" % (self.num_nodes, self.proc_per_node), #
            "-l walltime=%s" % self.max_time,                               #
            "-l cput=%s" % self.max_cpu_time,                               #
            "-l mem=%s" % self.max_mem,                                     #
            self.sh]                                                        #
        #####################################################################

        self._mdrun=os.path.join(settings.GROMACS_PATH, "mdrun_mpi")
        self.command = [self.sh]

    def make_script(self, workdir, options):
        '''PBS must load some modules in each node by shell scripts
        options is a list with all the options'''
        sh = open(self.sh, "w")
        sh.write("#!/bin/bash\n")
        sh.write("cd %s\n" % os.path.join(os.getcwd(), workdir))
        sh.write("module load gromacs\n")
        sh.write("mpirun %s %s -v &>mdrun.log\n" \
            % (self.mdrun, " ".join(options)))
        sh.close()
        os.chmod(self.sh, 0755)

        return True

class Svgd(Queue):
   '''Queue for the PBS system at svgd.cesga.es'''
    def __init__(self, *args, **kwargs):
        super(PBS, self).__init__(self, *args, **kwargs)
        '''Setting the command to run mdrun in pbs queue with mpi'''
        self._mdrun=os.path.join(settings.GROMACS_PATH, "mdrun")
        self.command = [self.sh]

    def make_script(self, workdir, options):
        '''PBS must load some modules in each node by shell scripts
        options is a list with all the options'''
        sh = open(self.sh, "w")
        sh.write("#!/bin/bash\n")
        sh.write("cd %s\n" % os.path.join(os.getcwd(), workdir))
        sh.write("module load acml\n")
        sh.write("module load gromacs/4.0.7\n")
        sh.write("%s %s -v &>mdrun.log\n" \
            % (self.mdrun, " ".join(options)))
        sh.close()
        os.chmod(self.sh, 0755)

        return True

class Other(Queue):
    def __init__(self):
        pass
