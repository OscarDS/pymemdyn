import os

class Queue(object):
    def __init__(self, *args, **kwargs):
        self.num_proc = 8 #Default number of processors to be used
        self.max_time = "50:00:00"
        self.sh = "mdrun.sh"

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

        self._mdrun = "/opt/gromacs405/bin/mdrun"

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

        self._mdrun = "/opt/gromacs405/bin/mdrun_slurm"

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

class PBS(Queue):
    def __init__(self, *args, **kwargs):
        super(PBS, self).__init__(self, *args, **kwargs)
        '''Setting the command to run mdrun in pbs queue with mpi'''
        self.num_proc = 40
        self.max_time = "08:00:00"
        self.max_cpu_time = "320:00:00"
        self.max_mem = "12gb"
        self.sh = "mdrun.sh"

        self._mdrun="mdrun_mpi"

        self.command = ["qsub",
            "-nodes=%d:ppn=4:ib" % self.num_proc,
            "-walltime=%s" % self.max_time,
            "-cput=%s" % self.max_cpu_time,
            "-mem=%s" % self.max_mem,
            self.sh]
    
    def make_script(self, workdir, options):
        '''PBS must load some modules in each node by shell scripts
        options is a list with all the options'''
        sh = open(self.sh, "w")
        sh.write("#!/bin/bash\n")
        sh.write("cd $PBS_O_WORKDIR\n")
        sh.write("module load gromacs\n")
        sh.write("mpirun %s %s\n" % (self.mdrun, " ".join(options)))
        sh.close()
        os.chmod(self.sh, 0755)
        
        return True

class Other(Queue):
    def __init__(self):
        pass
