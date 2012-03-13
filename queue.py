class Queue(object):
    def __init__(self, *args, **kwargs):
        self.num_proc = 8 #Default number of processors to be used
        self.max_time = "50:00:00"

class Slurm(Queue):
    def __init__(self, *args, **kwargs):
        super(Slurm, self).__init__(self, *args, **kwargs)
        self.command = ["srun", "-n", str(self.num_proc), "-t", self.max_time]
        self._mdrun = "mdrun_slurm"

    def set_mdrun(self, value):
        '''Sets the md_run command'''
        self._mdrun = value
    def get_mdrun(self):
        return self._mdrun
    mdrun = property(get_mdrun, set_mdrun)

class Other(Queue):
    def __init__(self):
        pass
