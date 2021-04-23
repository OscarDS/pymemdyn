import os

BINDIR = "/home/gpcradmin/gpcr_modsim/gpcr_backend/"

#This dir now has to be the absolute path to the source of templates
# You can refer it relatively like this:
TEMPLATES_DIR = os.path.join(BINDIR, "templates")
# Or absolutely like this:
# TEMPLATES_DIR = "/path/to/your/templates"
# But this WILL FAIL: * TEMPLATES_DIR = "templates" *

# Choose a path to the gromacs binaries.
GROMACS_PATH = "/home/apps/apps/gromacs/4.6.7/bin/"                      #csb.bmc.uu.se
#GROMACS_PATH = "/home/apps/gromacs-4.6.5/bin/"                      #csb.bmc.uu.se

# Define a path to the clustalw binary.
#CLUSTAL_BIN = os.path.join(ROOT_DIR, ".bin/clustalw_mac")
CLUSTAL_BIN = os.path.join(BINDIR, ".bin/clustalw_linux")

# Choose which queuing system to use. Look inside queue.py.
QUEUE = ""
#QUEUE = "slurm"
#QUEUE = "pbs"
#QUEUE = "pbs_ib"
#QUEUE = "svgd"

# Choose how many nodes to use in parallel
QUEUE_NUM_NODES = 1

# Choose how many processor to use in parallel
QUEUE_NUM_PROCS = 8

# Choose the maximum alloted time for your run.
QUEUE_MAX_TIME = "47:59:00"

