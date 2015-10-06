import os

# This is the folder where pymemdyn git repo has been deployed,
# or to be more specific, where the settings.py file is located at.
ROOT_DIR = os.path.dirname(os.path.realpath(__file__))

# The following variable can be given as a full path to the place where
# the templates reside or it can also use a relative path to the previously
# defined ROOT_DIR.
TEMPLATES_DIR = os.path.join(ROOT_DIR, "templates")
# Or using the full path like so:
# TEMPLATES_DIR = "/path/to/your/templates"

# Define a path to the gromacs binaries.
#GROMACS_PATH = "/opt/applications/gromacs/4.0.5/gnu/ib/bin/"
#GROMACS_PATH = "/opt/applications/gromacs/4.0.5/gnu/gige/bin/"
#GROMACS_PATH = "/opt/cesga/gromacs-4.0.7/bin/"
#GROMACS_PATH = "/opt/gromacs405/bin/"                               #cuelebre.inv.usc.es
#GROMACS_PATH = "/home/apps/gromacs-4.6.7/bin/"                      #csb.bmc.uu.se
#GROMACS_PATH = "/software/apps/gromacs/4.6.3/g472/bin/"             #Triolith
#GROMACS_PATH = "/sw/bin/"                                           #Standalone in Mac Fink
GROMACS_PATH = "/Users/esguerra/software/gromacs-4.6.5/bin/"        #Standalone in Mac
#GROMACS_PATH = "/c3se/apps/Glenn/gromacs/4.6.3-p20130821-gcc48/bin" #Glenn at Chalmers
#GROMACS_PATH = "/c3se/apps/Glenn/gromacs/5.0.4-gcc48-cuda/bin/"     #Glenn GPU at Chalmers
#GROMACS_PATH = "/sw/apps/gromacs/4.6.3/tintin/bin"                  #Tintin
#GROMACS_PATH = "/lap/gromacs/4.6.5/bin"                             #Abisko

# Define a path to the clustalw binary.
CLUSTAL_BIN = os.path.join(ROOT_DIR, ".bin/clustalw_linux")

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
QUEUE_MAX_TIME = "60:00:00"

QUEUE_NUM_TASK = 8
QUEUE_NTS_NODE = 8
