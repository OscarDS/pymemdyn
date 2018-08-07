import os

#PID_FILE = "/home/gpcradmin/gpcr_modsim/run/PyroDaemon.pid"
#LOG_FILE = "/home/gpcradmin/gpcr_modsim/log/PyroDaemon.log"
#MODEL_LOG_FILE = "/home/gpcradmin/gpcr_modsim/log/model_file.log"
#PYTHON_BIN = "/home/apps/bin/python2.7"
#JAVA_PATH = "/usr/java/latest/bin/java"
BINDIR = "/home/gpcradmin/gpcr_modsim/gpcr_backend/"
#TEMPLATES_DIR = "/home/gpcradmin/gpcr_modsim/templates/"

# Used in pyro_objects.py
#SLURMOUT = "/home/gpcradmin/gpcr_modsim/slurm-out/"
#PYMEMPATH = "/home/gpcradmin/gpcr_modsim/pymemdyn/"

# STAMP related fields (used in stamp.py)
#STAMP_BIN = "/home/apps/stamp.4.4/bin/linux/"
#STAMP_DIR =  "/home/apps/stamp.4.4/defs/"
#STAMP_TMPDIR = "/home/gpcradmin/tmp/"
#CLUSTAL_BIN = "/home/apps/clustalw"
#TEMPLATE_PDB = "/home/gpcradmin/gpcr_modsim/templates/3eml.pdb"

# MOLPROBITY
#MOLPROBITY_PATH = "/home/apps/molprobity/lib/"


#SCRIPT_DIR = os.path.dirname(os.path.realpath(pymoldyn.__file__))

# This is the dir where pymoldyn git repo has been deployed,
# or to be more specific, this file settings.py is located
#ROOT_DIR = os.path.dirname(os.path.realpath(__file__))

#This dir now has to be the absolute path to the source of templates
# You can refer it relatively like this:
TEMPLATES_DIR = os.path.join(BINDIR, "templates")
# Or absolutely like this:
# TEMPLATES_DIR = "/path/to/your/templates"
# But this WILL FAIL: * TEMPLATES_DIR = "templates" *

# Choose a path to the gromacs binaries.
#GROMACS_PATH = "/opt/applications/gromacs/4.0.5/gnu/ib/bin/"
#GROMACS_PATH = "/opt/applications/gromacs/4.0.5/gnu/gige/bin/"
#GROMACS_PATH = "/opt/cesga/gromacs-4.0.7/bin/"
#GROMACS_PATH = "/opt/gromacs405/bin/"                               #cuelebre.inv.usc.es
GROMACS_PATH = "/home/apps/gromacs/4.6.7/bin/"                      #csb.bmc.uu.se
#GROMACS_PATH = "/home/apps/gromacs-4.6.5/bin/"                      #csb.bmc.uu.se
#GROMACS_PATH = "/software/apps/gromacs/4.6.3/g472/bin/"             #Triolith
#GROMACS_PATH = "/sw/bin/"                                           #Standalone in Mac Fink
#GROMACS_PATH = "/Users/esguerra/software/gromacs-4.6.5/bin/"        #Standalone in Mac
#GROMACS_PATH = "/c3se/apps/Glenn/gromacs/4.6.3-p20130821-gcc48/bin" #Glenn at Chalmers
#GROMACS_PATH = "/sw/apps/gromacs/4.6.3/tintin/bin"                  #Tintin
#GROMACS_PATH = "/lap/gromacs/4.6.5/bin"                             #Abisko

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

