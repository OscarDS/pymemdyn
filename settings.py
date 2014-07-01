#import pymoldyn
import os

#SCRIPT_DIR = os.path.dirname(os.path.realpath(pymoldyn.__file__))

# This is the dir where pymoldyn git repo has been deployed,
# or to be more specific, this file settings.py is located
ROOT_DIR = os.path.dirname(os.path.realpath(__file__))

#This dir now has to be the absolute path to the source of templates
# You can refer it relatively like this:
REPO_DIR = os.path.join(ROOT_DIR, "templates")
# Or absolutely like this:
# REPO_DIR = "/path/to/your/templates"
# But this WILL FAIL: * REPO_DIR = "templates" *

#GROMACS_PATH = "/opt/applications/gromacs/4.0.5/gnu/ib/bin/"
#GROMACS_PATH = "/opt/applications/gromacs/4.0.5/gnu/gige/bin/"
#GROMACS_PATH = "/opt/cesga/gromacs-4.0.7/bin/"
#GROMACS_PATH = "/opt/gromacs405/bin/"
#GROMACS_PATH = "/home/apps/gromacs-4.6.5/bin/"
#GROMACS_PATH = "/software/apps/gromacs/4.6.3/g472/bin/"
GROMACS_PATH = "/sw/bin/"

#Set the queue to use. Look inside queue.py.
QUEUE = ""
#QUEUE = "slurm"
#QUEUE = "pbs"
#QUEUE = "pbs_ib"
#QUEUE = "svgd"

QUEUE_NUM_PROCS = 8

QUEUE_MAX_TIME = "23:50:00"
