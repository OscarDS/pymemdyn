import argparse
import logging
import os
import textwrap

import complex
import gromacs
import membrane
import protein
import queue
import recipes
import settings

class Run(object):
    def __init__(self, *args, **kwargs):
        self.own_dir = kwargs.own_dir or ""
        self.repo_dir = kwargs.repo_dir or ""
        self.pdb = kwargs.pdb or ""
        self.ligand = kwargs.ligand or ""
        self.queue = kwargs.queue or ""
        self.debug = kwargs.debug or False

        if self.pdb: protein.Monomer(pdb = self.pdb)
        if self.ligand:
            protein.Ligand(pdb = self.ligand + ".pdb",
                           itp = self.ligand + ".itp")

        membr = membrane.Membrane()

        prot_complex = protein.ProteinComplex(
            monomer = self.monomer,
            ligand = self.ligand)
        full_complex = complex.MembraneComplex()
        full_complex.complex = prot_complex
        full_complex.membrane = membr

        self.g = gromacs.Gromacs(membrane_complex = full_complex)
        if self.queue:
            if self.queue == "slurm":
                my_queue = queue.Slurm()
            self.g.queue = my_queue

    def moldyn(self):
        '''Runs all the dynamics'''
        self.g.run_recipe()
        if self.ligand:
            self.g.recipe = recipes.LigandMinimization()
        else:
            self.g.recipe = recipes.BasicMinimization()
        self.g.run_recipe()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description = textwrap.dedent('''\
    == This script runs a Molecular Dynamic with a PDB. ==
    '''))

    parser.add_argument('-b', dest = "own_dir",
                        help = "Working dir if different from actual dir",
                        default = os.getcwd())
    parser.add_argument('-r', dest = "repo_dir",
                        help = "Path to templates of fixed files",
                        default = settings.REPO_DIR)
    parser.add_argument('-p',
                        dest = "pdb",
                        required = True,
                        help = "Name of the pdb to insert into MD (mandatory)")
    parser.add_argument('-l', dest = "ligand",
                        help = "Name of the ligand, without extension. Two \
                                files must be present along with the molecule \
                                pdb: the ligand and its itp.")
    parser.add_argument('-q', dest = "queue",
                        help = "Queue system to use (SLURM supported)",
                        default = settings.QUEUE)
    parser.add_argument('--debug', action="store_true")
    args = parser.parse_args()

    if args.debug:
        logging.basicConfig(filename='md_file.log',level=logging.DEBUG)
    else:
        logging.basicConfig(filename='md_file.log',
                            format='%(asctime)s %(message)s',
                            datefmt='%m/%d/%Y %I:%M:%S',
                            level=logging.ERROR)

    run = Run(own_dir = args.own_dir,
              repo_dir = args.repo_dir,
              pdb = args.pdb,
              ligand = args.ligand,
              queue = args.queue,
              debug = args.debug)
    run.moldyn()

