import argparse
import logging
import os
import shutil
import sys
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
        self.own_dir = kwargs["own_dir"] or ""
        self.repo_dir = kwargs["repo_dir"] or ""
        self.pdb = kwargs["pdb"] or ""
        self.ligand = kwargs["ligand"] or ""
        self.waters = kwargs["waters"] or False
        self.ions = kwargs["ions"] or False
        self.cho = kwargs["cho"] or False
        self.queue = kwargs["queue"] or ""
        self.debug = kwargs["debug"] or False

        if self.pdb:
            self.pdb = protein.Monomer(pdb = self.pdb)
        if self.ligand:
            self.ligand = protein.Ligand(pdb = self.ligand + ".pdb",
                                         itp = self.ligand + ".itp")

        if self.waters:
            self.waters = protein.CrystalWaters()

        if self.ions:
            self.ions = protein.Ions()

        if self.cho:
            self.cho = protein.Cholesterol()

        self.membr = membrane.Membrane()

        prot_complex = protein.ProteinComplex(
            monomer = self.pdb,
            ligand = self.ligand or None,
            waters = self.waters or None,
            ions = self.ions or None,
            cho = self.cho or None)
        full_complex = complex.MembraneComplex()
        full_complex.complex = prot_complex
        full_complex.membrane = self.membr

        self.g = gromacs.Gromacs(membrane_complex = full_complex)
        if self.queue:
            if self.queue == "slurm":
                my_queue = queue.Slurm()
            self.g.queue = my_queue

    def clean(self):
        '''Removes all previously generated files'''
        to_unlink = ["#index.ndx.1#", "#ligand_ha.ndx.1#", "#mdout.mdp.1#",
            "#mdout.mdp.2#", "#mdout.mdp.3#", "#output.pdb.1#",
            "#output.tpr.1#", "#popc.pdb.1#",
            "#posre.itp.1#", "#proteinopls.pdb.1#", "#proteinopls.pdb.2#",
            "#proteinopls.pdb.3#", "#proteinopls.pdb.4#","#protein.top.1#",
            "#protpopc.pdb.1#", "#protpopc.pdb.2#", "#tmp.pdb.1#",
            "#topol.top.1#", "#topol.top.2#", "#topol.top.3#",
            "#topol.tpr.1#", "#topol.tpr.2#", "#topol.tpr.3#",
            "#topol.tpr.4#", "#water.pdb.1#", "protein_ca200.itp",
            "ffoplsaabon_mod.itp", "ffoplsaa_mod.itp", "ffoplsaanb_mod.itp",
            "GROMACS.output", "genion.log", "hexagon.pdb", "index.ndx",
            "ligand_ha.ndx", "mdout.mdp", "min.pdb", "output.pdb",
            "output.tpr", "popc.pdb", "popc.itp", "posre.itp",
            "posre_lig.itp", "protein.itp", "protein.top",
            "protein_ca200.itp", "proteinopls.pdb", "proteinopls-ligand.pdb",
            "protpopc.pdb", "steep.mdp", "traj.xtc", "tmp.pdb", "topol.top",
            "topol.tpr", "tmp_proteinopls.pdb", "Y1_min-his.pdb", "water.pdb"]

        dirs_to_unlink = ["Rmin", "eq", "eqCA"]

        for target in to_unlink:
            if os.path.isfile(target): os.unlink(target)

        for target in dirs_to_unlink:
            if os.path.isdir(target): shutil.rmtree(target)
    
        return True

    def moldyn(self):
        '''Runs all the dynamics'''
        self.g.run_recipe(debug = self.debug) #Basic recipe

        self.g.recipe = recipes.LigandMinimization(debug = self.debug)
        self.g.run_recipe()
        if self.ligand:
            self.g.recipe = recipes.LigandEquilibration(debug = self.debug)
        else:
            self.g.recipe = recipes.BasicEquilibration(debug = self.debug)
        self.g.run_recipe()

        if self.ligand:
            self.g.recipe = recipes.LigandRelax(debug = self.debug)
        else:
            self.g.recipe = recipes.BasicRelax(debug = self.debug)
        self.g.run_recipe()
        self.g.recipe = recipes.CAEquilibrate(debug = self.debug)
        self.g.run_recipe()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description = textwrap.dedent('''\
    == This script runs a Molecular Dynamic with a PDB. ==
    '''))

    parser.add_argument('-b',
                        dest = "own_dir",
                        help = "Working dir if different from actual dir",
                        default = os.getcwd())
    parser.add_argument('-r',
                        dest = "repo_dir",
                        help = "Path to templates of fixed files",
                        default = settings.REPO_DIR)
    parser.add_argument('-p',
                        dest = "pdb",
                        required = True,
                        help = "Name of the pdb to insert into MD (mandatory)")
    parser.add_argument('-l',
                        dest = "ligand",
                        help = "Name of the ligand, without extension. Two \
                                files must be present along with the molecule \
                                pdb: the ligand and its itp.")
    parser.add_argument('--waters',
                        action="store_true",
                        help = "Crystalized water molecules hoh.pdb file \
                                present.")
    parser.add_argument('--ions',
                        action="store_true",
                        help = "Crystalized ions ions_local.pdb and \
                                ions_loca.itp file present.")
    parser.add_argument('--cho',
                        action="store_true",
                        help = "Crystalized cholesterol molecules cho.pdb \
                                file present.")
    parser.add_argument('-q',
                        dest = "queue",
                        help = "Queue system to use (SLURM supported)",
                        default = settings.QUEUE)
    parser.add_argument('--debug',
                        action="store_true")
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
              waters = args.waters,
              ions = args.ions,
              cho = args.cho,
              queue = args.queue,
              debug = args.debug)
    run.clean()
    run.moldyn()

