import shutil
import os
import numpy as np

import complex
import gromacs
import membrane
import protein
import queue

import logging
run_logger = logging.getLogger('pymemdyn.run')

class Run(object):
    def __init__(self, pdb, *args, **kwargs):
        """
        A molecular dynamics *Run()*  MUST be given a *pdb* file.

        This class tries to initialize a full complex to send to simulation.
        Given a set of molecules (protein, ligand, other ligand, waters, ...),
        this class would try to build a full embedded-in-membrane complex.

        The complex is stored in self.g (a *Gromacs* object), and thus
        can be **run** through g.recipe and g.run_recipe procedure. See
        gromacs.py for more information.

        The queueing system is also created here to be used in certain steps.
        """

        self.pdb = pdb
        self.own_dir = kwargs.get("own_dir") or ""
        self.repo_dir = kwargs.get("repo_dir") or ""
        self.ligand = kwargs.get("ligand") or ""
        self.ligand_charge = kwargs.get("ligand_charge") or ""
        self.protein = 'protein.pdb'
        self.waters = kwargs.get("waters") or ""
        self.ions = kwargs.get("ions") or ""
        self.restraint = kwargs.get("restraint") or ""
        self.queue = kwargs.get("queue") or ""
        self.debug = kwargs.get("debug") or False
        self.debugFast = kwargs.get("debugFast")
        self.logger = logging.getLogger('pymemdyn.run.Run')

        self.logger.debug('Run arguments initialized')

        self.logger.debug('self.ligand = '+str(self.ligand))
        self.logger.debug('self.ligand_charge = '+str(self.ligand_charge))


        if self.pdb:
            self.logger.debug('self.pdb dectected.')                            
            self.logger.debug('Splitting pdb file.')
            protein.System(pdb=self.pdb).split_system(ligand=self.ligand,
                                                      waters=self.waters,
                                                      ions=self.ions)

            self.logger.debug('Checking nr of chains of '+ str(self.protein))
            self.protein = protein.Protein(pdb=self.protein).check_number_of_chains()
            try:
                self.protein_center = protein.Protein(pdb=pdb).calculate_center()
            except:
                self.logger.warning("Cannot calculate center of protein. Please check alignment manually.")
            self.logger.debug('nr of chains ok!')

        sugars = {"ligand": "Ligand",
                  "allosteric": "Alosteric",
                  "waters": "CrystalWaters",
                  "ions": "Ions",
                  "cho": "Cholesterol"}
        self.logger.debug("dict 'sugars' initialized")

        self.logger.debug('starting sugar prep')

        protein.Sugar_prep.__init__(self)
        self.logger = logging.getLogger('pymemdyn.run.Run')

        for sugar_type, class_name in sugars.items():
            if getattr(self, sugar_type):
                base_name = getattr(self, sugar_type)
                setattr(self,
                        sugar_type,
                        getattr(protein, class_name)(
                            pdb=base_name + ".pdb",
                            itp=base_name + ".itp",
                            ff=base_name + ".ff"))

        self.membr = membrane.Membrane()

        if self.allosteric:
            nr_allosteric = self._n_alo
        else:
            nr_allosteric = None

        try:
            self.check_dist()
        except:
            self.logger.warning("Cannot check distance between protein and compound. \
                                Please check alignment manually.")
        

        prot_complex = protein.ProteinComplex(
            monomer=self.protein,
            ligand=self.ligand or None,
            allosteric=self.allosteric or None,
            nr_alo = nr_allosteric,
            waters=self.waters or None,
            ions=self.ions or None,
            cho=self.cho or None)

        full_complex = complex.MembraneComplex()

        full_complex.complex = prot_complex
        full_complex.membrane = self.membr

        self.g = gromacs.Gromacs(membrane_complex=full_complex)

        # NOTE: If not provided in command line, self.queue is set to
        # NoQueue
        if self.queue:
            if self.queue == "slurm":
                my_queue = queue.Slurm()
            elif self.queue == "pbs":
                my_queue = queue.PBS()
            elif self.queue == "pbs_ib":
                my_queue = queue.PBS_IB()
            elif self.queue == "svgd":
                my_queue = queue.Svgd()
        else:
            my_queue = queue.NoQueue()

        self.g.queue = my_queue

    def clean(self):
        """
        Removes all previously generated files
        """
        to_unlink = ["#index.ndx.1#", "#index.ndx.2#", "#index.ndx.3#", 
                     "#index.ndx.4#", "#index.ndx.5#", "#output.pdb.1#", 
                     "#proteinopls.pdb.1#", "#proteinopls.pdb.2#", 
                     "#proteinopls.pdb.3#", "#proteinopls.pdb.4#", 
                     "#protpopc.pdb.1#", "#protpopc.pdb.2#" ,"#tmp.pdb.1#", 
                     "#topol.top.1#", "#topol.top.2#", "#topol.top.3#", 
                     "#topol.tpr.1#", "#topol.tpr.2#", "#topol.tpr.3#", 
                     "#topol.tpr.4#", "disre.itp", "ener_EQ.edr", 
                     "ffoplsaa_mod.itp", "ffoplsaabon_mod.itp", 
                     "ffoplsaanb_mod.itp", "index.ndx", 
                     "ions.itp", "ligand_ha.ndx", "MD_output.tgz", "mdout.mdp", 
                     "mdrun.sh", "min.pdb", "output.pdb", "popc.gro", 
                     "popc.itp", "popc.pdb", "posre.itp", "posre_lig.itp", 
                     "posre_alo.itp","pressure.log", "pressure.xvg", 
                     "pressure2.log", "pressure2.xvg", "protein.itp", 
                     "protein.top", "protein_ca200.itp", "proteinopls.fasta", 
                     "proteinopls.pdb", "proteinopls_bw.aln",
                     "proteinopls_CA.pdb", "protpopc.pdb",
                     "rmsd-all-atom-vs-start.xvg", 
                     "rmsd-backbone-vs-start.xvg", "rmsd-calpha-vs-start.xvg", 
                     "rmsf-per-residue.xvg", "spc.itp", "steep.mdp", 
                     "temp.log", "temp.xvg", "temp2.log", "temp2.xvg", 
                     "tmp.pdb", "tmp_proteinopls.pdb", "topol.top", 
                     "topol.tpr", "tot_ener.log", "tot_ener.xvg", 
                     "tot_ener2.log", "tot_ener2.xvg", "traj_EQ.xtc", 
                     "traj_pymol.xtc", "volume.log", "volume.xvg", 
                     "volume2.log", "volume2.xvg", "water.gro", "water.pdb"]

        dirs_to_unlink = ["Rmin", "eq", "eqProd"]

        for target in to_unlink:
            if os.path.isfile(target): os.unlink(target)

        for target in dirs_to_unlink:
            if os.path.isdir(target): shutil.rmtree(target)

        return True

    def moldyn(self):
        """
        Run all steps in a molecular dynamics simulation of a membrane protein
        """
        if self.restraint == "bw":
            steps = ["Init", "Minimization", "Equilibration", "Relax", 
                     "BWRelax", "BWCollectResults"]
        elif self.restraint == "ca":
            steps = ["Init", "Minimization", "Equilibration", "Relax", 
                     "CARelax", "CACollectResults"]

        for step in steps:
            self.logger.info('\n\n[{}/{}]: {}\n'.format(steps.index(step)+1, len(steps), step))
            self.g.select_recipe(stage=step, debugFast=self.debugFast)
            self.g.run_recipe(debugFast=self.debugFast)

    def light_moldyn(self):
        """
        This is a function to debug a run in steps
        """
        steps = ["BWRelax", "BWCollectResults"]
#        steps = ["CACollectResults"]

        for step in steps:
            self.g.select_recipe(stage=step, debugFast=self.debugFast)
            self.g.run_recipe(debugFast = self.debugFast)

    def check_dist(self):
        """Check distance between protein and all possible ligands (if any).
        Raise warning is dist > 50
        """
        possible_ligands = [self.ligand, self.allosteric, self.ions, self.waters, self.cho]
        prot_centre = self.protein_center

        for thing in possible_ligands:
            if thing == "":
                continue
            else:
                d = np.linalg.norm(prot_centre - thing.center)
                if d > 45:
                    self.logger.warning("""Centre of {} and {} are unusually far apart. 
                                        Did you correctly align both?""".format(self.pdb.pdb, thing.pdb))
        return True

