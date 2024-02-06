"""In this file multiple checks are defined for the protein.
"""
import os 
import math 

import logging
check_logger = logging.getLogger('pymemdyn.checks')

try:
    import modeller 
    from modeller import automodel
except:
    check_logger.warning("""!! WARNING !! : No installation of MODELLER was found.
        Missing loops and/or sidechains cannot be remodelled and will cause errors.""")


from aminoAcids import AminoAcids
from Bio.PDB import PDBParser, PDBIO, Select
from collections import Counter

import logging
check_logger = logging.getLogger('pymemdyn.checks')

class CheckProtein():
    def __init__(self, *args, **kwargs):
        self.pdb = kwargs['pdb']
        self.chains = kwargs['chains']
        self.loop_fill = kwargs['loop_fill']
        self.logger = logging.getLogger('pymemdyn.checks.checkProtein')
        self.aa = AminoAcids()
        self.nr_missing = 0
        
    def find_missingLoops(self):
        """Check if the residue numbering is continuous, write sequence.
        If loops are missing, return missingLoc.
        """
        with open(self.pdb, 'r') as inf:
            lines = inf.readlines()

        # Initialize params
        missingLoc = {}
        resID_prev = 0
        prev_chain = ''
        first_res = True
        first_chain_res = False
        pdbseq = ''
        aa_d = self.aa.codes321 # 3-letter code to 1 letter code dict is now aa_d

        # Loop over protein.pdb
        for i, line in enumerate(lines):
            splitted = line.split() # Careful here, sometimes columns are stuck together!
            
            # Determine if line is useful
            if not splitted[0] == 'ATOM':
                if splitted[0] == 'TER': # End of chain
                    first_chain_res = True
                    self.logger.debug('found TER at line {}'.format(i))
                    pdbseq += '/'
                continue # Skip to next line
            
            else: # Useful line
                try: 
                    resID = int(line[22:26])

                    # Save ID of the first residue
                    if first_res:
                        self.first_res_ID = int(line[22:26])
                        self.logger.debug('first resID: {}'.format(self.first_res_ID))
                        resID_prev = resID - 1
                        first_res = False
                        
                    chainID = line[21]
                    
                except:
                    raise Exception("Cannot read resID or chain ID in the following line: \n{}".format(line))
                
                # Check for non-std res
                if line[17:20] not in list(aa_d.keys()):
                    self.logger.exception("Residue {}, {} is not a standard amino acid. Non-standard amino acids are not supported.".format(resID, splitted[3]))
                    raise Exception("Residue {}, {} is not a standard amino acid. Non-standard amino acids are not supported.".format(resID, splitted[3]))
                
                resIDchain = chainID + str(resID)

                if (resID != resID_prev) and (resID != resID_prev + 1): # missing
                    if first_chain_res:
                        first_chain_res = False
                        self.logger.debug('Ignoring residue {} for missing residue check(first residue of chain)'.format(resID))
                    else:
                        missingLoc[prev_chain+str(resID_prev)] = resIDchain

                        self.nr_missing = resID - resID_prev -1
                        pdbseq += '-'*self.nr_missing # for missing residue
                
                if resID != resID_prev: # next res
                    pdbseq += aa_d[line[17:20]] # 1 letter code



                resID_prev = resID
                prev_chain = chainID
        
        # End of for-loop                
        self.logger.debug('length missingLoc: {}'.format(len(missingLoc)))
        
        # Generate fasta file (not used, but nice for debugging)
        seq = open("seq_original.fasta", "w") 
        seq.write(">Generated automatically from {}\n".format(self.pdb))   
        seq.write(pdbseq)
        
        self.seq = pdbseq
        
        return missingLoc

    def find_missingSideChains(self):
        self.structure = PDBParser(QUIET=True).get_structure('protein', self.pdb)

        missing_sideChains = []
        residues_to_exclude = []
        residues_to_keep = []

        for model in self.structure:
            for chain in model:
                prev_residue = None
                for residue in chain:
                    atom_count = Counter()
                    current_atoms = [atom.element for atom in residue if atom.element != 'H']  # Exclude hydrogen
                    for atom in current_atoms:
                        atom_count[atom] += 1

                    amino = residue.resname   
                    self.logger.debug('Residue {}: {} atom counts ({}) | reference atom counts ({}).'.format(prev_residue, amino, self.aa.sideChains[amino], atom_count))              
                    if self.aa.sideChains[amino] != atom_count:
                        self.logger.info('Residue {}: {} atom counts ({}) does not match with reference atom counts ({}). \
                                         It will be deleted from your pdb and replaced with MODELLER'.format(prev_residue, amino, self.aa.sideChains[amino], atom_count))
                        missing_sideChains.append((prev_residue, amino))
                        residues_to_exclude.append(residue)
                    else:
                        residues_to_keep.append(residue)

        class ResiduesSelect(Select):
            def __init__(self, residues):
                self.residues = residues

            """A custom selector that excludes specific residues."""
            def accept_residue(self, residue):
                return residue in self.residues

        io = PDBIO()
        io.set_structure(self.structure)
        io.save(self.pdb, ResiduesSelect(residues_to_keep))

        return missing_sideChains

    def make_ml_pir(self, **kwargs):
        """
        make_ml_pr: Modify missing regions and create a MODELLER alignment file (.pir)

        :kwarg work_dir: working directory
        :kwarg tgt1: alignment.pir
        """
        with open(os.path.join(kwargs["work_dir"], self.pdb), 'r') as in_pdb:
            lines_pdb = in_pdb.readlines()

        missingLoc = self.find_missingLoops()
        if len(missingLoc) == 0:
            self.logger.info('No missing loops found')
        else:
            self.logger.info('Missing loops found in locations {}. Will be replaced with poly-ala.'.format(list(missingLoc.keys())))
        
        missing_SideChains = self.find_missingSideChains()
 
        if len(missing_SideChains) == len(missingLoc) == 0:
            return False

        pdbseq = self.seq
        self.logger.debug('pdbseq: {}'.format(pdbseq))

        if pdbseq.endswith("/"):
            pdbseq = pdbseq[:-1]  # looses trailing '/'

        mod_seq = pdbseq
        tmpl_seq = pdbseq

        for sc in missing_SideChains:
            loc = sc[0] - self.first_res_ID
            code = self.aa.codes321[sc[1]]
            tmpl_seq = tmpl_seq[:loc] + '-' + tmpl_seq[loc+1:]
            mod_seq = mod_seq[:loc] + code + mod_seq[loc+1:]

        target = ''

        self.logger.debug(list(missingLoc.keys()))

        first_res_ID = self.first_res_ID

        for loopstart in list(missingLoc.keys()):
            loop = int(loopstart[1:])
            chain = loopstart[0]
            for line in lines_pdb: 
                if line[:4] == "ATOM":    
                    if int(line[22:26]) == loop and \
                        "".join(line[11:17]).strip() == "C" and line[21] == chain:

                        start = [float(line[30:38]), float(line[38:46]), float(line[46:54])]
                        target = int(missingLoc[loopstart][1:])
                        self.logger.debug('loop: {}\n'.format(loop))
                    if int(line[22:26]) == target and \
                        "".join(line[11:17]).strip() == "N" and line[21] == chain:

                        end = [float(line[30:38]), float(line[38:46]), float(line[46:54])]
                        self.logger.debug('target: {}\n'.format(target))

            aa_dist = float(self.loop_fill)           
            x = abs(start[0] - end[0])
            y = abs(start[1] - end[1])
            z = abs(start[2] - end[2])
            dist = math.sqrt(x*x + y*y + z*z)
            num_aa = math.ceil(dist / aa_dist)     
            self.logger.debug('\tnum_aa: {}\n'.format(num_aa))

            # Determine which chain we are in (and how many '/'-signs we need to ignore)
            ignore_chain_characters = self.chains.index(chain)
            
            start_seq_index = loop - first_res_ID +1 - ignore_chain_characters  # (we ignored the beginloop)
            end_seq_index = target - first_res_ID +1 - ignore_chain_characters
            self.logger.debug('\tstart: {}'.format(start_seq_index))
            self.logger.debug('\tend:   {}\n'.format(end_seq_index))

            loop_seq_dash = pdbseq[start_seq_index:end_seq_index-1]   # '----'
            before = pdbseq[start_seq_index-2:start_seq_index]        # 'AB'
            after = pdbseq[end_seq_index-1:end_seq_index+1]           # 'BA'
            loop_seq = before + loop_seq_dash + after                 # 'AB----BA'
            self.logger.debug('\tbefore: {}, dash: {}, after: {}'.format(before, loop_seq_dash, after))
            self.logger.debug('\tloop_seq: {}\n'.format(loop_seq))

            gaps = before + "A" * num_aa + after        # 'ABAABA'
            loop_replace = before + num_aa *'-' + after # 'AB--BA'
            
            self.logger.debug('\tgaps:     {}\n'.format(gaps))
            self.logger.debug('\tloop_rep: {}\n'.format(loop_replace))
            
            mod_seq = mod_seq.replace(loop_seq, gaps)
            tmpl_seq = tmpl_seq.replace(loop_seq, loop_replace)
            
        mod_seq += "*"
        tmpl_seq += "*"

        self.logger.debug(mod_seq)
        self.logger.debug(tmpl_seq)

        self.logger.debug("these lengths should be equal: {}, {}".format(len(mod_seq), len(tmpl_seq)))


        # Write to .pir file

        # tgt is the alignment.pir file
        tgt1 = open(os.path.join(kwargs['work_dir'], kwargs["tgt1"]), "w") 
        
        tgt1.write(f'>P1;{self.pdb}\n')
        
        firstResID_tmpl = int(first_res_ID) + (len(tmpl_seq) - len(tmpl_seq.lstrip('-')))
        finalResID = int((len(pdbseq) + 1) / len(self.chains))
        tgt1.write(f'structure:{self.pdb}:{firstResID_tmpl}:{self.chains[0]}:{finalResID+first_res_ID-1}:{self.chains[-1]}::: 2.10:-1.00\n')

        new_line = 0
        for aa in tmpl_seq:
            tgt1.write(aa)            
            new_line += 1
            if (new_line % 60) == 0:
                tgt1.write("\n")

        tgt1.write("\n>P1;refined\n")
        tgt1.write("sequence:::::::::\n")

        new_line = 0
        for aa in mod_seq:
            tgt1.write(aa)  
            new_line += 1
            if (new_line % 60) == 0:
                tgt1.write("\n")

        tgt1.close()
        return True

    def refine_protein(self, **kwargs):
        """
        Refine protein structure using MODELLER

        :kwarg knowns: the .pdb file to be refined
        :returns: .pdb file of refined protein
        """
        first_res = self.first_res_ID
        self.logger.debug('self.chains: {}{}'.format(self.chains, type(self.chains[0])))
        chains = self.chains

        class MyModel(automodel.AutoModel):
            def special_patches(self, aln):
                # Renumber residues in chains
                self.rename_segments(segment_ids=chains,      
                                    renumber_residues=[first_res, first_res, first_res, first_res])

        env = modeller.Environ()
        env.io.atom_files_directory = ['.', '../atom_files']

        a = MyModel(env,
                    alnfile  = 'alignment.pir',  
                    knowns   = kwargs["knowns"],     
                    sequence = 'refined')

        a.starting_model= 1              
        a.ending_model  = 1              
        a.md_level = None                 

        a.make()

        # Get new file name
        refined_pdb = a.get_model_filename(root_name='refined', id1=9999, id2=1, file_ext='.pdb')

        
        return refined_pdb
