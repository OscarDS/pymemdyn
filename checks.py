"""In this file multiple checks are defined for the protein.
"""
import os 
import math 

import modeller 
from modeller import automodel

from aminoAcids import AminoAcids
import logging
check_logger = logging.getLogger('pymemdyn.checks')

class CheckProtein():
    def __init__(self, *args, **kwargs):
        self.pdb = kwargs['pdb']
        self.chains = kwargs['chains']
        self.loop_fill = kwargs['loop_fill']
        self.logger = logging.getLogger('pymemdyn.checks.checkProtein')
        
    def find_missingLoops(self, tgt):
        """Check if the residue numbering is continuous.
        If loops are missing, return missingLoc.
        """
        with open(self.pdb, 'r') as inf:
            lines = inf.readlines()

        seq_dict = {}
        missingLoc = {}
        resID_prev = 0
        prev_line = "ATOM      0  N   XXX     0     000.000 000.000 000.000  0.00000.00           N" # Dummy line
        first_res = True

        aa = AminoAcids()
        aa_d = aa.codes321 # 3-letter code to 1 letter code dict is now aa_d

       
        for i, line in enumerate(lines):
            splitted = line.split() # Careful here, sometimes columns are stuck together!
            if not splitted[0] == 'ATOM':
                continue
            else:
                try: 
                    resID = int(line[22:26])
                    if first_res: # Save ID of the first residue
                        self.first_res_ID = int(line[22:26])
                        first_res = False
                    chainID = line[21]
                except:
                    raise Exception("Cannot read resID or chain ID in the following line: \n{}".format(line))
                if chainID not in list(seq_dict.keys()):
                    seq_dict[chainID] = {}
                
                if splitted[3] not in list(aa.codes321.keys()):
                    self.logger.exception("Residue {}, {} is not a standard amino acid. Non-standard amino acids are not supported.".format(resID, splitted[3]))
                    raise Exception("Residue {}, {} is not a standard amino acid. Non-standard amino acids are not supported.".format(resID, splitted[3]))
                
                seq_dict[chainID][resID] = splitted[3] # e.g. {A: {40: ASP}}
                
                resIDchain = chainID + str(resID)
                chain_prev = prev_line[21]

                if (resID != resID_prev) and (resID != resID_prev + 1):
                    if int(splitted[1]) == 1:
                        resID_prev = resID
                        prev_line = line
                        continue # ignore missing starting loop  

                    missingLoc[chain_prev+str(resID_prev)] = resIDchain
                    # for_modeller.append('{}:{}, {}:{}'.format(resID_prev, chain_prev, resID, chainID))

                resID_prev = resID
                prev_line = line
        self.logger.debug('length missingLoc: {}'.format(len(missingLoc)))

        # End of for-loop


        # # Create missingLoops.txt (Not used)
        # with open(tgt, 'w') as outf:
        #     outf.write("{:<6} |\n".format("resID"))
        #     outf.write("-"*8+"\n")
        #     for l in missingLoc.keys():
        #         outf.write("{:<6} |\n".format(l))
                
        #         k = missingLoc[l]
        #         outf.write("{:<6} |\n".format(k))
        #         outf.write("-"*8+"\n")

        

        # Generate 1 letter code sequence

        

        pdbseq = ''

        seq = open("seq_original.fasta", "w") # seq_original.fasta is not used, but self.seq is.
        seq.write(">Generated automatically from {}\n".format(self.pdb))   

        for c in list(seq_dict.keys()):
            ids = list(seq_dict[c].keys())
            ids.sort()
            self.logger.debug('chain {}:'.format(c))
            self.logger.debug('{} ids'.format(len(ids)))
            finalID = ids[-1]
            for i in range(self.first_res_ID, finalID+1):
                if i in ids:
                    assert seq_dict[c][i] in list(aa_d.keys()), 'non-standard amino acid found in residue {}: {}'.format(i, seq_dict[c][i])
                    seq.write(aa_d[seq_dict[c][i]])
                    pdbseq += aa_d[seq_dict[c][i]]

                else:
                    seq.write('-')
                    pdbseq += '-'
                    
            seq.write('/')
            pdbseq += '/'
        
        self.seq = pdbseq
        
        return missingLoc


    def make_ml_pir(self, **kwargs):
        """
        make_ml_pr: Modify missing regions and create a MODELLER alignment file (.pir)

        :param pdb: pdb file of protein
        :param tgt1: alignment.pir
        # :param tgt3: alignment.txt
        """
        with open(os.path.join(kwargs["work_dir"], self.pdb), 'r') as in_pdb:
            lines_pdb = in_pdb.readlines()

        missingLoc = self.find_missingLoops(tgt = 'missing_loops.txt')
        if len(missingLoc) == 0:
            self.logger.info('No missing loops found')
            return False
        else:
            self.logger.info('Missing loops found in locations {}. Will be replaced with poly-ala.'.format(list(missingLoc.keys())))


        pdbseq = self.seq
        pdbseq = pdbseq[:-1] # looses trailing '/'

        mod_seq = pdbseq
        tmpl_seq = pdbseq

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
            
            start_seq_index = loop - first_res_ID +1  # (we removed the beginloop)
            end_seq_index = target - first_res_ID +1

            loop_seq_dash = pdbseq[start_seq_index:end_seq_index-1]       # '----'
            before = pdbseq[start_seq_index-2:start_seq_index]                # 'AB'
            after = pdbseq[end_seq_index-1:end_seq_index+1]           # 'BA'
            loop_seq = before + loop_seq_dash + after   # 'AB----BA'
            
            # dashes_remain = len(loop_seq_dash) - num_aa
            gaps = before + "A" * num_aa + after        # 'ABAABA'
            loop_replace = before + num_aa *'-' + after # 'AB--BA'

            self.logger.debug('\tloop_seq: {}\n'.format(loop_seq))
            self.logger.debug('\tgaps:     {}\n'.format(gaps))
            
            mod_seq = mod_seq.replace(loop_seq, gaps)
            tmpl_seq = tmpl_seq.replace(loop_seq, loop_replace)
            # print('\t\tmodseq: \n\t\t{}\n\n'.format(mod_seq))
            # tmpl_seq = tmpl_seq[:loop[0]] + len(loop_mod) * "-" + tmpl_seq[loop[0]:]

            while mod_seq.startswith('-'):
                mod_seq = mod_seq[1:]
            
            while tmpl_seq.startswith('-'):
                tmpl_seq = tmpl_seq[1:]
            
            self.logger.debug('\tloop_rep: {}\n'.format(loop_replace))
            
        mod_seq += "*"
        tmpl_seq += "*"

        self.logger.debug("these lengths should be equal: {}, {}".format(len(mod_seq), len(tmpl_seq)))

        self.logger.debug(pdbseq)

        tgt1 = open(os.path.join(kwargs['work_dir'], kwargs["tgt1"]), "w") 
        # tgt3 = open(os.path.join(self.work_dir, kwargs["tgt3"]), "w")
        
        tgt1.write(f'>P1;{self.pdb}\n')
        finalResID = int((len(pdbseq) + 1) / len(self.chains))
        tgt1.write(f'structure:{self.pdb}:{first_res_ID}:{self.chains[0]}:{finalResID+first_res_ID-1}:{self.chains[-1]}::: 2.10:-1.00\n')

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

        refined_pdb = a.get_model_filename(root_name='refined', id1=9999, id2=1, file_ext='.pdb')

        
        return refined_pdb
