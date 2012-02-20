class ProteinComplex:
    def __init__(self):
        self.cres = 0 #Central residue
        self.n_wats = 0 #Number of experimental waters

    def _setCharge(self):
        '''Calculate the total charge of the compound'''
        self.charge = 0
        positives = ["HIP", "LYS", "ARG"]
        negatives = ["ASP", "GLU"]

        for line in open(self.pdb_file, "r"):
            if (len(line.split()) > 3 and line.split()[0] == "ATOM"):
                if line.split()[2] == "CA":
                    if line.split()[3] in positives:
                        self.charge += 1
                    elif line.split()[3] in negatives:
                        self.charge -= 1
        return True

    def _setHist(self):
        '''Touch the Histidine in the source protein, generating a new PDB'''
        tgt = open(self.pdb.replace(".pdb", "-his.pdb"), "w")

        for line in open(self.pdb_file, "r"):
            if len(line.split()) > 3:
                if line.split()[3] == "HIE":
                    tgt.write(line.replace('HIE ','HISB'))
                elif line.split()[3] == "HID":
                    tgt.write(line.replace('HID ','HISA'))
                elif line.split()[3] == "HIP":
                    tgt.write(line.replace('HIP ','HISH'))
                else:
                    tgt.write(line)
            else:
                tgt.write(line)
        tgt.close()

        return True

class Monomer(ProteinComplex):
    def __init__(self):
        super(Monomer, self).__init__()
        self.pdb = "" #The ABSOLUTE path to the pdb file
        self.itp = ""

class Dimer(ProteinComplex):
    def __init__(self):
        super(Dimer, self).__init__()

class Ligand(ProteinComplex):
    def __init__(self):
        super(Ligand, self).__init__()
        self.pdb = "" #The ABSOLUTE path to the pdb file
        self.itp = ""
