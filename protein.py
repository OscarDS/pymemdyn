class ProteinComplex(object):
    def __init__(self, *args, **kwargs):
        self.cres = 0 #Central residue
        self.n_wats = 0 #Number of experimental waters

        if "monomer" in kwargs.keys():
            self.setMonomer(kwargs["monomer"])
        if "ligand" in kwargs.keys():
            self.setLigand(kwargs["ligand"])

    def setMonomer(self, value):
        '''Sets the monomer object'''
        self.monomer = value
    def getMonomer(self):
        return self.monomer
    property(getMonomer, setMonomer)

    def setLigand(self, value):
        '''Sets the ligand object'''
        self.ligand = value
    def getLigand(self):
        return self.ligand
    property(getLigand, setLigand)

class Monomer(object):
    def __init__(self, *args, **kwargs):
        self.pdb = kwargs["pdb"]
        self._setHist()

    def _setHist(self):
        '''Touch the Histidine in the source protein, generating a new PDB'''
        tgt = open(self.pdb.replace(".pdb", "-his.pdb"), "w")
        self.pdb_his = tgt.name

        for line in open(self.pdb, "r"):
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

class Dimer(Monomer):
    def __init__(self):
        super(Dimer, self).__init__(self, *args, **kwargs)

class Ligand(object):
    def __init__(self, *args, **kwargs):
        self.pdb = kwargs["pdb"]
        self.itp = kwargs["itp"]
