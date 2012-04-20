class ProteinComplex(object):
    def __init__(self, *args, **kwargs):
        self.cres = 0 #Central residue
        self.trans = [0, 0, self.cres] #Module for translating complex
        self.n_wats = 0 #Number of experimental waters

        if "monomer" in kwargs.keys():
            self.setMonomer(kwargs["monomer"])
        if "ligand" in kwargs.keys():
            self.setLigand(kwargs["ligand"])
        if "waters" in kwargs.keys():
            self.setWaters(kwargs["waters"])
        if "ions" in kwargs.keys():
            self.setIons(kwargs["ions"])
        if "cho" in kwargs.keys():
            self.setCho(kwargs["cho"])

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

    def setWaters(self, value):
        '''Sets the crystal waters object'''
        self.waters = value
    def getWaters(self):
        return self.waters
    property(getWaters, setWaters)

    def setIons(self, value):
        '''Sets the ions object'''
        self.ions = value
    def getIons(self):
        return self.ions
    property(getIons, setIons)

    def setCho(self, value):
        '''Sets the cholesterol object'''
        self.cho = value
    def getCho(self):
        return self.cho
    property(getCho, setCho)


    def set_nanom(self):
        '''Set some meassurements to nanometers, as GROMACS wants'''
        NANOM = 10
        self.gmx_prot_xy = self.prot_xy / NANOM
        self.gmx_prot_z = self.prot_z / NANOM

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

class CrystalWaters(object):
    def __init__(self, *args, **kwargs):
        self.pdb = "hoh.pdb"
        self.itp = "hoh.itp"
        self.posre_itp = "posre_hoh.itp"
        self._n_wats = self.count_waters()

    def setWaters(self, value):
        '''Sets the crystal waters'''
        self._n_wats = value
    def getWaters(self):
        '''Get the crystal waters'''
        return self._n_wats
    number = property(getWaters, setWaters)

    def count_waters(self):
       '''Count and set the number of waters in the pdb'''
       return len([x for x in open(self.pdb, "r") if "OW" in x])

    def _setITP(self):
        '''Create the itp to this structure'''
        s = "\n".join([
            "; position restraints for crystallographic waters (resn HOH)",
            "[ position_restraints ]",
            ";  i funct       fcx        fcy        fcz",
            "   1    1       1000       1000       1000"])

        tgt = open(self.posre_itp, "w")
        tgt.writelines(s)
        tgt.close()

class Ions(object):
    def __init__(self, *args, **kwargs):
        self.pdb = "local_ions.pdb"
        self.itp = "local_ions.itp"

        self._n_ions = self.count_ions()

    def setIons(self, value):
        '''Sets the crystal ions'''
        self._n_ions = value
    def getIons(self):
        '''Get the crystal ions'''
        return self._n_ions
    number = property(getIons, setIons)

    def count_ions(self):
       '''Count and set the number of ions in the pdb'''
       ions = ["NA", "CA", "MG", "CL", "ZN"]
       ion_count = 0
       for line in open(self.pdb, "r"):
           if len(line.split()) > 2:
               if line.split()[2] in ions:
                   ion_count += 1
       return ion_count

class Cholesterol(object):
    def __init__(self, *args, **kwargs):
        self.pdb = "cho.pdb"
        self.itp = "cho.itp"
