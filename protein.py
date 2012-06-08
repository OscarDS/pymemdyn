import os

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
        if not os.path.isfile(self.pdb):
            raise IOError("File '{0}' missing".format(self.pdb))
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

class Compound(object):
    '''This is a super-class to provide common functions to added compounds'''
    def __init__(self, *args, **kwargs):
        self.check_files(self.pdb, self.itp)

    def check_files(self, *files):
        '''Check if files passed as *args exist'''
        for src in files:
            if not os.path.isfile(src):
                raise IOError("File {0} missing".format(src))

class Ligand(Compound):
    def __init__(self, *args, **kwargs):
        self.pdb = kwargs["pdb"]
        self.itp = kwargs["itp"]
        super(Ligand, self).__init__(self, *args, **kwargs)

        self.force_field = kwargs["ff"]

        self.check_forces()

    def check_forces(self):
        '''A force field must give a set of forces that matches every atom in
        the pdb file. This showed particularly important to the ligands, as they
        may vary along a very broad range of atoms'''

        #The itp matches each residue in the ligand pdb with the force field
        atoms_def = False
        molecules = {}
        for line in open(self.itp, "r"):
            if "[ atoms ]" in line:
                atoms_def = True
            if "[ bonds ]" in line:
                atoms_def = False
            if atoms_def and not line.startswith(";"):
                data = line.split()
                if len(data) > 6:
                    if data[3] not in molecules.keys(): molecules[data[3]] = {}
                    #{"LIG": {"C1": "TC1"},}
                    molecules[data[3]][data[4]] = data[1]

        atoms = {}
        #The force field matches each atom in the pdb with one line
        for line in open(self.force_field, "r"):
            if not line.startswith(";"):
                if (len(line.split()) > 6):
                    #{"TC1": "C1"}
                    atoms[line.split()[0]] = line.split()[1]

        #The pdb have the name of the atom in the third position.
        #Here we cross-check all three files to match their harvested values
        for line in open(self.pdb, "r"):
            data = line.split()
            if len(data) > 6:
                if molecules[data[3]][data[2]] not in atoms.keys():
                    #Some atoms in the pdb has no definition in force field
                    # TODO : add a guessing function
                    print "Atom {0} have no field definition".format(data[1])
                    #return False
                if atoms[molecules[data[3]][data[2]]] not in\
                    molecules[data[3]].keys():
                    print "Atom {0} have a wrong field definition".format(
                        data[1])
                    #return False

        return True

class CrystalWaters(Compound):
    def __init__(self, *args, **kwargs):
        self.pdb = "hoh.pdb"
        self.itp = "hoh.itp"
        super(CrystalWaters, self).__init__(self, *args, **kwargs)

        self.posre_itp = "posre_hoh.itp"
        self._setITP()
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

class Ions(Compound):
    def __init__(self, *args, **kwargs):
        self.pdb = "ions_local.pdb"
        self.itp = "ions_local.itp"
        super(Ions, self).__init__(self, *args, **kwargs)

        self.posre_itp = "posre_ion.itp"
        self._setITP()

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

    def _setITP(self):
        '''Create the itp to this structure'''
        s = "\n".join([
            "; position restraints for ions (resn HOH)",
            "[ position_restraints ]",
            ";  i funct       fcx        fcy        fcz",
            "   1    1       1000       1000       1000"])

        tgt = open(self.posre_itp, "w")
        tgt.writelines(s)
        tgt.close()

class Cholesterol(Compound):
    def __init__(self, *args, **kwargs):
        self.pdb = "cho.pdb"
        self.itp = "cho.itp"
        super(Cholesterol, self).__init__(self, *args, **kwargs)

        #self.posre_itp = "posre_cho.itp"
        #self._setITP()

        self._n_cho = self.count_cho()

    def setCho(self, value):
        '''Sets the crystal cholesterol'''
        self._n_cho = value
    def getCho(self):
        '''Get the crystal cholesterols'''
        return self._n_cho
    number = property(getCho, setCho)

    def count_cho(self):
       '''Count and set the number of cho in the pdb'''
       cho_count = 0
       for line in open(self.pdb, "r"):
           if len(line.split()) > 2:
               if line.split()[3] == "CHO":
                   cho_count += 1
       return cho_count / 74 #Each CHO has 74 atoms

    def _setITP(self):
        '''Create the itp to this structure'''
        s = "\n".join([
            "; position restraints for cholesterol (resn CHO)",
            "[ position_restraints ]",
            ";  i funct       fcx        fcy        fcz",
            "   1    1       1000       1000       1000"])

        tgt = open(self.posre_itp, "w")
        tgt.writelines(s)
        tgt.close()

class Lipids(Compound):
    def __init__(self, *args, **kwargs):
        self.pdb = "lip.pdb"
        self.itp = "lip.itp"
        super(Lipids, self).__init__(self, *args, **kwargs)

        #self.posre_itp = "posre_lip.itp"
        #self._setITP()

        self._n_lip = self.count_lip()

    def setLip(self, value):
        '''Sets the crystal lipids'''
        self._n_lip = value
    def getLip(self):
        '''Get the crystal lipidss'''
        return self._n_lip
    number = property(getLip, setLip)

    def count_lip(self):
       '''Count and set the number of lipids in the pdb'''
       lip_count = []
       for line in open(self.pdb, "r"):
           if len(line.split()) > 2:
               if (line.split()[3] == "LIP" and
                   line.split()[4].isdigit() and
                   line.split()[4] not in lip_count):
                   #Lipid + number + new
                   lip_count.append(line.split()[4])
       return len(lip_count)

    def _setITP(self):
        '''Create the itp to this structure'''
        s = "\n".join([
            "; position restraints for cholesterol (resn CHO)",
            "[ position_restraints ]",
            ";  i funct       fcx        fcy        fcz",
            "   1    1       1000       1000       1000"])

        tgt = open(self.posre_itp, "w")
        tgt.writelines(s)
        tgt.close()

class Alosteric(Compound):
    def __init__(self, *args, **kwargs):
        self.pdb = kwargs["pdb"]
        self.itp = kwargs["itp"]
        super(Alosteric, self).__init__(self, *args, **kwargs)

        self.force_field = kwargs["ff"]

