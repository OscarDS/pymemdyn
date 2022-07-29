class MembraneComplex(object):
    def __init__(self, *args, **kwargs):
        self.box_height = 3.5  # Min xy
        self.box_width  = 1.2

    def setMembrane(self, membrane):
        """
        Set the membrane pdb file
        """
        self.membrane = membrane
        
    def getMembrane(self):
        return self.membrane
    property(getMembrane, setMembrane)

    def setComplex(self, complex):
        """
        Set the complex object
        """
        self.complex = complex
        
    def getComplex(self):
        return self.complex
    property(getComplex, setComplex)