import protein

import unittest

class TestMonomerProtein(unittest.TestCase):

    def setUp(self):
        self.pdb = "tests/monomer.pdb"

    def test_monomer_load(self):
        '''Check the correct loading of a simple monomer
        '''
        monomer = protein.Monomer(pdb = self.pdb)


if __name__ == "__main__":
    unittest.main()
