import os

import gromacs


class Membrane(object):
    """
    Set the characteristics of the membrane in the complex.
    """
    def __init__(self, *args, **kwargs):
        if "pdb" not in kwargs.keys():
            self.pdb = os.path.join(
                gromacs.Wrapper().repo_dir, "x4bilayer.pdb")
        else:
            self.pdb = kwargs["pdb"]

        self.bilayer_x = 61.842  # xbilayer
        self.bilayer_y = 61.842  # ybilayer
        self.bilayer_zw = 154.053  # hbilayerw i.e. with waters
        self.bilayer_z = 57.78   # zbilayer i.e. without waters
        self.set_nanom()

    def set_nanom(self):
        """
        Convert some measurements to nanometers to comply with GROMACS units.
        """
        nanometer = 10
        self.gmx_bilayer_x = self.bilayer_x / nanometer
        self.gmx_bilayer_y = self.bilayer_y / nanometer
        self.gmx_bilayer_z = self.bilayer_z / (nanometer - 1) #Loose the lipids
        self.gmx_bilayer_zw = self.bilayer_zw / nanometer

        return True