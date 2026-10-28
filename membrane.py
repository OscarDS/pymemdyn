import os

import forcefield_defaults
import gromacs


class Membrane(object):
    """
    Set the characteristics of the membrane in the complex.
    """
    def __init__(self, *args, **kwargs):
        lipidforcefield = kwargs.get("lipidforcefield") or "berger"
        lipid_template = forcefield_defaults.LIPIDFF_TEMPLATES.get(lipidforcefield.lower())
        if lipid_template:
            repo_dir = lipid_template["templates_dir"]
        else:
            repo_dir = kwargs.get("repo_dir") or gromacs.Wrapper().repo_dir

        # Select bilayer PDB based on lipid forcefield
        bilayer_map = {
            "berger": "x4bilayer.pdb",
            "amberlipid21": "lipid21bilayer.pdb"
        }
        bilayer_file = bilayer_map.get(lipidforcefield.lower(), "x4bilayer.pdb")
        
        if "pdb" not in kwargs.keys():
            self.pdb = os.path.join(repo_dir, bilayer_file)
        else:
            self.pdb = kwargs["pdb"]

        #self.bilayer_x = 61.842  # xbilayer berger
        #self.bilayer_y = 61.842  # ybilayer berger
        #self.bilayer_zw = 154.053  # hbilayerw i.e. with waters berger
        #self.bilayer_z = 57.78   # zbilayer i.e. without waters (max - min coordinate of bilayer in the Z axis) Berger. Actual membrane thickness: 39.71
        
        self.bilayer_x = 64.329  # xbilayer amber 63.908
        self.bilayer_y = 64.329  # ybilayer amber 63.908
        self.bilayer_zw = 153.242  # hbilayerw i.e. with waters amber 147.557
        self.bilayer_z = 56.36   # zbilayer i.e. without waters amber. Membrane thickness: 38.76
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