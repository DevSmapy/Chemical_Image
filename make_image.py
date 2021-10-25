from rdkit.Chem import Draw
from rdkit import Chem
import matplotlib.pyplot as plt
import pybel
class molto2D:
    def __init__(self,inf):
        self._inf = inf
    def pmol(self):
        # Make Mol by using pybel

        # 1. Check SMILES #
        try:
            mol = pybel.readstring("smi",self._inf)
        except:
            print("It is not SMILES string")
            pass
        # 2. Check PDB #
        try:
            mol = pybel.readstring("smi",list(pybel.readfile("pdb",self._inf))[0].write("can").split()[0])
        except:
            print("It is not PDB file")
            pass
        # 3. Check SDF #
        try:
            mol = pybel.readstring("smi",list(pybel.readfile("sdf",self._inf))[0].write("can").split()[0])
        except:
            print("It is not SDF file")
            return 0
        return mol
    def rmol(self):
        # Make Mol by RD-kit

        # 1. Check SMILES #
        try:
            mol = Chem.MolFromSmiles(self._inf)
        except:
            print("It is not SMILES string")
            pass
        # 2. Check PDB #
        try:
            mol = Chem.MolFromSmiles(Chem.MolToSmiles(Chem.MolFromPDBFile(self._inf)))
        except:
            print("It is not PDB file")
            pass
        # 3. Check SDF #
        try:
            mol = Chem.MolFromSmiles(Chem.MolToSmiles(next(Chem.SDMolSupplier(self._inf))))
        except:
            print("It is not SDF file")
            return 0
        return mol

    ###########################
    # Make Image using RD-kit #
    ###########################
    def rmol_img_mem(self):
        mol = self.rmol()
        try:
            img = Draw.MolToImage(mol,size=(800,800),kekulize=True,fitImage=True)
        except:
            return 0
        return img
    def rmol_img_file(self,fname,fPath):
        img = self.rmol_img_mem()
        # customizing figure #
        plt.figure(figsize=(5,5))
        plt.imshow(img)
        plt.tight_layout()
        plt.axis("off")
        # notation #
        plt.text(0,700,"* %s"%fname.strip(),fontsize=12.0,fontstyle="oblique")
        plt.text(0,750,"%s"%self._inf.strip(),fontsize=10.0)

        if fPath is None or fPath == "":
            fPath = "./"
        else:
            pass
        #plt.show()
        # save image #
        plt.savefig(fPath +"/" +  fname + ".png",dpi=100)
    ##########################
    # Make Image using Babel #
    ##########################
    def pmol_img_file1(self,fname,fPath):
        mol = self.pmol()
        try:
            mol.draw(show=False,filename=fPath + "/" + fname + ".raw.png")
        except:
            return 0
        try:
            img = plt.imread(fPath + "/" + fname + ".raw.png")
        except:
            return 0
        return img
    def pmol_img_file2(self,fname,fPath):
        img = self.pmol_img_file1(fname,fPath)
        # customizing figure
        plt.figure(figsize=(5,5))
        plt.imshow(img)
        plt.tight_layout()
        plt.axis("off")
        # notation
        plt.text(0,225,"* %s"%fname.strip(),fontsize=12.0,fontstyle="oblique")
        plt.text(0,250,"%s"%self._inf.strip(),fontsize=10.0)

        if fPath is None or fPath == "":
            fPath = "./"
        else:
            pass
        #plt.show()
        # save image
        plt.savefig(fPath + "/" + fname + ".png",dpi=100)
#############
# Code Test #
#############
if __name__ == "__main__":
    smi = "./example/example_chemical.sdf" #"C3c5cc(Oc1nc2ccccc2(cc1))ccc5(OCC3Cc4cnccc4)"
    pp = molto2D(smi)
    pp.rmol_img_file("rmol_sdf","./outputs")