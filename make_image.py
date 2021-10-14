from rdkit.Chem import Draw
from rdkit import Chem
import matplotlib.pyplot as plt
class temp_class:
    def __init__(self,inf):
        self._inf = inf
    def rmol(self):
        try:
            mol = Chem.MolFromSmiles(self._inf)
        except:
            return 0
        return mol
    def rmol_img_mem(self):
        mol = self.rmol()
        try:
            img = Draw.MolToImage(mol,size=(800,800),kekulize=True,fitImage=True)
        except:
            return 0
        return img
    def rmol_img_file(self,fname,fPath):
        img = self.rmol_img_mem()
        # customizing figure
        plt.figure(figsize=(5,5))
        plt.imshow(img)
        plt.tight_layout()
        plt.axis("off")
        # notation
        plt.text(0,700,"* %s"%fname.strip(),fontsize=12.0,fontstyle="oblique")
        plt.text(0,750,"%s"%self._inf.strip(),fontsize=10.0)

        if fPath is None or fPath == "":
            fPath = "./"
        else:
            pass
        #plt.show()
        # save image
        plt.savefig(fPath + fname + ".png",dpi=100)

#############
# Code Test #
#############
if __name__ == "__main__":
    smi = "C3c5cc(Oc1nc2ccccc2(cc1))ccc5(OCC3Cc4cnccc4)"
    pp = temp_class(smi)
    img = pp.rmol_img_file("test","")