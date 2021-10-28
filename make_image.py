from rdkit.Chem import Draw
from rdkit import Chem
import matplotlib.pyplot as plt
import pybel
import pymol
import os
class molto2D:
    def __init__(self,inf):
        self._inf = inf
    def pmol(self):
        # Make Mol by using pybel
        mol = None
        # 1. Check SMILES #
        try:
            mol = list(pybel.readfile("smi",self._inf))[0]
        except:
            pass
        # 2. Check PDB #
        if mol is None:
            try:
                mol = pybel.readstring("smi",list(pybel.readfile("pdb",self._inf))[0].write("can").split()[0])
            except:
                pass
            if mol is None:
                # 3. Check SDF #
                try:
                    mol = pybel.readstring("smi",list(pybel.readfile("sdf",self._inf))[0].write("can").split()[0])
                except:
                    pass
                if mol is None:
                    try:
                        mol = pybel.readstring("smi",self._inf)
                    except:
                        pass
                    if mol is None:
                        return 0
                    else:
                        print("It is SMILES String")
                        return mol
                else:
                    print("It is SDF File")
                    return mol
            else:
                print("It is PDB File")
                return mol
        else:
            print("It is SMILES File")
            return mol
    def rmol(self):
        # Make Mol by RD-kit
        mol = None
        # 1. Check SMILES #
        try:
            mol = Chem.MolFromSmiles(self._inf)
        except:
            pass
        if mol is None:
            # 2. Check PDB #
            try:
                mol = Chem.MolFromSmiles(Chem.MolToSmiles(Chem.MolFromPDBFile(self._inf)))
            except:
                pass
            if mol is None:
                # 3. Check SDF #
                try:
                    mol = Chem.MolFromSmiles(Chem.MolToSmiles(next(Chem.SDMolSupplier(self._inf))))
                except:
                    pass
                if mol is None:
                    try:
                        with open(self._inf,"r") as F:
                            for line in F.readlines():
                                tline = line.strip()
                        mol = Chem.MolFromSmiles(tline)
                    except:
                        pass
                    if mol is None:
                        return 0
                    else:
                        print("It is SMILES File")
                        return mol
                else:
                    print("It is SDF File")
                    return mol
            else:
                print("It is PDB File")
                return mol
        else:
            print("It is SMILES String")
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

class molto3D:
    def __init__(self,inf):
        self._inf = inf
        pymol.cmd.bg_color("white")
        pymol.cmd.set("label_size",10)
        pymol.cmd.set("label_font_id",5)
        pymol.cmd.set("depth_cue",0)
    def load_chemical(self,fn,fPath):
        t = 1
        try:
            pymol.cmd.load(self._inf,"infile")
            t =1
        except:
            t = 0
        if t == 0:
            try:
                os.system("obabel -ismi %s -O %s/%s.sdf --gen3d"%(self._inf,fPath,fn))
            except:
                print("Error : Generate SMILES to 3D SDF")
                return t
            try:
                pymol.cmd.load(fPath + '/' + fn + ".sdf","infile")
                t = 1
            except:
                t = 0
        else:
            return t
        return t
    def stick_img(self,fname,fPath):
        re = self.load_chemical(fname,fPath)
        if re == 0:
            return
        else:
            pass
        pymol.cmd.show("stick")
        pymol.cmd.color("atomic")
        pymol.cmd.png(fPath + "/" + fname + ".png")
        pymol.cmd.save(fPath + "/" + fname + ".pse")



#############
# Code Test #
#############
if __name__ == "__main__":
    smi = "./example/example_chemical.txt" #"C3c5cc(Oc1nc2ccccc2(cc1))ccc5(OCC3Cc4cnccc4)"
    pp = molto3D(smi)
    pp.stick_img("test_sdf","./outputs")