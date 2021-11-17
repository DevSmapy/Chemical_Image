from rdkit.Chem import Draw
from rdkit import Chem
import matplotlib.pyplot as plt
import pybel
import pymol
import os
class molto2D: # several chemical format to 2D
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
            mol = Chem.MolFromSmiles(self._inf,sanitize=False)
        except:
            pass
        if mol is None:
            # 2. Check PDB #
            try:
                mol = Chem.MolFromSmiles(Chem.MolToSmiles(Chem.MolFromPDBFile(self._inf),sanitize=False),sanitize=False)
            except:
                pass
            if mol is None:
                # 3. Check SDF #
                try:
                    mol = Chem.MolFromSmiles(Chem.MolToSmiles(next(Chem.SDMolSupplier(self._inf)),sanitize=False),sanitize=False)
                except:
                    pass
                if mol is None:
                    try:
                        with open(self._inf,"r") as F:
                            for line in F.readlines():
                                tline = line.strip()
                        mol = Chem.MolFromSmiles(tline,sanitize=False)
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

##################
# Chemical to 3D #
##################
class molto3D: # several chemical format to 3D image
    def __init__(self,inf):
        self._inf = inf
        pymol.cmd.bg_color("white")
        pymol.cmd.set("label_size",10)
        pymol.cmd.set("label_font_id",5)
        pymol.cmd.set("depth_cue",0)
        pymol.cmd.color("atomic")
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
        pymol.cmd.hide("everything")
        pymol.cmd.show("stick")

        pymol.cmd.png(fPath + "/" + fname + ".stick.png")
        pymol.cmd.save(fPath + "/" + fname + ".stick.pse")
    def around_img(self,chain1,chain2,fname,fPath):
        re = self.load_chemical(fname,fPath)
        if re == 0:
            return
        else:
            pass
        pymol.cmd.hide("everything")
        pymol.cmd.select("chain %s within 5.0 of chain %s"%(chain1,chain2))
        pymol.cmd.show("surface","sele")
        pymol.cmd.set("transparency",0.3,"sele")
        pymol.cmd.show("stick","chain %s"%chain2)

        pymol.cmd.orient("visible") # Image location mediation

        pymol.cmd.png(fPath + "/" + fname + ".surface.png")
        pymol.cmd.save(fPath + "/" + fname + ".surface.pse")
    def All_cartoon_img(self,fname,fPath):
        re = self.load_chemical(fname,fPath)
        if re == 0:
            return
        else:
            pass
        pymol.cmd.show("cartoon")
        pymol.cmd.color("atomic")

        pymol.cmd.orient("visible")

        pymol.cmd.png(fPath + "/" + fname + ".A_cartoon.png")
        pymol.cmd.save(fPath + "/" + fname + ".A_cartoon.pse")
    def annotation(self,annot,fname,fPath):
        fn = fname.split(".")[0]
        if fPath is None or fPath == "":
            fPath = "./"
        else:
            pass
        try:
            img = plt.imread(fPath + "/" + fname)
        except:
            return 0
        plt.figure(figsize=(5,5))
        plt.imshow(img)
        plt.tight_layout()
        plt.axis("off")
        plt.text(0,500,"* %s"%self._inf.strip(),fontsize=12.0,fontstyle="oblique")
        plt.text(0,525,"%s"%annot,fontsize=10.0)

        plt.savefig(fPath + "/" + fn + ".annot.png")


#############
# Code Test #
#############
if __name__ == "__main__":
    smi = "./example/example_chemical.pdb" #"C3c5cc(Oc1nc2ccccc2(cc1))ccc5(OCC3Cc4cnccc4)"
    pp = molto3D(smi)
    pp.stick_img("test_stick","./outputs")
    #pp.around_img("A","X","test_surface","./outputs")
    pp.annotation("test annotation","test_stick.stick.png","./outputs")