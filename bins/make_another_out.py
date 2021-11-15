import matplotlib.pyplot as plt
import glob,os

class tiled:
    def __init__(self,fPath):
        try:
            self.f_list = sorted(glob.glob(fPath + "/*.png"))
        except:
            self.f_list = sorted(glob.glob(fPath + "/*.jpg"))
    def divide_list(self,n):
        li = self.f_list
        if n == 0: num = 1
        else: num = n

        for i in range(0,len(li),num): yield li[i:i+num]
    def make_tild_pic(self,a,b,fname,fPath):
        tf_list = list(self.divide_list(int(a)*int(b)))
        alp = 0
        for t in tf_list:
            alp += 1
            rows = int(a)
            cols = int(b)

            fig = plt.figure(figsize=(rows * 6, cols * 6), constrained_layout=True)
            i = 1
            for f in t:
                img = plt.imread(f)
                ax = fig.add_subplot(rows,cols,i)
                ax.imshow(img)
                ax.set_xticks([]),ax.set_yticks([])
                ax.axis("off")
                i += 1

            plt.tight_layout()
            plt.subplots_adjust(wspace=-.2,hspace=-.5)
            plt.savefig(fPath + "/" + str(alp).rjust(3,"0") + "_" + fname + ".tiled.png")
            plt.close()
if __name__== "__main__":
    p = tiled("../example/example_tiled/")
    p.make_tild_pic(5,5,"test","../outputs")