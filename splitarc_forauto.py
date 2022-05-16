# for split arc-file to single-DFT based by lasp-autolin
# last change by JamesBourbon in 20220423

import sys
import os
import shutil
from allstr_new import AllStr as AllStr_new


def split(arcfile, paradir=''):
    workname="para"
    ROOTDIR = os.getcwd()
    AllStr = AllStr_new() 
    AllStr.arcinit([0,0],arcfile, forfile="")
    str_count = len(AllStr)
    for i in range(str_count):
        if paradir:
            workdir="%s/%s%d"%(ROOTDIR,workname,i)
            os.mkdir(workdir)
            incar='%s/INCAR'%(paradir)
            potcar='%s/POTCAR'%(paradir)
            inputfile='%s/lasp.in'%(paradir)
            shutil.copy(incar, "%s/INCAR"%workdir)
            shutil.copy(potcar, "%s/POTCAR"%workdir)
            shutil.copy(inputfile, "%s/lasp.in"%workdir)
            os.chdir(workdir)
            AllStr[i].outPOSCAR("POSCAR")
            AllStr[i].genKPOINTS("KPOINTS")
            AllStr.gen_arc([i], "lasp.str")
            os.chdir(ROOTDIR)
        else:
            AllStr.gen_arc([i], "outstr_%d.arc"%i)
    return


if __name__ == "__main__":
    arcfile, paradir = sys.argv[1:]
    split(arcfile, paradir)
