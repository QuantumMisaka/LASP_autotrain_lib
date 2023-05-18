# for split arc-file to single-DFT based by lasp_pylib
# last change by JamesBourbon in 20220815, 

import sys
import os
import shutil
from allstr_new import AllStr as AllStr_new


def split(arcfile, paradir='', workname="para"):
    ROOTDIR = os.getcwd()
    AllStr = AllStr_new() 
    AllStr.read_arc(arcfile)
    str_count = len(AllStr)
    for i in range(str_count):
        if paradir and os.path.isdir(paradir):
            workdir="%s/%s-%d"%(ROOTDIR,workname,i)
            if not os.path.isdir(workdir):
                os.mkdir(workdir)
            incar='%s/INCAR'%(paradir)
            inputfile='%s/lasp.in'%(paradir)
            shutil.copy(incar, "%s/INCAR"%workdir)
            # shutil.copy(potcar, "%s/POTCAR"%workdir)
            shutil.copy(inputfile, "%s/lasp.in"%workdir)
            os.chdir(workdir)
            # genPOTCAR need to specifiy abs-path or relative-path of paradir(setting)
            if os.path.isdir(paradir):
                AllStr[i].genPOTCAR(paradir, "POTCAR")
            else:
                print("abs path of paradir not find, try to find it in rootdir")
                abs_paradir = f"{ROOTDIR}/{paradir}"
                if os.path.isdir(abs_paradir):
                    AllStr[i].genPOTCAR(abs_paradir, "POTCAR")
                    print("find POTCAR!")
                else:
                    print("paradir NOT Find Error in genPOTCAR!")
            AllStr[i].outPOSCAR("POSCAR")
            AllStr[i].genKPOINTS("KPOINTS")
            AllStr.gen_arc([i], "lasp.str")
            os.chdir(ROOTDIR)
        else:
            # for other like nn
            workdir="%s/%s-%d"%(ROOTDIR,workname,i)
            os.mkdir(workdir)
            os.chdir(workdir)
            AllStr.gen_arc([i], "lasp.str")
            os.chdir(ROOTDIR)
    return


if __name__ == "__main__":
    if len(sys.argv) == 2:
        arcfile = sys.argv[1]
        split(arcfile)
    elif len(sys.argv) == 4:
        arcfile, paradir, workname = sys.argv[1:]
        split(arcfile, paradir, workname)
    else:
        arcfile, paradir = sys.argv[1:]
        split(arcfile, paradir)
