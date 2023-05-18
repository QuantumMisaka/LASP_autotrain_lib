# for parallelly running collect allstr in auto.py
# last change by JamesBourbon in 20220903 due to gen_arc update

from allstr_new import AllStr as AllStr_new
import sys

def collect_allstr(workdirs,nbadneed):
    print("collect allstr from SSW nodejob.py")
    AllStr = AllStr_new() 
    AllStr.read_arc('%s/allstr.arc'%(workdirs), "") 
    str_count = len(AllStr)
    out_file = '%s/outstr.arc'%workdirs
    if (str_count == 0): 
        return
    elif str_count <= nbadneed:
        print('collect nbad not enough, but print outstr.arc to add')
        AllStr.gen_arc(range(str_count), out_file)
    else:
        print('collect Enough nbad, print outstr.arc add to VASP_DFT')
        AllStr.random_arange(200)
        AllStr.gen_arc(range(nbadneed), out_file)
    return

if __name__ == "__main__":
    workdir, nbadneed = sys.argv[1:]
    nbadneed = int(nbadneed)
    collect_allstr (workdir, nbadneed)
