# for parallelly running collect allstr in auto.py
# last change by JamesBourbon in 20220405

from allstr_new import AllStr as AllStr_new
import sys

def collect_allstr(workdirs,nbadneed):
    AllStr = AllStr_new() 
    AllStr.arcinit([0,0],'%s/allstr.arc'%(workdirs)) 
    str_count = len(AllStr)
    out_file = '%s/outstr.arc'%workdirs
    if (str_count == 0): 
        return
    elif str_count <= nbadneed:
        print('counting nbad: not enough')
        AllStr.gen_arc(list(range(str_count)), out_file,2)
    else:
        print('get enough nbad, gen-arc')
        AllStr.random_arange(200)
        AllStr.gen_arc(list(range(nbadneed)), out_file,2)
    return

if __name__ == "__main__":
    workdir, nbadneed = sys.argv[1:]
    nbadneed = int(nbadneed)
    collect_allstr (workdir, nbadneed)