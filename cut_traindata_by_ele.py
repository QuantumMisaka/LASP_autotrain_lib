# cut traindata by specific element conbination
# JamesBourbon in 20220730

import sys
import os
from allstr_new import AllStr
from structure_new import Str

ELEMENT_CONB = "Si O Sn"
STRFILE="TrainStr.txt"
FORFILE="TrainFor.txt" # can be none

def main(element_conb: set, strfile: str, forfile:str=""):
    # read TrainStr and TrainFor
    traindata_raw = AllStr()
    if not os.path.exists(strfile):
        print("TrainStr NOT exist and program will exit !!")
        return
    if os.path.exists(forfile):
        traindata_raw.train_data_init(strfile, forfile)
    else:
        traindata_raw.train_data_init(strfile)
    raw_count = len(traindata_raw)
    # cut to other allstr
    traindata_cut = AllStr()
    for struc in traindata_raw:
        struc: Str
        if set(struc.ele_nameList) == element_conb:
            traindata_cut.append(struc)
    # print out
    element_conb_sort = sorted(list(element_conb))
    ele_conb_tick = ''.join(element_conb_sort)
    str_out = f"{ele_conb_tick}_{strfile}"
    for_out = f"{ele_conb_tick}_{forfile}"
    cut_count = len(traindata_cut)
    printlist = range(cut_count)
    traindata_cut.gen_data_str(printlist, str_out)
    if forfile:
        traindata_cut.gen_data_for(printlist, for_out)
    print(f"---- Cut {cut_count} Structure from {raw_count} Train-data ----")
    print(f"---- OutPut StrFile as {str_out} ----")
    print("DONE!")

if __name__ == "__main__":
    
    HELPTXT = '''
    Cut Traindata by Element Conbination like 'Fe Si O'
    Written by James.Misaka.Bourbon.Liu in 20220730
    Usage: 
    1. Edit ELEMENT_CONB, STRFILE, FORFILE in python script
    2. Run by "python [exec] [element1] [element2] ... to editor element conbination
    Notice:
    1. Confirm your TrainStr.txt and TrainFor.txt Name
    2. Confirm your TrainStt.txt and TrainFor.txt from same Train-data
    3. Element conbination should be separated by SPACE
    4. Only cut the Train-data which have ALL the element in provided conbination
    '''
    
    argnum = len(sys.argv)
    if argnum == 2 and '-h' in sys.argv[1]:
        print(HELPTXT)
        exit()
    elif argnum == 1:
        print("---- Use Argments Written to Run Cutting ----")
        ele_comb_set = set(ELEMENT_CONB.split())
    elif argnum >= 1:
        print("---- Notice: Only Elements Conbination can be specified ----")
        ele_comb_set = set(sys.argv[1:])
    main(ele_comb_set, STRFILE, FORFILE)

        