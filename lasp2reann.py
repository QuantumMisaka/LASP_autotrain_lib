# lasp arc/Train file to REANN configuration file for training
# JamesMisaka in 20230503

import sys
from allstr_new import AllStr


def lasp2reann(strname, forname, filetype=1):
    '''
    lasp allstr.arc/allfor.arc or TrainStr.txt/TrainFor.txt 
    to REANN

    Inputs:
        strname: input str file (.arc or .txt will be identified)
        forname: input for file
        filetype: 0 for allstr/for.arc 1 for TrainStr/For.txt
    '''
    allstr = AllStr()
    if filetype==0 or strname[-4:] == '.arc':
        allstr.arcinit([0,0], strname, forname)
    elif filetype==1 or strname[-4:] == '.txt':
        allstr.train_data_init(strname, forname)
    else:
        print("filetype should be correct !")

    allstr.gen_REANN_configuration(range(len(allstr)))

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Input Format: lasp2reann.py <strfile> <forfile>")
    elif len(sys.argv) == 3:
        strname, forname = sys.argv[1], sys.argv[2]
        lasp2reann(strname, forname)
    else:
        strname, forname, filetype = sys.argv[1:]
        lasp2reann(strname, forname, filetype)





