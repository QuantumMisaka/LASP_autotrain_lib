# screen from SSW-NN to vasp-DFT
# JamesBourbon based on LASP-pylib in 20220502

import sys
from allstr_new import AllStr as AllStr_new
from allstr_new import BadStr

def screen_data(strfile, forcefile, nmax = 999999, maxF=10):
    '''for screen structure in strfile and force if exist
        key function for VASP-run str_file setting
        
    if finished, outstr.arc will be print-out (if force input, outfor.arc be print-out)
        
    noting by JamesBourbon in 20220402
    
    '''        
    AllStr = AllStr_new()
    if forcefile:
        AllStr.arcinit([0,0],strfile,forcefile)
    else :
        AllStr.arcinit([0,0],strfile)
    # Here can set HighE,MaxAngle,MinAngle
    if len(AllStr)==0: return
    # for screen structure/force data
    b=BadStr()
    #   b.HighE=-3.0
    b.MaxFor = maxF
    b.MaxLat = 40
    b.MinLat = 2.2
    #
    AllStr = AllStr.filter(b)

    # main filter function use allstr_new.py and structure_new.py
    if(len(AllStr) > nmax+50):
        AllStr.random_arange(200)
        AllStr = AllStr_new(AllStr[:(nmax+50)])

    print('All Str:',len(AllStr))
    #print 'present force',AllStr[0].Lfor

    if len(AllStr)==0: return
    if(len(AllStr) > nmax):
        AllStr.random_arange(200)
        AllStr = AllStr_new(AllStr[:(nmax)])

    #   AllStr.sort_by_energy()

    print('Final Str Num:',len(AllStr))
    if len(AllStr) >0:
        AllStr.gen_arc(list(range(len(AllStr))),'outstr.arc',2)
        if AllStr[0].Lfor: 
            AllStr.gen_forarc(list(range(len(AllStr))),'outfor.arc',2)
            
if __name__ == "__main__":
    if len(sys.argv) == 3:
        strfile, forfile = sys.argv[1:]
        screen_data(strfile, forfile)
    if len(sys.argv) == 4:
        strfile, forfile, nmax = sys.argv[1:]
        screen_data(strfile, forfile, nmax)
    if len(sys.argv) == 5:
        strfile, forfile, nmax, maxF = sys.argv[1:]
        screen_data(strfile, forfile, nmax, maxF)