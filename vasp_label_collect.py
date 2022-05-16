# collect and screen data from LASP-VASP-DFT to NN_train
# by JamesBourbon based by LASP_pylib, last change in 20220425
# sys.argv: 1 for file-set-name 2 for nmax to get

__author__ = "JamesBourbon"

from allstr_new import AllStr as AllStr_new
from allstr_new import BadStr
import sys
import os
import numpy as np
import glob


def collect_data(file_set_name="para", NELM_limit=120):
    '''collect VASP-DFT data
    
    Args:
        ncycle; cycle number
        nstr: structure number
        
    Returns:
        give 1_allstr.arc and 2_allfor.arc file to ROOTDIR
    '''
    strfile = "1_allstr.arc"
    forfile = "2_allfor.arc"
    os.system('rm -rf %s'%(strfile))
    os.system('rm -rf %s'%(forfile))
    workdirs=glob.glob("%s*"%(file_set_name))
    for workdir in workdirs:
        if glob.glob('%s/OSZICAR'%(workdir)):
            icontrol= int(os.popen('cat %s/OSZICAR | wc -l'%(workdir)).readline().strip())
            if(icontrol < NELM_limit): # if DFT-single end properly
                os.system('cat %s/allstr.arc >> %s'%(workdir, strfile))
                os.system('cat %s/allfor.arc >> %s'%(workdir, forfile))
    
    

def screen_data(strfile, forcefile=False, nmax=999999, maxF=500):
    '''for screen structure in strfile (but for what?) and force if exist
            key function for VASP-run str_file setting

        if finished, outstr.arc will be print-out (if force input, outfor.arc be print-out)

        noting by JamesBourbon in 20220402
    '''
    AllStr = AllStr_new()
    if forcefile:
        AllStr.arcinit([0, 0], strfile, forcefile)
    else:
        AllStr.arcinit([0, 0], strfile)
    # Here can set HighE,MaxAngle,MinAngle
    if len(AllStr) == 0:
        return
    # for screen structure/force data
    b = BadStr()
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

    print('All Str:', len(AllStr))
    #print 'present force',AllStr[0].Lfor
    if len(AllStr) == 0:
        print("no str input!")
        return
    if(len(AllStr) > nmax):
        AllStr.random_arange(200)
        AllStr = AllStr_new(AllStr[:(nmax)])

        #   AllStr.sort_by_energy()
    print('Final Dump Str:', len(AllStr))
    if len(AllStr) > 0:
        AllStr.gen_arc(list(range(len(AllStr))), 'outstr.arc', 2)
        if AllStr[0].Lfor:
            AllStr.gen_forarc(list(range(len(AllStr))), 'outfor.arc', 2)


def compare(vasp_allstr="1_allstr.arc", ssw_nn_allstr="allstr.arc-1"):
    '''compare NN_single and VASP-DFT result to get RMSE
            aimed to judge NN_train can finish or not 
        
        noted by JamesBourbon in 20220317
    '''
    # before VASP
    nn = AllStr_new()
    nn.readfile(ssw_nn_allstr)
    # After VASP
    vasp = AllStr_new()
    vasp.readfile(vasp_allstr)
    # compare energy
    if (len(nn)!= len(vasp)):
        print('some str failed to dft cal')
    else :
        err_all = 0
        for i in range(len(nn)):
            err_all += np.square(nn[i].Energy-vasp[i].Energy)
        rmse= np.sqrt(err_all/len(nn))
        print('----------- rmsE %14.8f eV ---------'%rmse)


def arc2train_data():    
    '''transfer outstr/outfor to TrainStr/TrainFor format file
        
    Returns:
        nadd: number of Structures in TrainStr.txt for NNtrain 
        put TrainStr.txt and TrainFor.txt in ROOTDIR
        
    noted by JamesBourbon in 20220403
    '''
    dir = os.path.abspath('.')
#   base = self.base
    tmpall= AllStr_new() 
    tmpall.readfile('%s/outstr.arc'%dir, '%s/outfor.arc'%dir) 
    # tmpall.shuffle(200) # random arange (not useful)
    # for str in tmpall: str.add_charge(autobase) # charge use default zero
    tmpall.gen_data_str(range(len(tmpall)), 'TrainStr.txt')
    tmpall.gen_data_for(range(len(tmpall)), 'TrainFor.txt')
    return len(tmpall) # num of str add to train



if __name__ == "__main__":
    ROOTDIR=os.getcwd()
    NN_allstr="allstr.arc-1"
    
    if len(sys.argv) == 2:
        collect_data(file_set_name=sys.argv[1])
        screen_data("1_allstr.arc", "2_allfor.arc")
    elif len(sys.argv) == 3:
        collect_data(file_set_name=sys.argv[1])
        screen_data("1_allstr.arc", "2_allfor.arc", sys.argv[2])
    else:  
        collect_data()
        screen_data("1_allstr.arc", "2_allfor.arc")
        
    if glob.glob("%s/%s"%(ROOTDIR,NN_allstr)):
        compare()
    nadd = arc2train_data()
    
