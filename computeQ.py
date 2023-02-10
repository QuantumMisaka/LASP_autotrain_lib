
__author__ = 'JamesBourbon learn from zpliu-LASP_pythonlib'

# modifed from screenminimum -- only for compute Q
import time
import cmath as m
import os
import numpy as np
import random
import sys
#import template as Tp
import string
# sys.path.append('/home11/Liugroup/tools/')
from multiprocessing import Pool

from structure_new import Str, ParaWrap_SteinhartQ_cal ,ParaWrap_JudgeShape 
from allstr_new import BadStr 
from allstr_new import AllStr as allstr_new


if __name__ == "__main__":
    allstr = allstr_new()
    allstr.arcinit([0,0],'all.arc')    # structure file name of all.arc
    print('All Str:',len(allstr))
    Ncore = 56

    print(' python job started at', time.strftime('%d/%m/%Y %H:%M:%S'))

# parallel version
    #Ncore=4
#   stage 1
  
    f= ParaWrap_SteinhartQ_cal
    results = allstr.para_run(f,Ncore)
    record=0
    b=[]
    for x in results:
        for y1,y2,y3,y4,y5,y6 in x.get():
            record +=1
            allstr[record-1].Q.append(y1)
            allstr[record-1].Q.append(y2)
            allstr[record-1].Q.append(y3)
            allstr[record-1].Q.append(y4)
            allstr[record-1].Q.append(y5)
            allstr[record-1].Q.append(y6)


#   if len(allstr) >0:
#       allstr.Gen_arc(range(len(allstr)))


    print("Summary of Stein-Q value Minimum")
    print("                     Q2          Q4          Q6   ")
    for i in range(len(allstr)):
        print("%5d   %8.4f  %8.4f  %8.4f  %8.4f "%(i,allstr[i].energy,float(allstr[i].Q[0]),float(allstr[i].Q[1]),float(allstr[i].Q[2])))

    print('Finish ---------------------------------------------')


    print("Summary of Distance-weighted Stein-Q value Minimum")
    print("                     Q2          Q4          Q6   ")
    for i in range(len(allstr)):
        print("%5d   %8.4f  %8.4f  %8.4f  %8.4f  "%(i,allstr[i].energy,float(allstr[i].Q[3]),float(allstr[i].Q[4]),float(allstr[i].Q[5])))

    print('Finish ---------------------------------------------')
