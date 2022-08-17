# get shrinked AlllStr from init AllSTr/TrainStr
# for refine LASP Train-data
# JamesBourbon update in 20220815

import sys
import os
import random as rd
import copy
from allstr_new import AllStr
from structure_new import Str
import numpy as np

STRFILE = "TrainStr_init.txt"
FILETYPE = 1
SHRINK_NUM = 20
RATIO_RANGE = [0.68, 0.86] # shrick ratio limit
# only 1st str will be used

def shrink(init_str: Str, ratio_range=RATIO_RANGE):
    output_str = copy.deepcopy(init_str)
    # you can figure out why use deep-copy
    for index, vector in enumerate(output_str.abc[:3]):
        # output_str.set_coord()
        frac = output_str.FracCoord()
        ratio = rd.uniform(ratio_range[0], ratio_range[1])
        output_str.abc[index] = vector * ratio
        output_str.Latt = output_str.abc 
        # abc and Latt problem should find a time to sort out -- tips
        output_str.Cell = output_str.Latt2Cell()
        # get all shrinked cart coordinate
        for i,atom in enumerate(output_str.atom):
            atom.xyz = np.dot(frac[i], output_str.Cell)
        output_str.set_coord()
    return output_str
        
    
def main(strfile, shrink_num = SHRINK_NUM, file_type = FILETYPE):
    print("---- Running Shrinking Script ---- ")
    allstr_raw = AllStr()
    allstr_shrink = AllStr()
    if strfile[-4:] == '.txt':
        print("input file is TrainStr file, file_type = 1") 
        file_type = 1
    elif strfile[-4:] == '.arc': 
        print("input file is ARC file, file_type = 0")
        file_type = 0
    else: 
        print("input file type not detect")
        print(f"use default setting file_type = {file_type}")
    if file_type == 0:
        allstr_raw.arcinit(strfile=strfile) 
    elif file_type == 1:
        allstr_raw.train_data_init(strfile)
    
    init_str = allstr_raw[0]
    for i in range(shrink_num):
        shrink_str = shrink(init_str)
        allstr_shrink.append(shrink_str)
    # print as shrink_allstr.arc
    output_file = 'shrink_allstr.arc'
    allstr_shrink.gen_arc(range(shrink_num), output_file)
    print(f"---- Generated {shrink_num} shrinked sturcture in {output_file} ----")
    print("DONE!")
    
if __name__ == "__main__":
    argnum = len(sys.argv)
    if argnum == 1:
        print("---- Use Default Setting ----")
        main()
    elif argnum == 2:
        print("---- Target File is Specified ----")
        main(strfile=sys.argv[1])
    elif argnum == 3:
        print("---- Target File and Shrink_num is Specified ----")
        main(strfile=sys.argv[1], shrink_num=eval(sys.argv[2]))
    else:
        print("---- All Parameter is Specified! ----")
        main(sys.argv[1], eval(sys.argv[2]), eval(sys.argv[3]))
    
        
    


