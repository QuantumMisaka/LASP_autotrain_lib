# get shrinked AlllStr from init AllStr/TrainStr
# for refine LASP Train-data
# JamesBourbon update in 20220827: soft/hard method can use together

import sys
import random as rd
import copy
from allstr_new import AllStr
from structure_new import Str
import numpy as np

STRFILE = "TrainStr_init.txt"
FILETYPE = 1
SOFT_NUM = 6
HARD_NUM = 6
RATIO_RANGE = [0.7, 0.9] # shrick ratio limit
# only 1st str will be used

def shrink_soft(init_str: Str, ratio_range=RATIO_RANGE):
    '''shrink init-structure by fractional coordination'''
    output_str = copy.deepcopy(init_str)
    # you can figure out why use deep-copy
    for index, vector in enumerate(output_str.abc[:3]):
        # output_str.set_coord()
        frac = output_str.FracCoord()
        ratio = rd.uniform(ratio_range[0], ratio_range[1])
        output_str.Latt[index] = vector * ratio
        # abc and Latt problem fixed
        output_str.Cell = output_str.Latt2Cell()
        # get all shrinked cart coordinate
        for i,atom in enumerate(output_str.atom):
            atom.xyz = np.dot(frac[i], output_str.Cell)
        output_str.set_coord()
    return output_str

def shrink_hard(init_str: Str, ratio_range=RATIO_RANGE):
    '''shrink init-structure by Cartesian coordination'''
    output_str = copy.deepcopy(init_str)
    for index, vector in enumerate(output_str.abc[:3]):
        # output_str.set_coord()
        ratio = rd.uniform(ratio_range[0], ratio_range[1])
        output_str.Latt[index] = vector * ratio
        # abc and Latt problem fixed
        output_str.Cell = output_str.Latt2Cell()
    return output_str
        
        
    
def main(strfile, soft_num = SOFT_NUM, hard_num = HARD_NUM ,file_type = FILETYPE):
    print("---- Running Shrinking Script ---- ")
    allstr_raw = AllStr()
    allstr_shrink_soft = AllStr()
    allstr_shrink_hard = AllStr()
    # find input-str type and read file
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
        allstr_raw.arcinit(strfile=strfile, forfile="") 
    elif file_type == 1:
        allstr_raw.train_data_init(strfile)
    # itered shrink it by soft/hard
    init_str = allstr_raw[0]
    for i in range(soft_num):
        soft_shrink_str = shrink_soft(init_str)
        allstr_shrink_soft.append(soft_shrink_str)
    for i in range(hard_num):
        hard_shrink_str = shrink_hard(init_str)
        allstr_shrink_hard.append(hard_shrink_str)
    # print as shrinked.arc
    output_file_soft = 'soft_shrinked.arc'
    output_file_hard = 'hard_shrinked.arc'
    allstr_shrink_soft.gen_arc(range(soft_num), output_file_soft)
    print(f"---- Generated {soft_num} shrinked sturcture in {output_file_soft} ----")
    allstr_shrink_hard.gen_arc(range(hard_num), output_file_hard)
    print(f"---- Generated {hard_num} shrinked sturcture in {output_file_hard} ----")
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
        print("---- Target File and soft_num is Specified ----")
        print("---- hard num will equal to soft num ----")
        main(strfile=sys.argv[1], soft_num=eval(
            sys.argv[2]), hard_num=eval(sys.argv[2]))
    elif argnum == 4:
        print("---- Target File and soft/hard_num is Specified ----")
        main(strfile=sys.argv[1], soft_num=eval(
            sys.argv[2]), hard_num=eval(sys.argv[3]))
    else:
        print("---- All Parameter is Specified! ----")
        main(sys.argv[1], eval(sys.argv[2]), 
                eval(sys.argv[3]), eval(sys.argv[4]))
    
        
    


