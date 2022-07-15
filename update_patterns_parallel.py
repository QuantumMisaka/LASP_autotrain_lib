# get coodination patterns from allstr.arc or TrainStr.txt
# can read exist pattern_db and update
# JamesBourbon in 20220714
# need parallel, not finished

import sys
import os

from allstr_new import AllStr, ParaWrap_CoorPatterns
# from structure_new import Str
from coordination_pattern import CoordinationPatterns
import time

NCORE=4


def get_patterns_parallel(strfile:str="TrainStr.txt",filemode:int=1, init_db=""):
    '''get all coordination patterns from structure
    
    Args:
        strfile: sturcture file name
        filemode: 0 for arc, 1 for TrainStr
    '''
    if not isinstance(filemode, int):
        print("Wrong FileMode! Program will Exit!")
        return 
    target_allstr = AllStr()
    if filemode == 0:
        target_allstr.arcinit(strfile=strfile)
    elif filemode == 1:
        target_allstr.train_data_init(strfile=strfile)
    else:
        print("Wrong FileMode! Program will Exit!")
        return
    print("---- Calculating coodination patterns of allstr/Trainstr ----")
    # parallel version   
    patterns_obj = CoordinationPatterns()
    # read from database
    if os.path.exists(init_db):
        print(f"---- Read {init_db} database for init-patterns ----")
        patterns_obj.read_coordination_json(init_db)
    patterns_set = target_allstr.para_run(ParaWrap_CoorPatterns, NCORE)
    patterns_obj.update_from_patterns(patterns_set)
    print("---- DONE! print single_view and database versions of coordination patterns ----")
    patterns_simple_str = patterns_obj.print_all_coordinations_simple()
    patterns_json = patterns_obj.print_all_coordinations_full()
    patterns_obj.save_file(patterns_simple_str, "patterns_simple.txt")
    patterns_obj.save_file(patterns_json, "patterns_db.json")
    
    
    
if __name__ == "__main__":
    print(INTRO)
    print("test mode")
    start_time = time.perf_counter()
    if len(sys.argv) == 2:
        get_patterns_parallel(sys.argv[1])
    elif len(sys.argv) == 3:
        strfile, filemode = sys.argv[1:]
        get_patterns_parallel(strfile, eval(filemode))
    elif len(sys.argv) == 4:
        strfile, filemode, init_db = sys.argv[1:]
        get_patterns_parallel(strfile, eval(filemode), init_db)
    else:
        get_patterns_parallel()
    end_time = time.perf_counter()
    used_time = end_time - start_time
    print(f"---- TOTAL TIME USED: {used_time} ----")
