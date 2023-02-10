# get coodination patterns from allstr.arc or TrainStr.txt
# can read exist pattern_db and update
# JamesBourbon in 20220714
# for large system, losts of time will consume
# need parallel test

import sys
import os

from allstr_new import AllStr, ParaWrap_CoorPatterns
# from structure_new import Str
from coordination_pattern import CoordinationPatterns
import time

INTRO = '''
---------------------- Introduction ----------------------
Generate Coordination Patterns of AllStr and update database (if input)
Python3 Script for LASP Arc or TrainStr.txt Strucuture file
Coding by James.Misaka.Bourbon.Liu, updated in 2022-0716

---------------------- Running Script ---------------------'''

def get_patterns(strfile:str="TrainStr.txt",filemode:int=1, init_db=""):
    '''get all coordination patterns from structure, 
    and update coordination patterns
    
    Args:
        strfile: sturcture file name
        filemode: 0 for arc, 1 for TrainStr
        init_db: target update database, default is None
        
    Returns:
        give patterns_db.json database in execute_dir
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
    # single node
    patterns_obj = CoordinationPatterns()
    # read from database
    if os.path.exists(init_db):
        print(f"---- Read {init_db} database for init-patterns ----")
        patterns_obj.read_coordination_json(init_db)
    patterns_set = target_allstr.gen_coordination_patterns()
    patterns_obj.update_from_patterns(patterns_set)
    print("---- DONE! print single_view and database versions of coordination patterns ----")
    patterns_simple_str = patterns_obj.print_all_coordinations_simple()
    patterns_json = patterns_obj.print_all_coordinations_full()
    patterns_obj.save_file(patterns_simple_str, "patterns_simple.txt")
    patterns_obj.save_file(patterns_json, "patterns_db.json")
    
    
if __name__ == "__main__":
    print(INTRO)
    start_time = time.perf_counter()
    if len(sys.argv) == 2:
        get_patterns(sys.argv[1])
    elif len(sys.argv) == 3:
        strfile, filemode = sys.argv[1:]
        get_patterns(strfile, eval(filemode))
    elif len(sys.argv) == 4:
        strfile, filemode, init_db = sys.argv[1:]
        get_patterns(strfile, eval(filemode), init_db)
    else:
        get_patterns()
    end_time = time.perf_counter()
    used_time = end_time - start_time
    print(f"---- TOTAL TIME USED: {used_time} ----")
