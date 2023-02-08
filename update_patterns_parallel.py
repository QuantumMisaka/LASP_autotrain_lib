# get coodination patterns from allstr.arc or TrainStr.txt
# can parallelly read exist pattern_db and update
# JamesBourbon in 20220807
# give the same result to single version, well used.
# simple parallel computation can be coding by multiprocessing package.

from multiprocessing.pool import MapResult
# from multiprocessing import Pool
import sys
import os

from allstr_new import AllStr
from structure_new import Str
from coordination_pattern import CoordinationPatterns
import time

NCORE=8

INTRO = '''
---------------------- Introduction ----------------------
Generate Coordination Patterns of AllStr and update database (if input)
Python3 Script for LASP Arc or TrainStr.txt Strucuture file
Parallel Version
Coding by James.Misaka.Bourbon.Liu, updated in 2022-0807

---------------------- Running Script ---------------------'''


def ParaWrap_Coor_Patterns(x):
    return x.coordination_pattern()

def get_patterns_parallel(strfile:str="TrainStr.txt",filemode:int=1, init_db="", ncore=NCORE):
    '''get all coordination patterns parallelly from structure
    
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
    # main for-loop running
    patterns_map_results = target_allstr.para_run(ParaWrap_Coor_Patterns, ncore)
    patterns_set = set()
    for result in patterns_map_results:
        result: MapResult
        for pattern_i in result.get():
            patterns_set.update(pattern_i)
    patterns_obj.update_from_patterns(patterns_set)
    # update and process
    print("---- DONE! print single_view and database versions of coordination patterns ----")
    patterns_simple_str = patterns_obj.print_all_coordinations_simple()
    patterns_json = patterns_obj.print_all_coordinations_full()
    patterns_obj.save_file(patterns_simple_str, "patterns_simple.txt")
    patterns_obj.save_file(patterns_json, "patterns_db.json")
    
    
    
if __name__ == "__main__":
    print(INTRO)
    start_time = time.perf_counter()
    if len(sys.argv) == 2:
        get_patterns_parallel(sys.argv[1])
    elif len(sys.argv) == 3:
        strfile, filemode = sys.argv[1:]
        get_patterns_parallel(strfile, eval(filemode))
    elif len(sys.argv) == 4:
        strfile, filemode, init_db = sys.argv[1:]
        get_patterns_parallel(strfile, eval(filemode), init_db)
    elif len(sys.argv) == 5:
        strfile, filemode, init_db, ncore = sys.argv[1:]
        get_patterns_parallel(strfile, eval(filemode), init_db, ncore)
    else:
        get_patterns_parallel()
    end_time = time.perf_counter()
    used_time = end_time - start_time
    print(f"---- TOTAL TIME USED: {used_time} ----")
