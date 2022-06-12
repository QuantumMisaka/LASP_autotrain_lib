# dynamic verlet coordination pattern choose structure
# JamesBourbon in 20220611

from fileinput import filename
from coordination_pattern import CoordinationPatterns
from allstr_new import AllStr as AllStr_new
from structure_new import Str
import os

project_name = "Au8Pd8O24_2010_new"
filename = "test/allstr.arc"
ACCEPT_COUNT = 20
REJECT_COUNT = 10
input_allstr = AllStr_new()
input_allstr.arcinit(strfile=filename, forfile="")
target_allstr = AllStr_new()
patterns_db = CoordinationPatterns(
        name=project_name, output=project_name+"_coor_run.log" )



def main_run():
    steps = 0
    stop = False
    steps += ACCEPT_COUNT
    update_count = 0
    dir_list = []
    # filter
    for struc in input_allstr:
        struc: Str
        if struc.Energy < -6000:
            input_allstr.remove(struc)
    step_limit = len(input_allstr)
    while not stop:
        # running sampling
        target_name = project_name+"_"+str(steps+1)
        target_str : Str = input_allstr[steps]
        target_patterns = target_str.coordination_pattern()
        updated = patterns_db.update_patterns_from_structure(
            target_patterns, target_name)
        if updated:
            update_count += 1
            steps += ACCEPT_COUNT
            target_allstr.append(target_str)
        else:
            steps += REJECT_COUNT
        if steps > step_limit:
            stop = True
    # stop here: just get patterns
    patterns_db.save_file(
        patterns_db.print_all_coordinations_simple(),
        patterns_db.name+"_simple.txt"
    )
    patterns_db.save_file(
        patterns_db.print_all_coordinations_full(),
        patterns_db.name+"_full.json"
    )
    patterns_db.print_operation_log()
  #   for struc in target_allstr:
    #     struc: Str
      #   str.outPOSCAR()
    
if __name__ == "__main__":
    main_run()
    
    