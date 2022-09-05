# for parallelly running collect allstr in auto.py
# can directly used to choose traindata_update structure
# last change by JamesBourbon in 20220903 due to gen_arc update

from allstr_new import AllStr as AllStr_new
# from structure_new import Str
from coordination_pattern import CoordinationPatterns
import sys
import os

# init setting
new_pattern_limit = 10
patterns_db = "../traindata_patterns.json" # db in workdir, run in SSWdir
chosen_patterns = CoordinationPatterns(name="train_data", limit=new_pattern_limit)
if os.path.exists(patterns_db):
    print("---- Read Exist Coor-Patterns Database ----")
    chosen_patterns.read_coordination_json(patterns_db)

def collect_allstr(workdir,nbadneed):
    choosing_volume = 20*nbadneed
    # volume give a larger choosing space
    print("collect allstr from SSW nodejob_coor.py")
    AllStr = AllStr_new()
    AllStr.arcinit([0,0],'%s/allstr.arc'%(workdir)) 
    AllStrGot = AllStr_new() 
    str_count = len(AllStr)
    out_file = '%s/outstr.arc'%workdir
    if (str_count == 0): 
        return
    elif str_count <= choosing_volume:
        print('counting nbad: not enough')
        # AllStr.gen_arc(list(range(str_count)), out_file)
    else:
        # use coor_pattern method to collect structure
        print('collect Enough nbad, print outstr.arc add to VASP_DFT')
        AllStr.random_arange(10)
        count = 0
        for struc in AllStr:
            # struc:Str
            one_str_patterns = struc.coordination_pattern()
            updating = chosen_patterns.update_patterns_from_structure(
                one_str_patterns)
            count += 1
            if updating:
                AllStrGot.append(struc)
                if len(AllStrGot) >= nbadneed:
                    print(f"---- effictive/total = {nbadneed}/{count} ----")
                    break
        chosen_patterns.print_operation_log() # print coor patterns log
        AllStrGot.gen_arc(list(range(nbadneed)), out_file)
        #  DO NOT update json coor-pattern database now: screen_data is needed backward
        # json_string = chosen_patterns.print_all_coordinations_full()
        # with open(patterns_db, 'w') as fo:
        #     fo.write(json_string)
        
    return

if __name__ == "__main__":
    workdir, nbadneed = sys.argv[1:]
    nbadneed = int(nbadneed)
    collect_allstr (workdir, nbadneed)
