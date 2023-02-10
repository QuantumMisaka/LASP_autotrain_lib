# dynamic verlet coordination pattern choose structure from arcfile
# JamesBourbon in 20220717

from coordination_pattern import CoordinationPatterns
from allstr_new import AllStr as AllStr_new
from structure_new import Str
import os
import time
import sys
# import shutil

# parameter and const
PROJECT_NAME = "Au8Pd8O24_5020"
FILENAME = "test/allstr.arc"
COOR_ACCEPT_RATE = 1
# rate(limit): >1 for new patterns num, <1 for new patterns ratio
ACCEPT_COUNT = 50
REJECT_COUNT = 20
ENERGY_LIMIT = -6000 # can disable filter by =-100000
# for generate dft-label project
# DFT_SET = ""
DFT_SET = "dft_setting" # setting dir name, empty for not use
OUTPUT_DIR = "sampled_out"
PROJECT_DIR = "dft_jobs"
COOR_DB_READ = "patterns_input.json" # read db if have
FATHER_DIR = os.getcwd()
# variable
SCR_DIR = "scratch"
input_allstr = AllStr_new()
input_allstr.arcinit(strfile=FILENAME, forfile="")
target_allstr = AllStr_new()
patterns_db = CoordinationPatterns(name=PROJECT_NAME, 
                output=PROJECT_NAME+"_coor_run.log", limit=COOR_ACCEPT_RATE)
if bool(COOR_DB_READ) and os.path.exists(COOR_DB_READ):
    patterns_db.read_coordination_json(COOR_DB_READ)


INTRO = '''
---------------------- Introduction ----------------------
Verlet-Choose Strutures Dataset by Coordination Patterns in Structure
Python3 Script for LASP SSW-sampled allstr.arc or all.arc
Coding by James.Misaka.Bourbon.Liu, updated in 2022-0716
Main Function is to choose unique structrue to add to database

---------------------- Running Script ---------------------'''


def main_run():
    print(INTRO)
    steps = 0
    stop = False
    steps += ACCEPT_COUNT
    update_count = 0
    target_dir_list = []
    # filter, can be muted
    if ENERGY_LIMIT > -10000:
        print("Structures having lower than" )
        print(f"{ENERGY_LIMIT} eV Energy will be filtered")
        for struc in input_allstr:
            struc: Str
            if struc.energy < ENERGY_LIMIT:
                input_allstr.remove(struc)
    step_limit = len(input_allstr)
    while not stop:
        # running sampling
        target_name = PROJECT_NAME+"_"+str(steps+1)
        target_str : Str = input_allstr[steps]
        target_patterns = target_str.coordination_pattern()
        updated = patterns_db.update_patterns_from_structure(
            target_patterns, target_name)
        if updated:
            update_count += 1
            steps += ACCEPT_COUNT
            target_dir_list.append(target_name)
            target_allstr.append(target_str)
        else:
            steps += REJECT_COUNT
        if steps >= step_limit:
            stop = True
    # print-out
    output_dir = FATHER_DIR+"/"+OUTPUT_DIR
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    os.chdir(output_dir)
    patterns_db.save_file(
        patterns_db.print_all_coordinations_simple(),
        patterns_db.name+"_simple.txt"
    )
    patterns_db.save_file(
        patterns_db.print_all_coordinations_full(),
        patterns_db.name+"_full.json"
    )
    patterns_db.print_operation_log()
    # print out
    print(f'---- choose {update_count} structure in {step_limit} allstr! ----')
    dft_setting = f"{FATHER_DIR}/{DFT_SET}"
    if os.path.isdir(dft_setting) and bool(DFT_SET):
        # print-out directly used for single-DFT
        print("--- DFT SETTING dir exists! ---")
        print(f"--- output sampled result to {OUTPUT_DIR}/{PROJECT_DIR} ---")
        if os.path.isdir(PROJECT_DIR):
            print("PROJECT_DIR exist!")
            if not os.path.isdir(SCR_DIR):
                os.mkdir(SCR_DIR)
            if bool(os.listdir(PROJECT_DIR)):
                print("PROJECT NOT EMPTY!")
                temp_dir = time.strftime("%d-%H-%M-%S", time.gmtime())
                print(f"All subjects will backup to {SCR_DIR}/{temp_dir}")
                os.mkdir(f"{SCR_DIR}/{temp_dir}")
                os.system(f"mv {PROJECT_DIR}/* {SCR_DIR}/{temp_dir}/")
        else:    
            os.mkdir(PROJECT_DIR)
        pro_dir = f"{FATHER_DIR}/{OUTPUT_DIR}/{PROJECT_DIR}"
        os.chdir(pro_dir)
        # get dft-jobs
        for index, dirname in enumerate(target_dir_list):
            os.mkdir(dirname)
            os.system(f'cp {dft_setting}/lasp.in {dft_setting}/INCAR {dirname}')
            os.chdir(dirname)
            struc = target_allstr[index]
            # if structure is cluster, KPOINTS need refine to gamma
            struc.genKPOINTS("KPOINTS") 
            struc.outPOSCAR("POSCAR")
            # update 20220720: print lasp.str from Str._Gen_arc
            struc._Gen_arc(coord=struc.Coord, fname='lasp.str')
            struc.genPOTCAR(sourcedir=dft_setting, outfile="POTCAR")
            os.chdir(pro_dir)
    else:
        # print as outstr.arc
        print("--- DFT SETTING dir NOT exists! ---")
        target_allstr.gen_arc(range(update_count))
        print("--- sampled result print as outstr.arc ! ---")
    print("---- DONE!!! ----")
    
if __name__ == "__main__":
    main_run()
    
    