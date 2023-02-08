# analysis traindata
# JamesBourbon update in 20220725
# update function for analysis elements conbination of all struc
import numpy as np
import sys
from allstr_new import AllStr
# from structure_new import Str


# setting
strfile = "TrainStr.txt"
forfile = "TrainFor.txt"

# reading TrainFile
allstr_raw = AllStr()
if len(sys.argv) == 2:
    allstr_raw.train_data_init(sys.argv[1])
elif len(sys.argv) > 2:
    allstr_raw.train_data_init(sys.argv[1], sys.argv[2])
else:
    allstr_raw.train_data_init(strfile, forfile)    
# out of range： TrainStr not related to TrainFor

# get Traindata info
size = len(allstr_raw)
check_25 = round(size/4)-1
check_50 = round(size/2)-1
check_75 = round(size * 3/4)-1
check_90 = round(size * 0.9)
max_for = 0
max_stress = 0
all_elements = allstr_raw.get_all_element()
natom_stat_by_ele = {}
# do statistic by sort
allstr_sort_by_natom = allstr_raw.sort_by_natom()
natom_stat = [allstr_sort_by_natom[0].natom, allstr_sort_by_natom[check_25].natom,
            allstr_sort_by_natom[check_50].natom, allstr_sort_by_natom[check_75].natom, 
            allstr_sort_by_natom[check_90].natom, allstr_sort_by_natom[size-1].natom]
# natom of each element
for ele_name in all_elements:
    allstr_sort_by_ele_natom = allstr_raw.sort_by_element_natom(ele_name=ele_name)
    ele_natom_stat = [allstr_sort_by_ele_natom[0].sp.get(ele_name, 0),
                    allstr_sort_by_ele_natom[check_25].sp.get(ele_name, 0),
                    allstr_sort_by_ele_natom[check_50].sp.get(ele_name, 0),
                    allstr_sort_by_ele_natom[check_75].sp.get(ele_name, 0), 
                    allstr_sort_by_ele_natom[check_90].sp.get(ele_name, 0),
                    allstr_sort_by_ele_natom[size-1].sp.get(ele_name, 0)]
    natom_stat_by_ele[ele_name] = ele_natom_stat
# max force
allstr_sort_by_force = allstr_raw.sort_by_force()
max_for = allstr_sort_by_force[size-1].maxF
# max stress
allstr_sort_by_stress = allstr_raw.sort_by_stress()
max_stress = max(allstr_sort_by_stress[size-1].stress)

# Traindata element conbination
elements_conbination_dict = {}
for stru in allstr_raw:
    # stru : Str
    ele_conb = tuple(stru.sp.keys())
    elements_conbination_dict[ele_conb] = elements_conbination_dict.get(ele_conb, 0) + 1

# print Traindata info
info_string = "---- Traindata Analysis Result ----\n"
info_string += f" Data Size: {size} (structures)\n"
info_string += f" Max Force: {max_for} (eV/Ang)\n"
info_string += f" Max Stress {max_stress} (GPa)\n"
info_string += "Table for Number of Atoms in Struc \n"
info_string += "Ele\t 0\t 25%\t 50%\t 75%\t 90%\t 100%\n"
info_string += "All\t %d\t %d\t %d\t %d\t %d\t %d\n"%(
    natom_stat[0], natom_stat[1], natom_stat[2], natom_stat[3], natom_stat[4], natom_stat[5])
for ele in all_elements:    
    info_string += "%s\t %d\t %d\t %d\t %d\t %d\t %d\n"%(
        ele, natom_stat_by_ele[ele][0], natom_stat_by_ele[ele][1], natom_stat_by_ele[ele][2], 
        natom_stat_by_ele[ele][3], natom_stat_by_ele[ele][4], natom_stat_by_ele[ele][5])
info_string += "Elements Conbinations in Train-data: \n"
for ele_conb, count in elements_conbination_dict.items():
    info_string += f"%s: %d \n"%(ele_conb, count)

info_string += "---- DONE! ----"

print(info_string)
