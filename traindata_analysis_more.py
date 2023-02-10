# analysis traindata
# JamesMisaka update in 20230208
# more infomation get from train-data
# refined and well used

from allstr_new import AllStr
from structure_new import Str
import numpy as np
import pandas as pd
import sys
import time



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
max_for = 0
max_stress = 0
all_elements = allstr_raw.get_all_element()
# do statistic by sort
# max force
allstr_sort_by_force = allstr_raw.sort_by_force()
max_for = allstr_sort_by_force[size-1].maxF
# max stress
allstr_sort_by_stress = allstr_raw.sort_by_stress()
max_stress = max(allstr_sort_by_stress[size-1].stress)



# forset all-in-one
# Specific infomation about Traindata element conbination and cell type
atom_conb_dict = {}
# get cell type: bulk, layer, cluster
elements_conb_dict = {}
# atom_conb_dict = {}
natom_dict = {}
cell_type_dict = {}
cell_and_atom_type_dict = {}


# can be parallel computation
start_time = time.perf_counter()
for stru in allstr_raw:
    stru : Str
    ele_conb = tuple(stru.sp.keys())
    elements_conb_dict[ele_conb] = elements_conb_dict.get(ele_conb, 0) + 1
    atom_conb = ""
    natom = stru.natom
    for atom, num in stru.sp.items():
        atom_conb += atom
        atom_conb += str(num)
    natom_dict[atom_conb] = natom_dict.get(atom_conb, 0) + 1
    cell_type = stru.get_basic_shape() # key time comsuming step
    cell_and_atom = (atom_conb, cell_type)
    atom_conb_dict[atom_conb] = atom_conb_dict.get(atom_conb, 0) + 1
    cell_type_dict[cell_type] = cell_type_dict.get(cell_type, 0) + 1
    cell_and_atom_type_dict[cell_and_atom] = cell_and_atom_type_dict.get(cell_and_atom, 0) + 1
end_time = time.perf_counter()
time_consumed = end_time - start_time

# print Traindata info
print(f" ---- Processing Done in {time_consumed} s")
print(" ---- Getting Traindata infomation ...")
info_string = "---- Traindata Analysis Result ----\n"
info_string += f" Data Size: {size} (structures)\n"
info_string += f" Max Force: {max_for} (eV/Ang)\n"
info_string += f" Max Stress {max_stress} (GPa)\n"
info_string += "Elements Conbinations in Train-data: \n"
for ele_conb, count in elements_conb_dict.items():
    info_string += "%s: %d \n"%(ele_conb, count)
info_string += "Atoms Conbinations in Train-data: \n"
for atom_conb, count in atom_conb_dict.items():
    info_string += "%s: %d \n"%(atom_conb, count)
info_string += "Structure Types in Train-data: \n"
for cell_type, count in cell_type_dict.items():
    info_string += f"{cell_type}: {count}\n"
info_string += "Atoms Conbinations and their Types in Train-data: \n"
for cell_and_atoms, count in cell_and_atom_type_dict.items():
    info_string += f"{cell_and_atoms}: {count}\n"
# for atom_conb, count in atom_conb_dict.items():
#     info_string += f"{atom_conb}: {count}\n"

info_string += "---- DONE! ----"

print(" ---- Printing Traindata infomation table ...")
# give a csv file for Traindata
atoms_and_type_info = {}
for atom_conb, count in atom_conb_dict.items():
    atoms_and_type_info[atom_conb] = {
        "natom": natom_dict.get(atom_conb, 0),
        "bulk": cell_and_atom_type_dict.get((atom_conb, "bulk"),0),
        "layer": cell_and_atom_type_dict.get((atom_conb, "layer"),0),
        "cluster": cell_and_atom_type_dict.get((atom_conb, "cluster"),0),
        "sum": count
    }
atoms_and_type_info["Total"] = {
    "natom": None,
    "bulk":cell_type_dict.get("bulk", 0),
    "layer":cell_type_dict.get("layer", 0),
    "cluster":cell_type_dict.get("cluster", 0),
    "sum":size
}
atoms_df = pd.DataFrame(atoms_and_type_info).T # need transpose
atoms_df.sort_index()
atoms_df.to_csv("Traindata_analysis_table.csv")

print(info_string)
