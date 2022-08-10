# cut TrainStr and TrainFor by Natom
# JamesBourbon in 20220617
# also can use to sort-by-natom
from allstr_new import AllStr
from structure_new import Str
import sys

# setting
MIN_NATOM = 30
MAX_NATOM = 1000
STR_OUT = f"TrainStr_{MIN_NATOM}_{MAX_NATOM}.txt"
FOR_OUT = f"TrainFor_{MIN_NATOM}_{MAX_NATOM}.txt"


# read TrainStr and TrainFor
allstr_raw = AllStr()
if len(sys.argv) == 2:
    allstr_raw.train_data_init(sys.argv[1])
elif len(sys.argv) > 2:
    allstr_raw.train_data_init(sys.argv[1], sys.argv[2])
else:
    allstr_raw.train_data_init('TrainStr.txt', "TrainFor.txt")    

# get slice
db_size = len(allstr_raw)
allstr_sorted = allstr_raw.sort_by_natom()
index_min = 0
index_max = 0
for index, struc in enumerate(allstr_sorted):
    struc: Str
    natom = struc.natom
    if (natom >= MIN_NATOM) and not index_min:
        index_min = index
    if (natom > MAX_NATOM):
        index_max = index
        break
else:
    index_max = db_size
# print slice
printlist = range(index_min, index_max)
allstr_sorted.gen_data_str(printlist, STR_OUT)
allstr_sorted.gen_data_for(printlist, FOR_OUT)
print(f"---- Cut {len(printlist)} Structure from {db_size} Train-data ----")
print(f"---- OutPut StrFile as {STR_OUT} ----")
print("DONE!")
