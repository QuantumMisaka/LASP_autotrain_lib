# for trainstr trainsfer to arc str file
# JamesBourbon in 20220726

import sys
from allstr_new import AllStr

str_file = "TrainStr.txt"
if len(sys.argv) >= 2:
    str_file = sys.argv[1]

target_str = AllStr()
target_str.train_data_init(str_file)
target_str.gen_arc(range(len(target_str)))

print("DONE!")
