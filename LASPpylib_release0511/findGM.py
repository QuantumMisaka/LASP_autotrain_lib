# to find global minimum 

from allstr_new import AllStr
import sys

strfile = "all.arc"
limit = 1000
outputfile = "mini_sorted_str.arc"

if __name__ == "__main__":
    nargvs = len(sys.argv)
    if nargvs == 2:
        strfile = sys.argv[1]
    if nargvs == 3:
        strfile, limit = sys.argv[1:]
    if nargvs == 4:
        strfile, limit, outputfile = sys.argv[1:]

    limit = int(limit)    
    all_str = AllStr()
    all_str.arcinit([0,0],strfile,"")
    all_str_sorted = all_str.sort_by_energy()
    all_str_sorted.gen_arc(range(limit), outputfile,)
    
