# to find global minimum 
# JamesBourbon in 20220725 update

from allstr_new import AllStr
import sys

strfile = "allstr.arc"
limit = 100
outputfile = "mini_sorted_str.arc"

if __name__ == "__main__":
    nargvs = len(sys.argv)
    if nargvs == 2:
        strfile = sys.argv[1]
    if nargvs == 3:
        strfile, limit = sys.argv[1:]
    if nargvs == 4:
        strfile, limit, outputfile = sys.argv[1:]
        
    all_str = AllStr()
    all_str.arcinit([0,0], strfile, forfile="")
    all_str_sorted = all_str.sort_by_energy()
    all_str_sorted.gen_arc(range(limit), outputfile,)
    
