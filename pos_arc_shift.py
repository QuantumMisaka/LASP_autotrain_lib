# pos2arc and arc2pos script

__author__ = 'refined by JamesBourbon, zpliu learn from hsd Arcread.py'
import sys
# sys.path.append('/home11/Liugroup/tools/')
from allstr_new import AllStr

POSCAR = "POSCAR"
ARC_INPUT = "all.arc"

def wrong_input_exit():
    print( "Please use following syntax :" )
    print( " pos2arc format " )
    print( " format : required argument, being arc2pos / pos2arc " )
    sys.exit()

if __name__ == "__main__":
    
    if (len(sys.argv)) < 2 or (sys.argv[1] not in ["pos2arc", "arc2pos"]):
        wrong_input_exit()
    else:
        if len(sys.argv) > 2:
            if sys.argv[1] == "pos2arc":
                POSCAR = sys.argv[2]
            elif sys.argv[1] == "arc2pos":
                ARC_INPUT = sys.argv[2]
            else:
                wrong_input_exit()
    allstr_read = AllStr()
    if sys.argv[1] == "pos2arc":
# POSCAR to arc
        allstr_read.read_coord_set_from_POSCAR(POSCAR)
        allstr_read.print_str_all('outstr.arc')
    if sys.argv[1] == "arc2pos":
#  arc to POSCAR: can be multiple
        allstr_read.arcinit([0,0],ARC_INPUT)  
        for i in range(len(allstr_read)): 
            allstr_read.gen_POSCAR_VASP(i)

