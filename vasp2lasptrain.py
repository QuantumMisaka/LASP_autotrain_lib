# only used when calculated Single-DFT label
# transfer OUTCAR, CONTCAR to TrainStr.txt and TrainFor.txt
# should read OSZICAR to judge VASP-done normally
# JamesBourbon update in 20220726

import numpy as np
from functools import reduce
import os
import sys

NLEM_filter = 200 # 120
# for check single-dft normally done        
        
Eledict = {'H': 1,     'He': 2,   'Li': 3,    'Be': 4,   'B': 5,     'C': 6,     'N': 7,     'O': 8,
           'F': 9,     'Ne': 10,  'Na': 11,   'Mg': 12,  'Al': 13,   'Si': 14,   'P': 15,    'S': 16,
           'Cl': 17,   'Ar': 18,  'K': 19,    'Ca': 20,  'Sc': 21,   'Ti': 22,   'V': 23,    'Cr': 24,
           'Mn': 25,   'Fe': 26,  'Co': 27,   'Ni': 28,  'Cu': 29,   'Zn': 30,   'Ga': 31,   'Ge': 32,
           'As': 33,   'Se': 34,  'Br': 35,   'Kr': 36,  'Rb': 37,   'Sr': 38,   'Y': 39,    'Zr': 40,
           'Nb': 41,   'Mo': 42,  'Tc': 43,   'Ru': 44,  'Rh': 45,   'Pd': 46,   'Ag': 47,   'Cd': 48,
           'In': 49,   'Sn': 50,  'Sb': 51,   'Te': 52,  'I': 53,    'Xe': 54,   'Cs': 55,   'Ba': 56,
           'La': 57,   'Ce': 58,  'Pr': 59,   'Nd': 60,  'Pm': 61,   'Sm': 62,   'Eu': 63,   'Gd': 64,
           'Tb': 65,   'Dy': 66,  'Ho': 67,   'Er': 68,  'Tm': 69,   'Yb': 70,   'Lu': 71,   'Hf': 72,
           'Ta': 73,   'W': 74,   'Re': 75,   'Os': 76,  'Ir': 77,   'Pt': 78,   'Au': 79,   'Hg': 80,
           'Tl': 81,   'Pb': 82,  'Bi': 83,   'Po': 84,  'At': 85,   'Rn': 86,   'Fr': 87,   'Ra': 88,
           'Ac': 89,   'Th': 90,  'Pa': 91,   'U': 92,   'Np': 93,   'Pu': 94,   'Am': 95,   'Cm': 96,
           'Bk': 97,   'Cf': 98,  'Es': 99,   'Fm': 100, 'Md': 101,  'No': 102,  'Lr': 103,  'Rf': 104,
           'Db': 105,  'Sg': 106, 'Bh': 107,  'Hs': 108, 'Mt': 109,  'Ds': 110,  'Rg': 111,  'Cn': 112,
           'Nh': 113, 'Fl': 114, 'Mc': 115, 'Lv': 116, 'Ts': 117, 'Og': 118}

class Str():
    def __init__(self):
        self.energy = 0
        self.stress_vasp = [] # vasp-format stress in XX/YY/ZZ/XY/YZ/ZX; KB
        self.stress_lasp = []  # lasp-format stress in XX/XY/ZX/YY/YZ/ZZ: ev/A^3
        self.all_force = []
        self.Coord = []
        self.frac = []
        self.all_element = []
        self.all_ele_ind = []
        self.ele_names = []
        self.ele_inds = []
        self.element_count = {}
        # self.max_force = 0
        self.Cell = []
        self.Latt = []
        self.Natom = 0
    
    def element_to_list(self, element_count: dict):
        '''transfer element_count dict to all_element list'''
        element_list = []
        ele_indexs = []
        for element,count in element_count.items():
            for i in range(count):
                element_list.append(element)
                ele_indexs.append(Eledict[element])
        self.all_element = element_list
        self.all_ele_ind = ele_indexs

    def stress_vasp_to_lasp(self, stress_vasp: list):
        '''transfer vasp-format stress to lasp-format'''
        XX,YY,ZZ,XY,YZ,ZX = stress_vasp[:]
        KB_to_eVA3 = np.float64(1602.1766208)
        stress_reord = np.array([XX,XY,ZX,YY,YZ,ZZ],dtype=np.float64)
        return stress_reord / KB_to_eVA3
    
    def Frac_Coord(self):
        '''use Cart Coord and Cell info to get Frac Coord'''
        cellr = np.linalg.inv(self.Cell)
        return np.dot(self.Coord, cellr)


    def read_OUTCAR(self, filename="OUTCAR"):
        '''read OUTCAR to get force, stress and energy'''
        stress_tick = "in kB"
        # lattice_tick = "length of vectors"
        force_tick = "TOTAL-FORCE"
        energy_tick = "free  energy   TOTEN"
        # max_force_all = 0
        # test
        try:
            f = open(filename, 'r')
        except:
            print("----No OUTCAR info----")
            return None
        else:
            f.close()
        # read OUTCAR, get max-force, energy and allfor.arc
        with open(filename, 'r') as fo:
            print("read OUTCAR")
            for line in fo:
                # get stress
                if stress_tick in line:
                    stress_list = line.strip().split()
                    self.stress_vasp = -1 * np.array(stress_list[-6:], dtype=np.float64)
                    # print(stress_vasp)
                    self.stress_lasp = self.stress_vasp_to_lasp(self.stress_vasp)
                if force_tick in line:
                    space_line = fo.readline()
                    # print(space_line)
                    data_list = fo.readline().split()
                    while len(data_list) == 6:
                        one_atom_force = np.array(data_list[-3:],dtype=np.float64)
                        # max_force_one = np.max(np.sqrt(np.square(one_atom_force)))
                        # max_force_all = np.max([max_force_one, max_force_all])
                        self.all_force.append(one_atom_force)
                        data_list = fo.readline().split()
                    # self.max_force = max_force_all
                if energy_tick in line:
                        energy_list = line.split()
                        self.energy = np.float64(energy_list[4])
        return
                        
    def read_CONTCAR(self, filename="CONTCAR"):
        '''read CONTCAR(POSCAR) to get cell and atom coord
        
        in single_calc, recommended to read POSCAR instead
        '''
        try:
            f = open(filename, 'r')
        except:
            print('----No CONTCAR info----')
        else:
            print("read CONTCAR")
            index = 0
            atom = 0
            coord_type_status = False
            Cart = 0
            for line in f:
                L_list = line.strip().split()
                index += 1
                if index > 2 and index < 6:
                    # read Cell info
                    self.Cell.append([np.float64(x) for x in L_list])
                if index == 6:
                    # read element name
                    # self.Latt = self.Cell2Latt(Cell=self.Cell)
                    elements = L_list
                if index == 7:
                    # read element count and get all_element list
                    for i in range(len(L_list)):
                        self.ele_names.append(elements[i])
                        self.ele_inds.append(Eledict[elements[i]])
                        self.element_count[elements[i]] = int(L_list[i])
                    self.element_to_list(self.element_count)
                # problem: something will add line
                if index > 7 and coord_type_status == False:
                    # read coord type
                    if (L_list[0][0]=="C" or L_list[0][0]=="c"): 
                        Cart=1 # Cartesian coord
                        coord_type_status = True
                    if (L_list[0][0]=="D" or L_list[0][0]=="d"): 
                        Cart=0# frac coord
                        coord_type_status = True
                if coord_type_status == True:
                    # read atom coord
                    if L_list == []: break # pick the end
                    if len(L_list) < 3: continue
                    if Cart == 1:
                        self.Coord.append([np.float64(x) for x in L_list[0:3]])
                        self.frac = self.Frac_Coord()
                    if Cart == 0:
                        self.frac.append([np.float64(x) for x in L_list[0:3]])
                        self.Coord.append(np.dot([np.float64(x) for x in L_list[0:3]], self.Cell))
                    atom += 1    
                
            self.Natom = atom    
            f.close()
        return    
        
        
    def gen_trainstr(self, fname="TrainStr.txt"):
        '''directly print one-str TrainStr.txt'''    
        with open(fname, 'w') as fout:
            fout.write(" Start one structure\n")
            fout.write("Energy is %12.6f eV\n"%self.energy)
            fout.write("total number of element is %5d\n"%(self.Natom))
            fout.write("element in structure:\n")
            fout.write("symbol %s\n" %reduce(lambda a,b:a+b , ["%4s"%s   for s   in self.ele_names] ))
            fout.write("No.    %s\n" %reduce(lambda a,b:a+b , ["%4d"%num for num in self.ele_inds] ))
            fout.write("number %s\n" %reduce(lambda a,b:a+b , ["%4d"%num for num in self.element_count.values()]))
            fout.write("weight    1.000    1.000\n") # weight for train?
            for lat in self.Cell:
                fout.write("lat %15.8f  %15.8f  %15.8f\n"%(lat[0], lat[1], lat[2]))
            for i in range(self.Natom):
                atom_ele_ind = self.all_ele_ind[i]
                atom_xyz = self.Coord[i]
                atom_charge = 0 
                fout.write("ele %4s %15.8f  %15.8f  %15.8f  %15.8f\n"%(
                    atom_ele_ind, atom_xyz[0], atom_xyz[1], atom_xyz[2], atom_charge))
            fout.write(" End one structure\n\n")
    
    def gen_trainfor(self, fname="TrainFor.txt"):
        '''directlt print one-str TrainFor.txt'''
        with open(fname, 'w') as fout:
            stress = self.stress_lasp
            fout.write(" Start one structure\n")
            fout.write("stress %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f\n"%(
                stress[0], stress[1], stress[2], stress[3], stress[4], stress[5]))
            for i in range(self.Natom):
                atom_ele_ind = self.all_ele_ind[i]
                atom_force = self.all_force[i]
                fout.write("force %4d %15.8f %15.8f %15.8f\n"%(
                    atom_ele_ind, atom_force[0], atom_force[1], atom_force[2]))
            fout.write(" End one structure\n\n")

def is_VASP_done(NLEM_limit=NLEM_filter):
    '''check VASP-calc. is done or not by read OSZICAR'''
    if os.path.exists("OSZICAR"):
        print("read OSZICAR")
        finish_line = os.popen("grep F OSZICAR").readline()
        icontrol = int(os.popen("cat OSZICAR | wc -l").readline().strip())
        # print(finish_line)
        # print(icontrol)
        if "E0" in finish_line:
            if icontrol <= NLEM_limit:
                return True
        print("VASP-calc. not normally done")
        return False
    else:
        print("OSZICAR not found")
        return False
    
if __name__ == "__main__":
    if len(sys.argv) >= 2:
        NLEM_filter = sys.argv[1]
    ROOTDIR = os.getcwd()
    print("Transfer VASP-OUTPUT to LASP-TrainData")
    print("need OSZICAR, OUTCAR, CONTCAR")
    print("Run in %s"%ROOTDIR)
    if is_VASP_done():
        one_str = Str()
        one_str.read_OUTCAR()
        one_str.read_CONTCAR()
        one_str.gen_trainfor()
        one_str.gen_trainstr()
        print("done")
    else:
        print("Transfer not properly run")
  

            

