# JamesBourbon in 20220814
# add coordination pattern calc function
from atom_k import S_atom 
from bond_k import bond # not useful in LASP-autotrain
import numpy as np
import PeriodicTable as PT
import math as m
import os
import ctypes
import re
from functools import reduce # in py3
import pandas as pd
# from coordination_pattern import CoordinationPatterns


# variable likes to be Str object
def wrapCalBondMatrix(Str):
    return Str.Bondmatrix()

def wrapSegmolecular(Str):
    return Str.Segmolecular()

def wrapCalFakebmx(Str):
    Lminstr, bmx,bondneed,surface = Str.CheckMin()
    return Lminstr,bmx,bondneed,surface


class Str(object):
    '''define Structure object
    '''
    def __init__(self):
        self.atom = [] # S_atom list for all atom in Str
        self.energy = 0 # Energy of a Str
        self.abc= [] # lattice parameter
        self.Latt = [] # lattice parameter
        self.Cell = [] # abc-in-xyz, POSCAR-abc, transfer-matrix
        self.natom = 0 # number of atoms in Str
        self.Ele_Name =[] # element name list for all atom in str
        self.Ele_Index  = [] # element index of all atom in str
        # Ele_Index are also named iza in ZPLiu group？
        self.Coord = [] # Coordination of each atom
        self.sp ={} # atom and their number
        self.sporder= {} # atom ordered by element
        self.Cell = [] # lattice abc in Cartesian-3D
        self.natompe= []
        self.eleList =[] # specific element index in Str
        self.ele_nameList = [] # specific element in Str
        
        # relevant to force
        self.Lfor = True # have allfor.arc to read or not
        self.For= [] # list contains force of Str in each-atom-3D
        self.stress = [] # stress list ?

        # self.serial_num = 0
        self.maxF = 0
        self.max_stress = 0
        self.frac = []
        self.centerf = {}
        self._cycle=[]
        self.cart = []
        self.Q= []

    def addatom(self, line, flag = 1 ):
        '''add one atom from structure format file line
        
        Args:
            line(string): one line in structure file
            flag(int): structure file format flag, default is 1
                1 for arc and MS-xsd,
                2 for ?,
                3 for cat file,
                4 for mol file,
                5 for QM9,
                6 for ? ,
        '''
        
        if flag == 1:
            #for arc and ms
            coord = [float(x) for x in line.split()[1:4]]
            ele_symbol = line.split()[0] # get element symbol
            ele_number = PT.Eledict[ele_symbol] # get element number using symbol
            single_atom = S_atom(coord,ele_number) # 
            # get charge
            try:
                single_atom.charge = float(line.split()[-2])
            except :
                single_atom.charge = 0.0
            # append to atom list
            self.atom.append(single_atom)
            return
        # not important for LASP-VASP training
        elif flag == 2:
            # for TrainStr
            coord = [float(x) for x in line.split()[2:5]]
            ele_number = int(line.split()[1])
            self.atom.append(S_atom(coord,ele_number))
        elif flag == 3:# for cat file
            coord = [float(x) for x in line.split()[1:4]]
            ele_symbol = str(line.split()[0])[0]
            ele_number = PT.Eledict[ele_symbol] 
            self.atom.append(S_atom(coord,ele_number))
        elif flag == 4:#for mol file
            coord = [float(x) for x in line.split()[0:3]]
            ele_symbol = line.split()[3]
            ele_number = PT.Eledict[ele_symbol]
            self.atom.append(S_atom(coord,ele_number))
        elif flag == 5: #for QM9
            coord = [float(x) for x in line.split()[1:4]]
            ele_symbol = line.split()[0]
            ele_number = PT.Eledict[ele_symbol]
            self.atom.append(S_atom(coord,ele_number))
        elif flag == 6:
            coord = [float(x) for x in line.split()[1:4]]
            ele_symbol = line.split()[-2]
            ele_number = PT.Eledict[ele_symbol]
            single_atom= S_atom(coord,ele_number)
            single_atom.charge =float(line.split()[-1])
            self.atom.append(single_atom)



    def add_force(self, line, serial_num, flag=1):
        if flag == 1:
            # for arc
            self.atom[serial_num].force = [float(x) for x in line.split()]
        elif flag == 2:
            self.atom[serial_num].force= [float(x) for x in line.split()[2:5]]

    def add_stress(self, line, flag=1): # type 1 => for arc file; 2 => Data file
        if   flag == 1: self.stress = [float(x) for x in line.split()]
        elif flag == 2: self.stress = [float(x) for x in line.split()[1:7]]

    def add_charge(self, base):
        for atom in self.atom:
            atom.charge = base[atom.ele_symbol]

    def element_to_list(self, element_count: dict):
        '''transfer element_count dict (Str.sp) to all_element list'''
        element_list = []
        ele_indexs = []
        for element,count in element_count.items():
            for i in range(count):
                element_list.append(element)
                ele_indexs.append(PT.Eledict[element])
        all_element = element_list
        all_ele_ind = ele_indexs 
        return all_element, all_ele_ind
    
    def build_coord_set_from_POSCAR(self, filename="POSCAR"):
        '''read POSCAR to get cell and atom coord to build a Str
    '''
        try:
            f = open(filename, 'r')
        except:
            print('----No POSCAR info----')
        else:
            print("read POSCAR")
            index = 0
            atom_id = 0
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
                        self.ele_nameList.append(elements[i])
                        self.eleList.append(PT.Eledict[elements[i]])
                        self.sp[elements[i]] = int(L_list[i])
                    # get Ele_Name and Ele_Index
                    self.Ele_Name, self.Ele_Index = self.element_to_list(self.sp)
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
                        cart_coord = [np.float64(x) for x in L_list[0:3]]
                        self.Coord.append(cart_coord)
                        self.frac = self.FracCoord()
                    if Cart == 0:
                        frac_coord = [np.float64(x) for x in L_list[0:3]]
                        cart_coord = np.dot([np.float64(x) for x in L_list[0:3]], self.Cell)
                        self.frac.append(frac_coord)
                        self.Coord.append(cart_coord)
                    self.atom.append(S_atom)
                    self.atom[atom_id].xyz = cart_coord
                    self.atom[atom_id].ele_num = self.Ele_Index[atom_id]
                    self.atom[atom_id].ele_symbol = self.Ele_Name[atom_id]
                    atom_id += 1    
            self.natom = atom_id    
            f.close()
        return  


    def get_max_force(self):
        maxf= 0
        for atom in self.atom:
            _tmpf =  max([abs(x) for x in atom.force])
            if maxf < _tmpf:
                maxf = _tmpf
        self.maxF = maxf
        self.maxF = maxf
        
    def get_max_stress(self):
        self.max_stress = max(self.stress)
        



    def get_atom_num(self, ):
        '''used in VASP.run and train_str_init
        
        for get: natom, eleList, natompe, nele, ele_nameList
        '''
        self.natom = len(self.atom)
        self.natom = self.natom
        self.eleList,self.natompe = list(np.unique([atom.ele_num for atom in self.atom],return_counts=True))
        self.ele_nameList = [PT.Eletable[ele_num-1] for ele_num in self.eleList]
        self.sp = {}
        self.sporder = {}
        for index, ele_name in enumerate(self.ele_nameList):
            order_count = 1
            self.sp[ele_name] = self.natompe[index]
            self.sporder[ele_name] = order_count
            order_count += 1
        self.nele = len(self.eleList)
        self.Ele_Index = [atom.ele_num for atom in self.atom]
        self.Ele_Name  = [atom.ele_symbol for atom in self.atom]

    def screen_upper_surf(self,):
        self.upper =0
        for atom in self.atom:
            if atom.xyz[2] > 0.9*self.abc[2]:
                self.upper =1
                break

    def sort_atom_by_element(self,):
        '''used in VASP.run'''
        self.atom.sort(key =lambda X: X.ele_num)

    def sort_atom_by_z(self,):
        self.atom.sort(key= lambda X: X.xyz[2])

    def calc_two_dim_coord(self, ):
        self.xa = np.array([atom.xyz for atom in self.atom])

    def calc_one_dim_coord(self, ):
        self.cal_two_dim_coord()
        self.xa_onedim = reduce(lambda a,b: a+b, [list(x) for x in self.xa])


    def add_atom_ID(self):
        for iatom,atom in enumerate(self.atom):
            atom.id = iatom
 
    def cdnt2fcnt(self, ):
        '''set frac Coord used by self.fdnt'''
        self.cal_two_dim_coord()
        latinv    = np.linalg.inv(self.Cell)
        self.fdnt = [list(x) for x in np.matmul(self.xa, latinv)]
        # self.fdnt 

    def set_coord(self):
        '''set self.Coord from self.atom'''
        self.Coord = []
        for atom in self.atom:
            # atom: S_atom
            self.Coord.append(atom.xyz)
        self.Coord = np.array(self.Coord)
            

    def calc_centroid(self):
        # centroid = mass center
        self.get_atom_num()
        xyzall =np.array([0,0,0])
        for atom in self.atom:
            xyzall= xyzall+np.array(atom.xyz)
        self.centroid =xyzall/self.natom

    def mv_str_to_boxcenter(self):
        #for Ortholat
        self.calc_centroid()
        center =np.array([self.abc[0]/2,self.abc[1]/2,self.abc[2]/2])
        mv = center - self.centroid
        for i,atom in enumerate(self.atom):
            atom.xyz = list(self.xa[i]+mv)

    def add_ortho_lat(self):
        #self.calAtomnum()
        if self.natom != 0:
            self.cal_two_dim_coord()
            max3d=np.amax(self.xa, axis=0)
            min3d=np.amin(self.xa, axis=0)
            delta =max(list(max3d-min3d))
            a = 10
            while ((a-delta) < 4):
                a=a+5
            self.abc =[a,a,a,90,90,90]
            

    def abc2lat(self):
        a, b, c  = self.abc[0:3]
        alpha, beta, gamma = [x*np.pi/180.0 for x in self.abc[3:]]

        bc2 = b**2 + c**2 - 2*b*c*np.cos(alpha)
        h1 = a
        #h1= checkzero(h1)
        h2 = b * np.cos(gamma)
        # h2 = checkzero(h2)
        h3 = b * np.sin(gamma)
        #  h3 = checkzero(h3)
        h4 = c * np.cos(beta)
        # h4 = checkzero(h4)
        h5 = ((h2 - h4)**2 + h3**2 + c**2 - h4**2 - bc2)/(2 * h3)
        # h5 = checkzero(h5)
        h6 = np.sqrt(c**2 - h4**2 - h5**2)
        # h6 = checkzero(h6)
        
        self.Cell = [[h1, 0., 0.], [h2, h3, 0.], [h4, h5, h6]]

    def lat2abc (self, flag=False, inlat=False):
        if not flag: lat = self.Cell
        else:        lat = inlat

        self.nlat = np.array(self.Cell)
        a = np.linalg.norm(lat[0])
        b = np.linalg.norm(lat[1])
        c = np.linalg.norm(lat[2])
        alpha = m.acos(np.dot(lat[1],lat[2]) / (b*c))*180.0/np.pi
        beta  = m.acos(np.dot(lat[0],lat[2]) / (a*c))*180.0/np.pi
        gamma = m.acos(np.dot(lat[0],lat[1]) / (a*b))*180.0/np.pi
        if not flag: self.abc = [a,b,c,alpha, beta, gamma]
        else:        return [a,b,c,alpha, beta, gamma]


    def outPOSCAR(self,outfile="POSCAR"):
        '''print structure to POSCAR file'''
        f=open(outfile,'w')
        f.write('system\n')
        f.write('1.000000000000\n')
        for item in self.Cell:
            f.write('%12.8f  %12.8f  %12.8f\n'%(item[0],item[1],item[2]))
        f.write("%s\n" %reduce(lambda a,b:a+b , ["%4s"%s   for s   in self.ele_nameList]))
        f.write("%s\n" %reduce(lambda a,b:a+b , ["%4d"%num for num in self.natompe]  ))
        f.write('Cart\n')
        for atom in self.atom:
            f.write('%12.8f  %12.8f  %12.8f\n'%(atom.xyz[0],atom.xyz[1],atom.xyz[2]))
        f.close()

    def genPOTCAR(self,sourcedir,outfile="POTCAR"):
        '''get needed POTCAR from sourcedir, POTCAR will be print named outfile'''
        os.system('rm -rf %s'%outfile)
        for ele in self.ele_nameList:
            os.system('cat %s/POTCAR.%s >> %s'%(sourcedir,ele,outfile))


    def genKPOINTS(self,outfile="KPOINTS", criteria=25):
        '''get needed kpoints, ka = kb = kc = criteria
        
        writed by JamesBourbon
        '''
        ka = int(m.ceil(criteria / float(self.abc[0])))
        kb = int(m.ceil(criteria / float(self.abc[1])))
        kc = int(m.ceil(criteria / float(self.abc[2])))
        
        with open(outfile, 'w') as fout:
            fout.write('mesh auto\n')
            fout.write('0\n')
            fout.write('G\n')
            fout.write('%d %d %d\n'%(ka,kb,kc))
            fout.write('0 0 0')
        
    def set_element_radius(self, dataset_name = ''):
        """Setting element radius data in structure, by JamesBourbon
        
        Default radii dataset from 
        Chem. Eur. J. 2009, 15, 186–197, DOI: 10.1002/chem.200800987
        
        Args:
            dataset_name (str): Defaults using Eleradii in PeriodicTable.
            
        Returns: 
            radius_dict (dict): Ele_Index for key and Ele_radii for value
        """
        radius_dict = {}
        if bool(dataset_name):
            # read from external csv file if given
            try:
                radius_dataset = pd.read_csv(dataset_name)
                for ele_index in set(self.Ele_Index):
                    radii = np.float64(radius_dataset[
                        radius_dataset["index"] == ele_index ]["radius(A)"])
                    radius_dict[ele_index] = radii
            except:
                for ele_index in set(self.Ele_Index):
                    radii = PT.Eleradii[ele_index-1] # notice!
                    radius_dict[ele_index] = radii
        else:
            for ele_index in set(self.Ele_Index):
                radii = PT.Eleradii[ele_index-1]
                radius_dict[ele_index] = radii
        return radius_dict
        
    
    
    def calc_dis(self, iatom1, iatom2):
        '''calc periodic atom distance, considered periodic in x-np.round(x)
        
        just consider one cell'''
        self.cdnt2fcnt() # get fractional coord
        vbond =np.array(self.fdnt[iatom1])- np.array(self.fdnt[iatom2])
        dis = np.linalg.norm(np.matmul( np.array([x-np.round(x) for x in vbond]), self.Cell))
        return dis
    
    def calc_neighbour(self,iatom):
        '''calc distance of one atom to all other atom
        
        cannot consider trans-cell coordination.'''
        dict ={}
        for i in range(self.natom):
            if (i!= iatom):
                dis= self.calc_dis(iatom,i)
                dict[i]= dis
        return dict
    
    def make_supercell(self,x_max:int,y_max:int,z_max:int,
                    x_min:int=0, y_min:int=0, z_min:int=0):
        '''make supercell by x, y, z, running well
        
        returns:
            supercell: list of (ele_index, coordinate)
        '''
        supercell = []
        for i in range(x_min,x_max+1):
            for j in range(y_min,y_max+1):
                for k in range(z_min,z_max+1):
                    dims = [i,j,k]
                    # get move vector: dims@self.Cell, moving along a,b,c
                    moved_str = ((self.Ele_Index[index], 
                                    coord+np.matmul(dims,self.Cell)) 
                                    for index, coord in enumerate(self.Coord))
                    supercell += moved_str
        return supercell
            
    # need to be more : this just calc one Str, coordination pattern analysis should in allstr or other
    # coding in 20220609
    
    def coordination_pattern(self, tol=1.21,):
        '''calc coordination pattern of Str, interface to coordination_pattern.py
        
        coor_patterns example: {(Pd, ((Au, 2.4),)), (Au, ((Pd, 2.4),))}
        
        Returns:
            coor_patterns: coordination patterns format set
        '''
        # setting
        radius_dict = self.set_element_radius()
        # judge each dim is thin or not
        thin_edge_list = [max(radius_dict.values())*i*2*tol+0.2 for i in range(1,5)]
        shell_k_list = [1.2, 0.6, 0.4, 0.3, 0.25] # test result
        shell_k = 1.2
        # last k for large lattice parameter, cutoff ref: 1.0, 0.5, 0.34, 0.25, 0.20
        # Judge min bond length
        min_bond = min(PT.Eleradii) * tol
        # 3x3x3 supercell
        lat_xyz = np.diag(self.Cell) # describe a shell
        supercell_333 = self.make_supercell(1,1,1,-1,-1,-1)
        # distance calculated from central cell
        coor_patterns = set()
        for ci,central_coord in enumerate(self.Coord):
            central_ele = self.Ele_Name[ci]
            coor_pattern_dict = {}
            for atom_info in supercell_333:
                # atom_info: (ele_index, atom_xyz)
                do_calc = True
                square_dist_array = []
                for dim, vector in enumerate(lat_xyz):
                    # to simplify calculation: not to far in any dim
                    for i, edge in enumerate(thin_edge_list):
                        if vector < edge:
                            shell_k = shell_k_list[i]
                            break
                    else:
                        shell_k = shell_k_list[-1]
                    square_dist_dim = np.power(central_coord[dim] - atom_info[1][dim], 2)
                    if square_dist_dim > np.power(vector*shell_k, 2): 
                        do_calc = False
                        break
                    else:
                        do_calc = True
                        square_dist_array.append(square_dist_dim)
                # calc coordination main
                if do_calc:
                    bond_dist = np.sqrt(np.sum(square_dist_array))
                    central_index = self.Ele_Index[ci]
                    central_radii = radius_dict[central_index]
                    coor_radii = radius_dict[atom_info[0]]
                    # coor-bond main judgement
                    max_bond = (central_radii + coor_radii) * tol
                    if (bond_dist > max_bond) or (bond_dist < min_bond):
                        continue
                    else:
                        coor_ele = PT.Eletable[atom_info[0]-1]
                        coor_dist = np.round(bond_dist,1)
                        coor_pair = (coor_ele, coor_dist)
                        # get CN
                        coor_pattern_dict[coor_pair] = coor_pattern_dict.get(coor_pair,0) + 1
                    # end of coor_bond calc.
            # transfer to Set(tuple) for coor-patterns format
            coor_atoms = []
            for coor_pair, count in coor_pattern_dict.items():
                coor_atoms.append((coor_pair[0], coor_pair[1], count))
            coor_atoms = tuple(coor_atoms)
            coor_pattern_i = (central_ele, coor_atoms)
            coor_patterns.add(coor_pattern_i)
        return coor_patterns # can be directly used in CoordinationPatterns.patterns        
        
    # end of coordination pattern method embedding by JamesBourbon

    def special_neighbour(self,iatom,spe_ele):
        '''calc distance of one atom to all atom of one element''' 
        dict ={}
        for i in range(self.natom):
            if (i!= iatom and self.atom[i].ele_num == spe_ele):
                dis= self.calc_dis(iatom,i)
                dict[i]= dis
        return dict

    def longest_bond(self,ele1,ele2,lim):
        '''calc longest bond in cell'''
        self.cdnt2fcnt()
        maxlist = []
        for i,atom in enumerate(self.atom):
            if atom.ele_num == ele1:
                result = self.special_neighbour(i,ele2)
                dis = [x for x in result.values() if x<lim]
                if len(dis) > 0:
                    maxlist.append(max(dis))
        longest = max(maxlist)
        return longest


    def shortest_bond(self,ele1,ele2,lim):
        '''calc shortest bond in cell'''
        self.cdnt2fcnt()
        minlist = []
        for i,atom in enumerate(self.atom):
            if atom.ele_num == ele1:
                result = self.special_neighbour(i,ele2)
                dis = [x for x in result.values() ]
                if len(dis) > 0:
                    minlist.append(min(dis))
        short = min(minlist)
        return short
    
    

        

    def simple_class(self,iatom):
        for i in range(self.natom):
            if (i != iatom):
                #dislist.append(self.cal_distance(iatom,i))
                #print dis
                dis2 = self.calc_dis(iatom,i)
                #print dis
                bondtest = bond(self.atom[iatom].ele_num,self.atom[i].ele_num,dis2)
                bondorder = bondtest.judge_bondorder()
                if bondorder != 0:
                    if self.atom[i].ele_num == 8:
                        self.atom[iatom].bondtype.append(101)
                    if self.atom[i].ele_num == 6:
                        self.atom[iatom].bondtype.append(11)
                    if self.atom[i].ele_num == 1:
                        self.atom[iatom].bondtype.append(1)

        self.atom[iatom].bondtype.sort( )
        return

    def bond_search(self,iatom,ele2,flag=1):
        bondlist =[]
        for i in range(self.natom):
            if (i != iatom) and self.atom[i].ele_num ==ele2:
                if flag ==2 and i<iatom and self.atom[iatom].ele_num ==ele2: continue
                #dislist.append(self.cal_distance(iatom,i))
                dis = self.calc_dis(iatom,i)
                bondtest = bond(self.atom[iatom].ele_num,self.atom[i].ele_num,dis)
                bondorder = bondtest.judge_bondorder()
                if bondorder!= 0:
                    bondlist.append(dis)
        return bondlist
                
    
    def determine_species(self,iatom,elelist):
        #dislist = []
        for i in range(self.natom):
            if (i != iatom):
                #dislist.append(self.cal_distance(iatom,i))
                dis = self.calc_dis(iatom,i)
                bondtest = bond(self.atom[iatom].ele_num,self.atom[i].ele_num,dis)
                bondorder = bondtest.judge_bondorder()
                if bondorder != 0 :#and self.atom[i].ele_num == 8:
                    self.atom[iatom].bondlist.append(bondorder)
                    #for id,iza in enumerate(elelist):
                    #    if self.atom[i].ele_num ==iza:
                    #        if id >0:
                    #            self.atom[iatom].bondtype.append(bondorder+10**id)
                    #        else:
                    #            self.atom[iatom].bondtype.append(bondorder)
                    #        break

                    if self.atom[i].ele_num == elelist[3]:
                        self.atom[iatom].bondtype.append(bondorder+1000)
                    if self.atom[i].ele_num == elelist[2]:
                        self.atom[iatom].bondtype.append(bondorder+100)
                    if self.atom[i].ele_num == elelist[1]:
                        self.atom[iatom].bondtype.append(bondorder+10)
                    if self.atom[i].ele_num == elelist[0]:
                        self.atom[iatom].bondtype.append(bondorder)
        self.atom[iatom].species= len(self.atom[iatom].bondlist)

        bond_all = sum(self.atom[iatom].bondlist)
        self.atom[iatom].bondtype.sort(reverse = True )
        if bond_all != 4:
            self.atom[iatom].Ctype = 'radical'
        else :
            if self.atom[iatom].bondtype[-1] == 103 :
                self.atom[iatom].Ctype = 'CO'
            elif self.atom[iatom].bondtype[-1] == 102 :
                if self.atom[iatom].bondtype[-2] == 101 :
                    self.atom[iatom].Ctype = 'acid'
                else:
                    self.atom[iatom].Ctype = 'ketone'
            elif self.atom[iatom].bondtype[-1] == 101 :
                self.atom[iatom].Ctype = 'alcohol'
            elif self.atom[iatom].bondtype[-1] == 13 :
                self.atom[iatom].Ctype = 'alkyne'
            elif self.atom[iatom].bondtype[-1] == 12 :
                self.atom[iatom].Ctype = 'alkene'
            elif self.atom[iatom].bondtype[-1] == 11 :
                self.atom[iatom].Ctype = 'alkane'
            elif self.atom[iatom].bondtype[-1] == 1 :
                self.atom[iatom].Ctype = 'CH4'

                #self.atom[iatom].species =self.atom[iatom].species+ bondorder
        #dislist.sort()
        #length = min(len(dislist),4)
        #for i in range(length):
        #    if dislist[i]> 1.7:
        #        break
        #self.atom[iatom].species = i + 1
        self.atom[iatom].bondall = bond_all
        return

    def calc_Ctypes(self, ):
        self.calc_one_dim_coord()
        cell = reduce(lambda a,b: a+b, self.Cell)
        #print cell
        self.c_natm = pointer(ctypes.c_int(self.natom))
        self.c_xa   = pointer((ctypes.c_double*len(self.xa_onedim))(*self.xa_onedim))
        self.c_rv   = pointer((ctypes.c_double*len(cell))(*cell))
        iza = [atom.ele_num for atom in self.atom]
        self.c_iza  = pointer((ctypes.c_int*self.natom)(*iza))
    

    def JudgeBond(self):
        return

    def calc_frag_charge(self):
        "charge here is fakecharge, actually group id"
        groupdict ={}
        for i,atom in enumerate(self.atom):
            #print atom.charge
            try:
                groupdict[atom.charge].append(atom)
            except:
                groupdict[atom.charge] =[]
                groupdict[atom.charge].append(atom)
        self.frag= groupdict
       
    def cal_group_frag(self):
        groupdict ={}
        for i,atom in enumerate(self.atom):
            #print atom.charge
            try:
                groupdict[self.group[i]].append(atom)
            except:
                groupdict[self.group[i]] =[]
                groupdict[self.group[i]].append(atom)
        self.frag= groupdict

 
    def determine_charge(self):
        self.cal_group_frag()
        for item in self.frag.values():
            _tmpcharge= 0
            _elelist =[]
            for atom in item:
                _elelist.append(atom.ele_num)
                _tmpcharge =_tmpcharge + atom.ele_num
            #print _elelist
            if _tmpcharge%2 != 0:
                for atom in item:
                    if atom.ele_num== 7:
                        if _elelist ==[7,8]:
                            return 15
                        if  _elelist ==[7,8,8]:
                            return 23
                        return -1
                return 1
        return 2

    def check_min(self,flag = 1):
        sqnatm = self.natom**2

        self.calc_Ctypes()
        program='/home7/kpl/pymodule/Lib/Lib_fillbond/checkminbond.so'
        Lminstr = pointer(ctypes.c_bool(0))
        bmatrix = pointer((ctypes.c_int*sqnatm)(*[0 for i in range(sqnatm)]))
        bondneed = pointer((ctypes.c_int*(self.natom))(*[0 for i in range(self.natom)]))
        surface = pointer((ctypes.c_int*(self.natom))(*[0 for i in range(self.natom)]))

        checkmin = ctypes.cdll.LoadLibrary(program)
        if flag == 0:
            checkmin.judgebond_(self.c_natm,self.c_iza,self.c_xa,self.c_rv,Lminstr,bmatrix,bondneed)
        elif flag ==1 :
            checkmin.judgebondsurface_(self.c_natm,self.c_iza,self.c_xa,self.c_rv,Lminstr,bmatrix,bondneed,surface)


        bmx = list(bmatrix.contents)
        #self.bmx2D = np.array(bmx).reshape(self.natom, self.natom)
        #self.bmx1D = bmx

        self.bondneed = list(bondneed.contents)
        self.Lminstr = bool(Lminstr.contents)
        return self.Lminstr,bmx,list(bondneed.contents),list(surface.contents)


    def bond_matrix(self, ):
        program='/home7/kpl/pymodule/Lib/Lib_bondmatrix/bondmatrix.so'
        self.calc_Ctypes()
        sqnatm = self.natom**2
        self.cdnt2fcnt()
        onedim_fa = reduce(lambda a,b: a+b, self.fdnt)
        self.c_fxa = pointer((ctypes.c_double*len(onedim_fa))(*onedim_fa))

        bmatrix = pointer((ctypes.c_int*sqnatm)(*[0 for i in range(sqnatm)]))
        bcal = ctypes.cdll.LoadLibrary(program)
        bcal.hbondmatrix_ (self.c_natm, self.c_fxa, self.c_iza, self.c_rv, bmatrix)
        return list(bmatrix.contents)

    def seg_molecular (self, recal=True):
        """ warning : bond matrix must first be calculated."""
        program='/home7/kpl/pymodule/Lib/Lib_bondmatrix/bondmatrix.so'
        if recal:
            self.calc_Ctypes()
            self.cdnt2fcnt()
            onedim_fa = reduce(lambda a,b: a+b, self.fdnt)
            self.c_fxa = pointer((ctypes.c_double*len(onedim_fa))(*onedim_fa))
        self.c_bmx = pointer((ctypes.c_int*len(self.bmx1D))(*self.bmx1D))
        group = pointer((ctypes.c_int*self.natom)(*[0 for i in range(self.natom)]))
        scal = ctypes.cdll.LoadLibrary(program)
        scal.hsegmentmol_ (self.c_natm, self.c_fxa, self.c_iza, self.c_rv, group, self.c_bmx)
        return list(group.contents)

    def get_pattern_atom(self,pattern):
        Atom1 =pattern[0]
        Atom2 =pattern[1]
        mode = int(pattern[-1])

        atomlist1 =[]
        design = re.findall(r"\d+\.?\d*",Atom1)
        if len(design)!= 0:
            for item in design:
                atomlist1.append(int(item)-1)
        else:
            for iatom,atom in enumerate(self.atom):
                if atom.ele_symbol == Atom1:
                    atomlist1.append(iatom)

        atomlist2 =[]
        design = re.findall(r"\d+\.?\d*",Atom2)
        if len(design)!= 0:
            for item in design:
                atomlist2.append(int(item)-1)
        else:
            atomlist2 =[]
            for iatom,atom in enumerate(self.atom):
                if atom.ele_symbol == Atom2:
                    atomlist2.append(iatom)

        return atomlist1,atomlist2,mode


    def get_all_bond(self):
        self.allbond=[]
        for i in range(self.natom):
            for j in xrange(i+1,self.natom):
                if self.bmx2D[i][j] > 0:
                    self.allbond.append([i,j,self.bmx2D[i][j]])

        self.nbond = len(self.allbond)

    def get_atom_info(self):
        for i in range(self.natom):
            self.atom[i].expbond= 0
            self.atom[i].imph= 0

        for i in range(self.natom):
            for j in range(i+1,self.natom):
                if self.bmx2D[i][j] > 0:
                    if self.atom[j].ele_num !=1:
                        self.atom[i].expbond = self.atom[i].expbond +1
                    else:
                        self.atom[i].imph = self.atom[i].imph +1
                    if self.atom[i].ele_num !=1:
                        self.atom[j].expbond = self.atom[j].expbond +1
                    else:
                        self.atom[j].imph = self.atom[j].imph +1


    def calc_unsaturated_number(self):
        nH =0
        needH = 2
        num_heavyatom = 0
        for atom in self.atom:
            if atom.ele_num == 1 :
                nH = nH +1
            else:
                num_heavyatom = num_heavyatom +1 
                needH= needH + (8-atom.ele_num)
        self.UnsNum =(needH -nH)/2
        self.nheavyatom = num_heavyatom


    def remove_metal_bond(self):
        self.surface_atom = [0 for i in range(self.natom)]
        self.bondneed = [0 for i in range(self.natom)]
        self.bmxsave =np.zeros((self.natom,self.natom),dtype =np.int)

        for i in range(self.natom):
            for j in range(self.natom):
                if self.bmx2D[i][j] > 0:
                    self.bmxsave[i][j] = self.bmx2D[i][j]
                    if self.atom[i].ele_num > 18 :#or self.atom[j].ele_num > 18:
                        self.bmx2D[i][j] =0
                        self.surface_atom[j]= 1
                    if self.atom[j].ele_num > 18:
                        self.bmx2D[i][j] =0
                        self.surface_atom[i]= 1

        self.bmx1D = []
        for line in self.bmx2D:
            self.bmx1D.extend(line)


        
        

#    def GetNeighbour(self,iatom):
#        for i in range(self.natom):
#            if self.bmx2D[iatom][i]>0:
#                self.atom[iatom].neighbour.append()
#
#    def GetAtomnodes(self):
#
#        for iatom,atom in enumerate(self):
#            self[i]

#
#"=============== the following part is designed for transfer ============="
#"must excute transfer before calling the different part function"
#
    def transfer_toP_XYZ_coordStr(self):
        self.natom = self.natom
        self.Latt = self.abc
        self.Cell = self.Latt2Cell()
        self.energy = self.energy
        self.Ele_Name =[]
        self.Ele_Index  = []
        self.Coord = []
        for atom in self.atom:
            self.Ele_Name.append(atom.ele_symbol)
            self.Ele_Index.append(atom.ele_num)
            self.Coord.append(atom.xyz)
    
        for iele,elesymbol in enumerate(self.ele_nameList):
            self.sp[elesymbol] = self.natompe[iele]
            self.sporder[elesymbol] = iele+1
    
    
        if self.Lfor :
            self.For= []
            for atom in self.atom:
                self.For.append(atom.force)
            self.get_max_force()
            self.maxF = self.maxF
            self.maxF = self.maxF
    
    
    def TransferToKplstr(self):
        # many many questions about why should do this function
        self.atom =[] # why define it again?
        for i in range(self.natom):
            self.atom.append(S_atom(self.Coord[i],self.Ele_Index[i]))
            if self.Lfor:
                self.atom[i].force= self.For[i]
        self.abc= self.Latt
        # self.energy= self.energy
        self.sort_atom_by_element()
        self.get_atom_num()
        self.abc2lat()

#
#"=============== the following part is copyed from zpliu XYZCoord.py ============="
#"must excute TransferToXYZcoordStr before calling the following function"
#

    def printCell(self):
        """ formated cell print---used to generate INPUT_DEBUG"""
        s = ""
        for x in self.Cell: s += "%10.4f %10.4f %10.4f\n  "%(x[0],x[1],x[2])
        return s

    def printCoord(self):
        """ used to generate INPUT_DEBUG"""
        s = ""
        for i,xa in enumerate(self.Coord):
            s += " %15.9f  %15.9f  %15.9f   %3d\n"%(xa[0],xa[1],xa[2],self.sporder[self.Ele_Name[i]])
        return s

    def myfilter(self,BadStrP):
        if( \
           self.energy > BadStrP.HighE*float(self.natom) or \
           self.energy < BadStrP.LowE*float(self.natom) or \
           min(self.Latt[:3]) < BadStrP.MinLat or \
           max(self.Latt[:3]) > BadStrP.MaxLat or \
           max(self.Latt[3:7]) > BadStrP.MaxAngle or \
           min(self.Latt[3:7]) < BadStrP.MinAngle or \
##          max(self.Latt[:3]) > 5.0*min(self.Latt[:3]) or \
           self.maxF > BadStrP.MaxFor):
           #print '111',self.maxF,BadStrP.MaxFor
           return False
        else:
          #print self.energy,True
         # return True
           if len(self.For) :
         #    #print self.maxF,BadStrP.MaxFor
              if self.maxF > BadStrP.MaxFor : 
                  print(self.maxF, BadStrP.MaxFor)
                  return False
              else: return True
           else:
              return True

    def myfilter_byd(self,Badd1,Badd2,type):

        L = True
        x1 = [x for x in Badd1.keys()]
        y1 = [y for y in Badd2.keys()]
#        print 'XXX',x1, y1
        for key,value in self.d.items():
#           print key, value, Badd1[key], Badd2[key]
            if key in x1:
                if type==1 and value < Badd1[key] : L = False ; continue
            if key in y1:
                if type==1 and value > Badd2[key] : L = False ; continue
            if key in x1 and key in y1:
                if type ==2 and value > Badd1[key] and value < Badd2[key] : L =  False
        return L


    def myfilter_byshortd(self,Badd1):

        L = True
        for key,value in self.d.items():
            if value < Badd1 : L = False ; continue
        return L



    def Latt2Cell(self):
        """ transform of (a,b,c,alpha,beta,gamma) to lattice vector of xyz
        
        tips1: lattice vector is POSCAR format Coord 
        
        tips2: Matrix of lattice vector is coord_translation matrix
        """
        # lower-triangle
        a,b,c,alpha,beta,gamma = self.Latt[:6]
        pi = 3.14159265358937932384626
        alpha,beta,gamma = alpha*pi/180.0,beta*pi/180.0,gamma*pi/180.0
        bc2 = b**2 + c**2 - 2*b*c*np.cos(alpha)
        h1 = a
        h2 = b * np.cos(gamma)
        h3 = b * np.sin(gamma)
        h4 = c * np.cos(beta)
        # h5 = ((h2 - h4)**2 + h3**2 + c**2 - h4**2 - bc2)/(2 * h3)
        h5 = c * (np.cos(alpha)-np.cos(beta)*np.cos(gamma))/np.sin(gamma)
        # modified in 20220422 by JamesBourbon
        # x = c**2 - h4**2 - h5**2
        h6 = np.sqrt(c**2 - h4**2 - h5**2)
        cell = [[h1, 0., 0.], [h2, h3, 0.], [h4, h5, h6]]
        return cell

    def Cell2Latt(self):
        """ transform of lattice vector to (a,b,c,alpha,beta,gamma)"""
        lat= self.Cell[0:3]
        a = np.linalg.norm(lat[0])
        b = np.linalg.norm(lat[1])
        c = np.linalg.norm(lat[2])
        pi = 3.141592653589793238
        alpha = m.acos(np.dot(lat[1],lat[2]) / (b*c))*180.0/pi
        beta  = m.acos(np.dot(lat[0],lat[2]) / (a*c))*180.0/pi
        gamma = m.acos(np.dot(lat[0],lat[1]) / (a*b))*180.0/pi
        return [a,b,c,alpha,beta,gamma]

    def Volume(self):
        """ calculate volume of structure"""
        a = self.Cell[0]
        b = self.Cell[1]
        c = self.Cell[2]
        return np.dot(np.cross(a,b),c)

    def ReciCell(self):
        """ reciprocal lattice """
        return np.transpose(np.linalg.inv(self.Cell))

    def FracCoord(self):
        """ get fractional coordinate """
        cellr = np.linalg.inv(self.Cell)
        return np.dot(self.Coord,cellr)

    def frac_module(self):
        frac = self.FracCoord()
        for i,x in enumerate(frac):
            frac[i]=map(lambda y:(y+1000.0) % 1.0,frac[i])
        return frac

    def centralize(self,n=1):
        if n==1 : frac = self.frac_module()
        else:
#           key= self.centerf.keys() ; key.sort()       #  sort by keys in dict
#           frac = np.array([self.centerf[i] for i in key])
            frac = [self.centerf[i] for i in sorted(self.centerf.keys())]
#       print frac
        cart = np.dot(frac,self.Cell)
        self.central_pos = np.dot([0.5,0.5,0.5],self.Cell)
        com = []
        for i in range(3):
            com.append(sum([x[i] for x in cart])/float(self.natom))
        self.cart = np.array([np.add(np.subtract(x,com),self.central_pos) for x in cart])
        self.frac = np.dot(self.cart,np.linalg.inv(self.Cell))
#       return self.frac

    def frac_center(self,at):
        # recursive to position all fractional coordinate
        if len(self._cycle)==0: self.centerf[at.keys()[0]]=at.values()[0]; self._currentlist =[]; self._currentlist.append(at.keys()[0])
        self._cycle.append(at.keys()[0])

        if len(self.centerf)==self.natom : return True
        else:
            neig = self.neighbor(at)
            # get frac of atom at  neig ={1:[XX,XX,XX],2:[XX,XX,XX]}
            for i in neig.keys():
                if i not in self.centerf.keys():
                    self.centerf[i]=neig[i]
                    self._currentlist.append(i)
                    if len(self.centerf)==self.natom : return True
            for i in neig.keys():
                if i not in self._cycle:  self.frac_center({i:self.centerf[i]})
                if len(self.centerf)==self.natom : return True
            if at.keys()[0]==0 : return False


    def chemical_formula(self):
        # recursive to position all fractional coordinate
        Nmol=0; fragment = {}
        while len(self.centerf) != self.natom:
            for i in range(self.natom) :
                if i not in self.centerf.keys():
                    self._cycle = [];self._currentlist =[]
                    L = self.frac_center({i:self.FracCoord()[i]})
                    Nmol +=1; elelist={}; fragment[Nmol]=""
                    #print self._currentlist, len(self._currentlist)
                    for j in self._currentlist:
                        if self.Ele_Name[j] in elelist.keys() : elelist[self.Ele_Name[j]] +=1
                        else: elelist[self.Ele_Name[j]] = 1
                    for j in elelist.keys():
                        fragment[Nmol]=fragment[Nmol]+j+str(elelist[j])
                    if L:
                        self.centralize(0)
                        for j in range(self.natom) : self.Coord[j][0:3]= self.cart[j][0:3]
                        return fragment, self.Coord, self.Cell
                    
    # update 0814, need refinement
    def get_basic_shape(self, vac=5):
        '''find the basic cell shape
        
        Returns: bulk, layer, cluster
        '''
        if (self.abc[0] == self.abc[1] == self.abc[2]
            ) and (self.abc[3] == self.abc[4] == self.abc[5]
            ) and self.abc[0] >= 10:
            return "cluster"
        else:
            # not enough, need to place the cell in center
            self.set_coord()
            max_coord = np.max(self.Coord, axis=0)
            min_coord = np.min(self.Coord, axis=0)
            coord_space = max_coord - min_coord
            is_vacc = (coord_space - vac - self.abc[:3]) >= 0
            sum_vac_dim = np.sum(is_vacc)
            if sum_vac_dim >= 3:
                return "cluster"
            elif sum_vac_dim >=1:
                return "layer"
            else:
                return "bulk"
            


    def judge_shape(self,cut=2.6):
        """ find the shortest bond between neighboring cell to judge solid shape"""
        if len(self.centerf)==0: L = self.frac_center({0:self.FracCoord()[0]})
        if not L: return {"fragments":-1}

        vect={}; vect_1={};com=[]
        cellconnect={}
        self.centralize(0)

        V=[];vect={}; Vcom={}
        for i in range(-1,2):
            for j in range(-1,2):
                for k in range(-1,2):
                    if (i==0 and j==0 and k==0): continue
                    V.append([i,j,k])
                    com = np.dot(np.add([0.5,0.5,0.5],[i,j,k]),self.Cell)
                    Vcom[str(i)+str(j)+str(k)]=np.subtract(com,self.central_pos)

        line = False ; layer = False; solid = False; count =0
        normalvector =[]; cellconnect ={}
        for i0,i1,i2 in V :
            if not str(i0)+str(i1)+str(i2) in cellconnect.keys():
                cellconnect[str(i0)+str(i1)+str(i2)] = self.cell_bondconnect(cut,[i0,i1,i2])
            if cellconnect[str(i0)+str(i1)+str(i2)] :
                line = True
                layer_hole_detect=[]
                for j0,j1,j2 in V :
                    if i0==j0 and i1==j1 and i2==j2: continue
                    angle=self.vect_angle(Vcom[str(i0)+str(i1)+str(i2)],Vcom[str(j0)+str(j1)+str(j2)])
                    if abs(angle-90) < 70:
                        if not str(j0)+str(j1)+str(j2) in cellconnect.keys():
                            cellconnect[str(j0)+str(j1)+str(j2)] = self.cell_bondconnect(cut,[j0,j1,j2])
                        if cellconnect[str(j0)+str(j1)+str(j2)]:
                            layer_hole_detect.append(angle)
                            layer = True
                            VectorC = np.cross(Vcom[str(i0)+str(i1)+str(i2)],Vcom[str(j0)+str(j1)+str(j2)])
                            for i in range(len(normalvector)-1):
                                for j in range(i+1,len(normalvector)):
                                    angle=self.vect_angle(normalvector[i],normalvector[j])
                                    if angle > 40:
                                        return {"solid":7777}   # early break if two intercected layer identified
                            angle=0.0
                            for i in range(len(normalvector)): angle=max(angle,self.vect_angle(normalvector[i],VectorC))
                            if len(normalvector)>1 and angle < 10: break
                            normalvector.append(VectorC)
                            hole_detect=[]
                            for k0,k1,k2 in V :
                                if k0==i0 and k1==i1 and k2==i2: continue
                                if k0==j0 and k1==j1 and k2==j2: continue
                                angle=self.vect_angle(Vcom[str(k0)+str(k1)+str(k2)],VectorC)
                                if angle < 70:
                                    if not str(k0)+str(k1)+str(k2) in cellconnect.keys():
                                        cellconnect[str(k0)+str(k1)+str(k2)] = self.cell_bondconnect(cut,[k0,k1,k2])
                                    if cellconnect[str(k0)+str(k1)+str(k2)]:
                                        hole_detect.append(angle)
                            #------THIRD LOOP END----
                            #print hole_detect
                            if len(hole_detect) ==0 :
                                return {"layer":7777}
                            elif (max(hole_detect)-min(hole_detect)) < 15:
                                return {"solid-largehole":9999}
                            elif (max(hole_detect)-min(hole_detect)) < 35:
                                return {"solid-smallhole":8888}
                            else : return {"solid-dense":6666}

                #        else: continue  ------SECOND LOOP END----
                # the following info on layer should not be met due to the early return above
                if len(layer_hole_detect) ==0:
                    return {"line":7777}
                elif (max(layer_hole_detect)-min(layer_hole_detect)) < 15:
                    return {"layer-largehole":9999}
                elif (max(layer_hole_detect)-min(layer_hole_detect)) < 35:
                    return {"layer-smallhole":8888}
                else : return {"layer-dense":6666}
        #   else: continue   ---- FIRST LOOP END------

        if layer : return {"?layer?":6666}   # should not be useful
        elif line : return {"?line?":6666}   # should not be useful
        else: return {"cluster":6666}

    def vect_angle(self,V0,V1):
        a = np.dot(V0,V1)/m.sqrt(np.dot(V0,V0))/m.sqrt(np.dot(V1,V1))
        angle = 180.0
        if abs(abs(a)-1.0) >1e-6: angle = m.acos(a)*180.0/3.14159265
        if angle > 90.0 : angle = 180.0 - angle
        return angle

    def cell_bondconnect(self,cut,cell1,cell0=[0,0,0]):

       frac1=[];VectorA = []; N=self.natom
       for i in range(N): frac1.append(np.add(self.frac[i],cell1))
       cart1=np.dot(frac1,self.Cell)

       if cell0[0]!=0 or cell0[1]!=0 or cell0[2]!=0 :
           frac0=[]
           for i in range(N): frac0.append(np.add(self.frac[i],cell0))
           cart0=np.dot(frac0,self.Cell)
           pos=[]
           pos.append(sum([y[0] for y in cart0])/float(self.natom))
           pos.append(sum([y[1] for y in cart0])/float(self.natom))
           pos.append(sum([y[2] for y in cart0])/float(self.natom))
       else:
           fract0 = self.frac
           cart0  = self.cart
           pos    = self.central_pos

       d0=999
       for i in range(N):
           for j in range(N):
               b = np.subtract(cart1[i],cart0[j])
               d = m.sqrt(np.dot(b,b))
               if d< d0: d0= d
       return  d0 < cut


    def neighbor(self,at,cut=2.6):
        """ find the shortest bond neighbors of atom  """
        frac = self.frac_module()
        cart = np.dot(frac,self.Cell)
        cart0= np.dot(at.values()[0],self.Cell)
        N=self.natom   #len(cart)
        neig={};neig2={}
        for j in range(N):
            if j==at.keys()[0]: continue
            b = np.subtract(cart0,cart[j])
            d = m.sqrt(np.dot(b,b))
            if d < cut:
                neig[j] = frac[j]
                neig2[j]=d
                #print d,neig
            #print bonddict
        cut1 = [-1,0,1]
        for i0 in cut1:
            for i1 in cut1:
                for i2 in cut1:
                    if i0==0 and i1 ==0 and i2==0: continue
                    #print i0,i1,i2
                    frac2=[]
                    for i in range(N):
                        frac2.append(np.add(frac[i],[i0,i1,i2]))
                    cart2=np.dot(frac2,self.Cell)
                    for j in range(N):
                        b = np.subtract(cart2[j],cart0)
                        d = m.sqrt(np.dot(b,b))
                        if d < cut:
                            if j in neig.keys():
                                if neig2[j] >d : neig[j] = frac2[j]; neig2[j]=d
                            else: neig[j] = frac2[j]; neig2[j]=d

                            #print d,neig
        #neig = sorted(neig,keys=neigh.valueskeys())
        return neig


    def Shortestbond(self):
        """ find the shortest bond """
        frac = self.frac_module()
        cart = np.dot(frac,self.Cell)
        N=self.natom   #len(cart)
        bonddict={}
        if N==1:
            bondname=self.Ele_Name[0]+"-"+self.Ele_Name[0]
            bonddict[bondname]=10.0
            return bonddict
        for i in range(N):
            for j in range(N):
                bondname=self.Ele_Name[i]+"-"+self.Ele_Name[j]
#               if self.sporder[self.EleNam[i]]>self.sporder[self.EleNam[j]]: bondname=self.EleNam[j]+"-"+self.EleNam[i]
                if PT.Eledict[self.Ele_Name[i]]>PT.Eledict[self.Ele_Name[j]] : bondname=self.Ele_Name[j]+"-"+self.Ele_Name[i]
                bonddict[bondname]=10.0
                #return bonddict
        for i in range(N):
            for j in range(i+1,N):
                b=np.subtract(cart[i],cart[j])
                d = m.sqrt(np.dot(b,b))  #(cart[i]-cart[j]),(cart[i]-cart[j])))
                bondname=self.Ele_Name[i]+"-"+self.Ele_Name[j]
                #print bondname
#               if self.sporder[self.EleNam[i]]>self.sporder[self.EleNam[j]]: bondname=self.EleNam[j]+"-"+self.EleNam[i]
                if PT.Eledict[self.Ele_Name[i]]>PT.Eledict[self.Ele_Name[j]] : bondname=self.Ele_Name[j]+"-"+self.Ele_Name[i]
                if bondname not in bonddict: bonddict[bondname]=d
                else: bonddict[bondname] =min(d,bonddict[bondname])
                #print bonddict
# for cluster
#       return bonddict
# for cluster
        cut1 = [-1,0,1]
        for i0 in cut1:
            for i1 in cut1:
                for i2 in cut1:
                    if i0==0 and i1 ==0 and i2==0: continue
                    #print i0,i1,i2
                    frac2=[]
                    for i in range(N): frac2.append(np.add(frac[i],[i0,i1,i2]))
                    cart2=np.dot(frac2,self.Cell)
                    for i in range(N):
                        for j in range(N):
                            b = np.subtract(cart2[i],cart[j])
                            d = m.sqrt(np.dot(b,b)) #(cart2[i]-cart[j]),(cart2[i]-cart[j])))
                            bondname=self.Ele_Name[i]+"-"+self.Ele_Name[j]
                           # if self.sporder[self.EleNam[i]]>self.sporder[self.EleNam[j]]: bondname=self.EleNam[j]+"-"+self.EleNam[i]
                            if PT.Eledict[self.Ele_Name[i]]>PT.Eledict[self.Ele_Name[j]] : bondname=self.Ele_Name[j]+"-"+self.Ele_Name[i]
                            if bondname not in bonddict: bonddict[bondname]=d
                            bonddict[bondname] =min (d,bonddict[bondname])
                            #print bonddict['Ge-Ge']
        return bonddict


    def _Gen_arc(self,coord, fname='outtmp.arc'):
        with open(fname, 'w') as fout:
            fout.write("!BIOSYM archive 2\nPBC=ON\n")
            energy = self.energy
            i=0
            if not self.Lfor : fout.write("\t\t\t\tEnergy\t%8d        -0.0000  %17.6f\n"%(i+1,energy))
            else: fout.write("\t\t\t\tEnergy\t%8d        %10.6f  %17.6f\n"%(i+1,self[i].maxF,energy))

            fout.write("!DATE\n")
            lat = self.Latt[0:6]
            fout.write("PBC  %12.6f%12.6f%12.6f%12.6f%12.6f%12.6f\n" %(lat[0],lat[1],lat[2],lat[3],lat[4],lat[5]) )
            for j in range(self.natom):
                ele = self.Ele_Name[j]
                xa =  coord[j]
                fout.write("%-2s%18.9f%15.9f%15.9f CORE %4d %-2s %-2s   0.0000 %4d\n"%\
                     (ele,xa[0],xa[1],xa[2],j+1,ele,ele,j+1))
            fout.write("end\nend\n")

    def hkl_dspacking(self,hkl):

        hklLib = ctypes.cdll.LoadLibrary('./crystal-hkl-subroutines.so')
        onedimrv = reduce(lambda a,b:a+b, self.Cell)
        rv    = pointer((ctypes.c_double*len(onedimrv))(*onedimrv))
        hkl0 = pointer((ctypes.c_int*len(hkl))(*hkl))
        d= 0.00
        d_c = pointer((ctypes.c_double(d)))

        hklLib.plane_dspacing2_(rv,hkl0,d_c)

        d = float(str(d_c.contents).split('(')[1].split(')')[0])

        return d


    def SteinhartQ_cal(self):
        '''calculate OP2 OP4 OP6 of structure
        
        use Fortran_lib calQ.so make by calQ.f90
        '''

        q = ctypes.cdll.LoadLibrary('calQ.so')
        xa = reduce(lambda a,b:a+b, self.Coord)
        cell = reduce(lambda a,b:a+b, self.Cell)

        natom = pointer(ctypes.c_int(self.natom))
        el    = pointer(ctypes.c_double(self.energy))
        atom  = pointer(ctypes.c_int(0))
        sym   = pointer(ctypes.c_int(int(0)))
        za    = pointer((ctypes.c_int*len(self.Ele_Index))(*self.Ele_Index))
        coord = pointer((ctypes.c_double*len(xa))(*xa))
        rv    = pointer((ctypes.c_double*len(cell))(*cell))
        qglobal = pointer((ctypes.c_double*4)(*[0.0,0.0,0.0,0.0]))
        q.get_order_parameter_(natom, za, rv, coord, atom, qglobal, el, sym)
        qval = (qglobal.contents)[1:4]

        q.get_dis_weight_order_parameter_(natom, za, rv, coord, atom, qglobal, el, sym)
        qval2 = (qglobal.contents)[1:4]

        return (qval[0], qval[1], qval[2],qval2[0], qval2[1], qval2[2])


class FooError(Exception):
    pass

class BadStr(object):
    def __init__(self):
        self.HighE    = -0
        self.LowE     =-9999999
        self.MinAngle = 50
        self.MaxAngle = 130
        self.MinLat   = 1.6
        self.MaxLat   = 50.0
        self.MaxFor   = 10.0

def ParaWrap_JudgeShape(x):
    return x.JudgeShape()
def ParaWrap_Shortestbond(x):
    return x.Shortestbond()
def ParaWrap_ChemicalFormula(x):
    return x.ChemicalFormula()
def ParaWrap_SteinhartQ_cal(x):
    return x.SteinhartQ_cal()
def ParaWrap_Coor_Patterns(x):
    return x.coordination_pattern()



if __name__=='__main__':
    a=1
    print(a)
