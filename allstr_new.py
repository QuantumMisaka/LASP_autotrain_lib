# noting by JamesBourbon
# last change in 20220730
# from atom_k import S_atom
from structure_new import Str
from structure_new import wrapCalBondMatrix, wrapSegmolecular,wrapCalFakebmx, BadStr, FooError, ParaWrap_JudgeShape, ParaWrap_Shortestbond
from structure_new import ParaWrap_ChemicalFormula
import random
#import bondall
from multiprocessing import Pool
import PeriodicTable as PT
# from babelfunc import *
# zpliu part
import numpy as np
import os
import string
# import time
# import cmath as m
from functools import reduce


class AllStr(list):
    # element:Str
    # def __init__(self):
    #     list.__init__(self)
    #     self.strlen = 0
    
    # arcfile read 1st
    def read_arc(self,strfile, forfile=False, allformat = 0):
        '''read allstr.arc (and allfor.arc) file
            used in VASP.run()
            
        noted by JamesBourbon in 20220402    
        '''
        f= open(strfile,'r')
        currentStr = -1
        for line in f:
            if ('Energy' in line\
                or 'React' in line\
                or 'TS' in line\
                or 'SSW-Crystl' in line\
                or 'Str' in line):
                self.append(Str())
                currentStr +=  1
                self[currentStr].Lfor = False
                try:
                    self[currentStr].energy = float(line.split()[-1])
                    try :
                        self[currentStr].maxF = float(line.split()[-2])
                        if self[currentStr].maxF<0: self[currentStr].maxF=0
                    except :
                        self[currentStr].maxF = 990

                    if self[currentStr].energy.is_integer():
                        self[currentStr].energy = float(line.split()[-2])
                        self[currentStr].maxF = float(line.split()[-3])
                except:
                    self[currentStr].maxF = 990
                    if '***' not in line:
                        self[currentStr].energy = float(line.split()[-2])
                    else: 
                        self[currentStr].energy = 990
            elif 'Materials Studio' in line:
                currentStr +=  1
                self[currentStr].Lfor = False
                self[currentStr].energy = 0
                self.maxF = 0
            elif ('CORE' in line) or ('XXXX' in line):
                self[currentStr].addatom(line, 1 )
            elif ('PBC' in line ) and ('ON' not in line):
                self[currentStr].Latt= [float(x) for x in line.split()[1:]]
        f.close()
        if forfile:
            f = open(forfile,'r')
            currentStr= -1
            for line in f:
                if "For" in line:
                    self[currentStr].Lfor = True
                    currentStr += 1
                    iatom = 0
                    for atom in self[currentStr].atom: atom.force = [0.0, 0.0, 0.0]
                elif len(line.split()) == 6:
                    self[currentStr].add_stress(line)
                elif len(line.split()) == 3:
                    if "****" not in line: self[currentStr].add_force(line, iatom)
                    else:                  self[currentStr].add_force('0.0 0.0 0.0', iatom)
                    iatom += 1
            f.close()
        # str normalization
        for str in self:
            str: Str
            str.sort_atom_by_element()
            str.set_strinfo_from_atom()
            str.set_coord()
            str.Cell = str.Latt2Cell()
            if str.Lfor:
                str.set_for_list()
        if allformat:
            for str in self:
                str: Str
                # str.TransferToXYZcoordStr() # not use

    # arcfile read 2nd
    def build_coord_set(self, onestr=[0, 0], filename='allstr.arc'):
        '''bulid coordinate to Str.atom from allstr.arc file
        
        update in 20220902 by JamesBourbon      
        '''
        print('Build_Coord_Set')
        ind_getstr = -1  # str found, can be range in 'onestr'
        index_str = -1  # str added-in
        for line in open(filename):
            item = line.split()[:]  # line: list
            # head line
            if ('Energy' in item
                or 'SSW-fixlat' in item
                or 'SSW-Crystal' in item
                or 'Str' in item
                    or 'React' in item):
                ind_getstr += 1
                if onestr[1] == 0 or ind_getstr in range(onestr[0], onestr[1]):
                    index_str += 1
                    append_str = Str()
                    if 'from' not in item:
                        try:
                            append_str.energy = float(line.split()[-1])
                            mk = -2
                            if append_str.energy.is_integer():
                                append_str.energy = float(
                                    line.split()[-2])
                                mk = -3
                            try:
                                append_str.maxF = float(line.split()[mk])
                            except:
                                append_str.maxF = 990
                        except:
                            try:
                                append_str.energy = float(
                                    line.split()[-2])
                                append_str.maxF = float(line.split()[-3])
                            except:
                                append_str.energy = 1
                                append_str.maxF = 990
                    else:
                        append_str.energy = float(line.split()[-4])
                    self.strlen = index_str + 1
                    append_str.Lfor = False
                    ele_order = 1
                    self.append(append_str)
            elif 'Materials Studio' in line:
                # read from MS-gen-ed CAR file
                ind_getstr += 1
                if onestr[1] == 0 or ind_getstr in range(onestr[0], onestr[1]):
                    index_str += 1
                    append_str = Str()
                    append_str.energy = 0
                    append_str.maxF = 0
                    self.strlen = index_str + 1
                    append_str.Lfor = False
                    ele_order = 1
                    self.append(append_str)
            # Latt line
            if 'PBC' in line and 'ON' not in line \
                    and (onestr[1] == 0 or ind_getstr in range(onestr[0], onestr[1])):
                try:
                    self[index_str].Latt = [float(x)
                                            for x in line.split()[1:7]]
                except:
                    self[index_str].Latt = [99, 99, 99, 90, 90, 90]
            # atom line: add coord to atom, not to self.Coord
            if  (('CORE' in line) or ('XXXX' in line)) \
                and (onestr[1] == 0 or ind_getstr in range(onestr[0], onestr[1])):
                atom_ind = self[index_str].natom
                try:
                    self[index_str].addatom(line)  # add info to atoms
                except:
                    self[index_str].addatom('NaN 0 0 0')
                self[index_str].natom += 1  # update natom and atom_ind

        
    def build_for_set(self, onestr=[0, 0], filename='allfor.arc'):
        '''bulid force set from allfor.arc file
            if allfor.arc not exist, return None
        
        update by JamesBourbon in 20220902'''

        print('Build_For_Set')
        index_str = -1
        ind_getstr = -1
        try:
            f = open(filename, 'r')
        except:
            print('--No force info--')
            return None
        line = f.readline()
        if line == None:
            print('--No force info--')
            return None
        while line:
            item = line.split()[:]
            if ('For' in item):
                #print index,self[index+1].energy
                if onestr[1] == 0 or ind_getstr+1 in range(onestr[0], onestr[1]):
                    if 'from' not in item:
                        err = float(line.split()[-1]) - \
                            self[index_str+1].energy
                    else:
                        err = float(line.split()[-4]) - \
                            self[index_str+1].energy
                #print float(item.split()[-1]), self[index+1].energy, e
                    if (err > 1e-5):
                        print(index_str+1, float(line.split()
                                                [-1]), self[index_str+1].energy, err)
                        raise FooError('Error: incompatible Str with For')
                    else:
                        ind_getstr += 1
                        index_str += 1
                        self[index_str].Lfor = True
                        stress_line = f.readline()
                        if '**' not in stress_line:
                            try:
                                self[index_str].add_stress(stress_line)
                            except:
                                self[index_str].stress = [99999, 99999,
                                                        99999, 99999, 99999, 99999]
                        else:
                            self[index_str].stress = [99999, 99999,
                                                    99999, 99999, 99999, 99999]
                        self[index_str].maxF = 0.0
                        for i in range(self[index_str].natom):
                            for_line = f.readline()
                            if '**' not in for_line:
                                try:
                                    self[index_str].add_force(for_line, i)
                                    # self[index_str].For.append(list(float(x) for x in strfor.split()[:3]))
                                except:
                                    self[index_str].add_force(
                                        "999999 999999 999999", i)
                                    # self[index_str].For.append([999999, 999999, 99999])
                            else:
                                self[index_str].add_force(
                                    "999999 999999 999999", i)
                                # self[index_str].For.append([999999, 999999, 99999])
                            self[index_str].maxF = max(self[index_str].maxF, max(
                                map(abs, [x for x in self[index_str].atom[i].force])))
            line = f.readline()


    def arcinit(self, ind=[0, 0], strfile='allstr.arc', forfile='allfor.arc', allformat=0):
        '''init arcfile" read allstr and allfor(if exist)'''
        self.build_coord_set(ind, strfile)  # build coordinate set from allstr.arc
        self.build_for_set(ind, forfile)  # build force set from allfor.arc
        for str in self:
            str: Str
            str.sort_atom_by_element()  # need to sort-arange
            str.set_strinfo_from_atom()
            str.set_coord()
            str.Cell = str.Latt2Cell()
            if str.Lfor:
                str.set_for_list()
        if allformat:
            for str in self:
                try:
                    str.TransferToKplstr()  # not need
                except:
                    # print(str.energy)# always turn to this
                    continue



    # train-data read
    def train_data_init(self, strfile, forcefile= False):
        '''for read train_format data to AllStr object'''
        print("read from TrainStr (and TrainFor)")
        list.__init__(self)
        self.read_str(strfile, forcefile)

    def read_str(self, strfile, forcefile = False):
        '''read TrainStr.txt and TrainFor.txt'''
        currentStr = -1
        for line in open(strfile):
            if 'Start' in line:
                currentStr = currentStr+1
                self.append(Str())
                #self[currentStr].serial_num = currentStr
            elif 'Energy' in line:
                self[currentStr].energy= float(line.split()[2])
            elif 'lat' in line:
                self[currentStr].Cell.append([
                        float(x) for x in line.split()[1:4]])
            elif 'ele' in line and 'element' not in line:
                self[currentStr].addatom(line, 2 )

        self.numstr = currentStr + 1
        for i in range(self.numstr):
            # self[i].sort_atom_by_element() # need to read force before
            self[i].set_strinfo_from_atom()
            self[i].Latt = self[i].Cell2Latt()
            self[i].set_coord()
            self[i].Lfor = False
        # print(self[i].natom)
        print("TrainStr has read!")        
        
        if forcefile:
            f = open(forcefile,'r')
            currentStr = -1
            for line in f:
                if 'Start' in line:
                    currentStr = currentStr+1
                elif 'stress' in line:
                    self[currentStr].stress= [float(x) for x in line.split()[1:]]
                    iatom = 0
                elif 'force' in line:
                    self[currentStr].add_force(line, iatom ,2)
                    iatom = iatom +1
            for i in range(self.numstr):
                self[i].get_max_force()
                self[i].Lfor = True
            print("TrainFor has read!")
            
        for i in range(self.numstr):
            self[i].sort_atom_by_element() # need to read force before

    
    # print file
    def print_str_all(self,outfile='outstr.arc'):
        '''print All-Str by str.atom to arcfile-format, output is outfile
        
        substitute for gen_arc'''
        f = open(outfile,'w')
        f.write('!BIOSYM archive 2\n')
        f.write('PBC=ON\n')
        istr = 0

        for str in self:
            istr = istr+1

            f.write('     Energy     %8d    %8d     %12.6f\n'%(istr,istr,str.energy))
            f.write('!DATE\n')
            f.write('PBC %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f\n'%
                    (str.Latt[0],str.Latt[1],str.Latt[2],str.Latt[3],str.Latt[4],str.Latt[5]))
            #f.write('PBC %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f\n'%(100,100,100,90,90,90))
            for i,atom in enumerate(str.atom):
                f.write('%-3s %15.9f %15.9f %15.9f CORE %4d %2s %2s %8.4f %4d\n' %
                        (atom.ele_symbol,atom.xyz[0],atom.xyz[1],atom.xyz[2],i+1,atom.ele_symbol,atom.ele_symbol,atom.charge,i+1))
            f.write('end\nend\n')


    def print_str(self, printlist, outfile):
        '''print selected Str by atom to outfile in arc_format
                select in printlist
        
        same function as gen_arc, but outfile can be designated
        '''
        f = open(outfile,'w')
        f.write('!BIOSYM archive 2\n')
        f.write('PBC=ON\n')
        for istr in printlist:
            str=self[istr]

            f.write('     Energy     %8d    %8d     %12.6f\n'%(istr,istr,str.energy))
            f.write('!DATE \n')
            f.write('PBC %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f\n'%
                    (str.Latt[0],str.Latt[1],str.Latt[2],str.Latt[3],str.Latt[4],str.Latt[5]))
            #f.write('PBC %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f\n'%(100,100,100,90,90,90))
            for i,atom in enumerate(str.atom):
                f.write('%-3s %15.9f %15.9f %15.9f CORE %4d %2s %2s %8.4f %4d\n' %
                        (atom.ele_symbol,atom.xyz[0],atom.xyz[1],atom.xyz[2],i+1,atom.ele_symbol,atom.ele_symbol,atom.charge,i+1))
            f.write('end\nend\n')
        return

    def print_for(self,printlist,outfile='outfor.arc', flag= False):
        '''print outfor.arc from AllStr(Str) by atom'''
        fout = open(outfile, "w")
        for istr in printlist:
            str = self[istr]
            if flag:
                fout.write(" For   %d  %d  SS-fixlat   %12.6f\n"%(str.nminimum, str.numstr, str.energy))
            else:
                fout.write(" For   0  0  SS-fixlat   %12.6f\n"%(str.energy))

            fout.write("%15.8f %15.8f %15.8f %15.8f %15.8f %15.8f\n"%(
                str.stress[0], str.stress[1], str.stress[2],
                str.stress[3], str.stress[4], str.stress[5]))
            for atom in str.atom:
                fout.write("%15.8f %15.8f %15.8f\n"%(atom.force[0], atom.force[1], atom.force[2]))
            fout.write("\n ")
        fout.close()



    def gen_data_str(self, printlist, fname):
        '''generate str data by str.atom for NNtrain
        
        used in VASP-DFT arc2train
        
        '''
        with open(fname,'w') as fout:
            for istr in printlist:
                str = self[istr]
                #str.sortAtom()
                fout.write(" Start one structure\n")
                fout.write("Energy is %12.6f eV\n"%(str.energy))
                fout.write("total number of element is %5d\n"%str.natom)
                fout.write("element in structure:\n")
                fout.write("symbol %s\n" %reduce(lambda a,b:a+b , ["%4s"%s   for s   in str.ele_nameList]))
                fout.write("No.    %s\n" %reduce(lambda a,b:a+b , ["%4d"%num for num in str.eleList]))
                fout.write("number %s\n" %reduce(lambda a,b:a+b , ["%4d"%num for num in str.natompe]  ))
                for lat in str.Cell:
                    fout.write("lat %15.8f  %15.8f  %15.8f\n"%(lat[0], lat[1], lat[2]))
                for atom in str.atom:
                    fout.write("ele %4s %15.8f  %15.8f  %15.8f  %15.8f\n"%(
                            atom.ele_num, atom.xyz[0], atom.xyz[1], atom.xyz[2], atom.charge))
                fout.write(" End one structure\n\n")


    def gen_data_for(self,printlist ,fname):
        '''generate force data by atom for NNtrain
        
        used in VASP-DFT arc2train
        '''
        with open(fname,'w') as fout:
            for istr in printlist:
                str = self[istr]
                fout.write(" Start one structure\n")
                fout.write("stress %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f\n"%(
                    str.stress[0], str.stress[1], str.stress[2],
                    str.stress[3], str.stress[4], str.stress[5]))
                for atom in str.atom:
                    fout.write("force %4d %15.8f %15.8f %15.8f\n"%(
                        atom.ele_num, atom.force[0], atom.force[1], atom.force[2]))
                fout.write(" End one structure\n\n")
                

    # problem will from not use Str.atom
    def gen_arc(self, set, filename='outstr.arc',):   # ,fmode=0, dirname='test', HydraPre=False):
        """ generate new arc file, set   => list (ot iter) of structure want to output
                    filename => name of output arc file
                        fmode => 0 is defaule format, 1 is allstr.arc type
                            press => default is zero, else will add pressure effect
                
            use it for get new_arc_file from AllStr(Str) (not from Str.atom)
            
            noted by JamesBourbon in 20220902
        """
        with open(filename, 'w') as fout:
            fout.write("!BIOSYM archive 2\nPBC=ON\n")
            for i in set:
                energy = self[i].energy
# PV adding but not use
#               if HydraPre:     energy = self[i].energy + self[i].Volume * HydraPre / 160.2176487
#               else:         energy = self[i].energy
                if not self[i].Lfor:
                    fout.write(
                        "\t\t\t\tEnergy\t%8d        -0.0000  %17.6f\n" % (i+1, energy))
                else:
                    fout.write("\t\t\t\tEnergy\t%8d        %10.6f  %17.6f\n" % (
                        i+1, self[i].maxF, energy))

# arc_file mode to give different head (not used)
#                if fmode==0 :  fout.write("\t\t\t\tEnergy\t%8d        %10.6f  %17.6f\n"%(i+1,self[i].maxF,energy))
#                elif fmode==2:
#                fout.write("\t\t\t\tEnergy\t%8d        %8.4f %17.6f\n"%(i, self[i].q2, energy))
#                elif fmode==1:
#                    fout.write("  Str\t%8d  %8d  SSW-fixlat   %12.6f\n"%(self[i].numsym,i,self[i].energy))
#                elif fmode==8:
#                    fout.write("  Str\t%8d  %8d  SSW-fixlat   %12.6f   from %s\n"%(self[i].numMinimum,self[i].numStr, self[i].energy, dirname))
                fout.write("!DATE\n")

                lat = self[i].Latt
                if not lat:
                    lat = self[i].Cell2Latt()
                fout.write("PBC  %12.6f%12.6f%12.6f%12.6f%12.6f%12.6f\n" %
                        (lat[0], lat[1], lat[2], lat[3], lat[4], lat[5]))
                try:
                    # default: read from Str.atom
                    for i,atom in enumerate(self[i].atom):
                        fout.write('%-3s %15.9f %15.9f %15.9f CORE %4d %2s %2s %8.4f %4d\n' %
                        (atom.ele_symbol, atom.xyz[0], atom.xyz[1], atom.xyz[2],
                            i+1, atom.ele_symbol, atom.ele_symbol, atom.charge, i+1))
                    # print("gen_arc from Str.atom")
                except:
                    # if order-ed by Str.sort_by_ele
                    print("Wronging: gen_arc from Coord but not Str.atom")
                    for j in range(self[i].natom):
                        ele_name = self[i].Ele_Name[j]
                        coord = self[i].Coord[j]
                        fout.write("%-2s%18.9f%15.9f%15.9f CORE %4d %-2s %-2s   0.0000 %4d\n" %
                                (ele_name, coord[0], coord[1], coord[2], j+1, ele_name, ele_name, j+1))
                fout.write("end\nend\n")

    def gen_forarc(self, iter, fname='outfor.arc'):
        """ generate new For arc file
        
        use it for get new_arc_file from AllStr(Str) (not from Str.atom)
            
        noted by JamesBourbon in 20220902    
        """
        with open(fname, 'w') as fout:
            for i in iter:
                fout.write("For %4d  %12.6f\n" % (i+1, self[i].energy))
                stres = self[i].stress[0:6]
                fout.write("   %15.9f %15.9f %15.9f %15.9f %15.9f %15.9f\n"
                            % (stres[0], stres[1], stres[2], stres[3], stres[4], stres[5]))
                try:
                    for atom in str.atom:
                        fout.write("   %15.9f %15.9f %15.9f\n"%(atom.force[0], atom.force[1], atom.force[2]))
                except:
                    for j in range(self[i].natom):
                        for_axis = self[i].For[j]
                        fout.write("   %15.9f %15.9f %15.9f \n" %(for_axis[0], for_axis[1], for_axis[2]))
                fout.write("\n")
                

    def gen_POSCAR_VASP(self, num):
        """ generate POSCAR, num is the num-th structure, fname is the name of output POSCAR
        
        refine by JamesBourbon by using str.outPOSCAR
        """
        fname = 'POSCAR_'+str(num)
        #fout = open(fname, 'w')
        self[num].outPOSCAR(fname)


    def gen_REANN_configuration(self, iter, fname="configuration"):
        """ generate configuration file for REANN

        Inputs:
            iter: iterable to indicate range which print
            fname: output filename, default "configuration"
        """
        print("== Generating Congifuration file for REANN ==")
        point = 0
        with open(fname, 'w') as fout:
            for i in iter:
                point += 1
                str = self[i]
                str: Str
                fout.write(f"point={point:>4}\n")
                for lat in str.Cell:
                    fout.write("%12.8f  %12.8f  %12.8f\n"%(lat[0], lat[1], lat[2]))
                fout.write("pbc 1  1  1\n") # for pbc in reann
                for atom in str.atom:
                    fout.write(f"{atom.ele_symbol:4s} {atom.mess:6.3f} {atom.xyz[0]:12.8f} {atom.xyz[1]:12.8f} {atom.xyz[2]:12.8f} {atom.force[0]:10.6f} {atom.force[1]:10.6f} {atom.force[2]:10.6f} \n")
                fout.write(f"abprop: {str.energy}\n")
        print(f"== Done with configuration file containing {point} structures ! ==")



    def gen_INPUT_SIESTA(self, i, dir='./'):
        """ generate INPUT_DEBUG, i is the i-th structure, dir is the path of INPUT_DEBUG.template"""
        fname = 'INPUT_DEBUG_'+str(i)
        dict = {}
        dict['NumberOfAtoms'] = str(self[i-1].natom)
        dict['NumberOfSpecies'] = str(len(self[i-1].sp))
        dict['cell'] = self[i-1].printCell()
        dict['coord'] = self[i-1].printCoord()
        #print dict
        with open('%s/INPUT_DEBUG.template' % (dir)) as ftemp:
            _data = string.Template(ftemp.read())
        with open(fname, 'w') as fout:
            fout.write(_data.safe_substitute(dict))
        #os.system('cp %s/*psf .'%(dir))

    def GenGIN_GULP(self, num):
        """ generate gin file, num is the num-th structure, fname is the name of output gin file"""
        fname = 'test.gin_'+str(num)
        with open(fname, 'w') as fout:
            fout.write("conp md\n")
            fout.write("vector\n")
            for x in self[num-1].Cell:
                fout.write("   %15.8f  %15.8f  %15.8f\n" % (x[0], x[1], x[2]))
            fout.write("cart\n")

            for i in range(self[num-1].natom):
                xa = self[num-1].Coord[i]
                ele = self[num-1].Ele_Name[i]
                fout.write("   %2s CORE  %12.6f  %12.6f  %12.6f\n" %
                            (ele, xa[0], xa[1], xa[2]))
            #fout.write(Tp.gulp_TiO2.safe_substitute())


    def read_coord_set_from_POSCAR(self, filename='POSCAR'):
        '''build coord set from POSCAR by using Str.build_coord_set_from_POSCAR()
        
        re-defined by JamesBourbon
        '''
        print('Build_Coord_Set_from_POSCAR')
        read_str = Str()
        read_str.build_coord_set_from_POSCAR(filename)
        self.append(read_str)


    def gen_coordination_patterns(self):
        '''generate coordination patterns set from all Str
        
        [Returns] set: coordination-patterns format set
        '''
        allstr_coordination_pattern = set()
        for struc in self:
            struc: Str
            allstr_coordination_pattern.update(struc.coordination_pattern())
        return allstr_coordination_pattern


    def get_all_element(self):
        '''get all element in all_str'''
        elements = set()
        for struc in self:
            elements.update(tuple(struc.sp.keys()))
        return elements


    def shuffle(self, shuffletime):
        """ shuffle structures"""
        for i in range(shuffletime): random.shuffle(self)


    def screen_upper(self):
        for str in self:
            str.screenuppersurf()


    def single_split(self):
        '''split allstr to each input dir and run SSW-NN'''
        workdirs=[]
        for i in range(0,len(self)):
            workdir='para%d'%(i)
            workdirs.append(workdir)
            os.mkdir(workdir)
            self.print_str([i],'%s/input.arc'%(workdir))
            os.system('cp ./input %s'%(workdir))
            os.system('ln -s %s/*.pot %s'%(os.getcwd(),workdir))
        return workdirs


    def out_mol_file(self, printlist=False,outfile='out.mol'):
        if not printlist: printlist = range(len(self))

        f = open(outfile,'w')
        for istr in printlist:

            str=self[istr]
            str.GetAllBond()
            f.write('%d  %d\n'%(istr,istr))
            f.write('     RDKit \n')# %15.9f\n'%(str.energy))
            f.write('\n')
            f.write('%3d%3d  0  0  1  0  0  0  0  0999 V2000\n'%(str.natom,str.nbond))
            for i,atom in enumerate(str.atom):
                f.write('%10.4f%10.4f%10.4f%3s 0  0  0  0  0  0  0  0  0%3d  0  0\n' %
                         (atom.xyz[0],atom.xyz[1],atom.xyz[2],atom.ele_symbol,i+1))
            for bond in str.allbond:
                f.write('%3d%3d%3d  0\n'%(bond[0]+1,bond[1]+1,bond[2]))
            f.write('M  END\n')
            f.write('$$$$\n')

        f.close()


    def para_run(self,f,num=4):
        '''parallel wrapper to run each structure set with f function
        
        Returns: [list] of all results from each Str
        '''
        pool = Pool(processes=num)
        a=[]
        a.append(pool.map_async(f, self))
        pool.close()
        pool.join()
        return a


    def random_arange(self, shuffletime=50):
        '''random arange allstr
        
        Args: shuffletime (defalut 50)
        '''
        for i in range(shuffletime): random.shuffle(self)


    def para_shortest_bond(self,Ncore=4):
        d_results = self.para_run(ParaWrap_Shortestbond,Ncore)
        record=0; d = {}; strlist ={};strlist2=[]
        for x in d_results:
            for y in x.get():
                record +=1
                self[record-1].d={}
                #print y
                shortd = min(y.values())
                for mm in range(len(y)):
#                   print record,y.keys()[mm],y.values()[mm],self[record-1].energy,shortd
                    if shortd == y.values()[mm] : shortb = y.keys()[mm]
                tt = str(shortd)+str(shortb)+str(self[record-1].energy)
#               print record, tt
#               if shortd < 0.3 and shortb =='H-H': continue
#               if shortd < 0.5 and shortb =='H-C': continue
#               if shortd < 0.5 and shortb =='H-O': continue
#               if shortd < 0.5 and shortb =='H-N': continue
#               if shortd < 0.8 and shortb =='O-O': continue
#               if shortd < 1.0 and shortb =='O-Rh': continue
                if tt not in strlist :
#                   if ( abs(self[strlist[shortd]].energy-self[record-1].energy) > 0.0001  or
#                        strlist['tt'+str(strlist[shortd])] != shortb ) :
#                       strlist2.append(record-1)
#                   else:
#                       print 'same distance but different energy--bond'
#                       print 'XXXX',shortd, strlist[shortd],record-1
#               else:
#                   strlist[str(shortd)+str(shortb)+(self[record-1].energy)]= record-1
#                   strlist['tt'+str(record-1)]= shortb
                    strlist[tt]= record-1
                    strlist2.append(record-1)
#               strlist[y.values()[0]]=record-1
                #print strlist
                #if record==1:
                #    for bondtype, bonddis in y.items():
                #        self[record-1].d[bondtype]=bonddis
                #        d[bondtype]=bonddis
                #else:
                for bondtype, bonddis in y.items():
                    self[record-1].d[bondtype]=bonddis
                    try :
                        d[bondtype]=min(bonddis,d[bondtype])
                    except:
                        d[bondtype] = bonddis
#       print 'remove duplicate'
#       for i in strlist.keys():
#           print i,strlist[i],self[strlist[i]].energy
# serial version
        #d = AllStr[0].Shortest_bond()
        #for i in range(1,len(self)):
        #    #print i,AllStr[i].bonds
        #    for x, y in AllStr[i].Shortest_bond().items():
        #        d[x]=min(y,d[x])
        return d, strlist2

    def filter(self,b):
        return AllStr(x for x in self if x.myfilter(b))

    def filter_byset(self,b):
        return AllStr(x for i,x in enumerate(self) if i in b)

    def filter_byd(self,d1,d2,type=1):
        return AllStr(x for x in self if x.myfilter_byd(d1,d2,type))

    def filter_byshortd(self,d1):
        return AllStr(x for x in self if x.myfilter_byshortd(d1))

    def sort_by_energy(self):
        return AllStr(sorted(self, key=lambda x: x.energy))
    
    def sort_by_stress(self):
        return AllStr(sorted(self, key=lambda x: max(x.stress)))
    
    def sort_by_element_natom(self,ele_name):
        return AllStr(sorted(self, key=lambda x: x.sp.get(ele_name, 0)))

    def filter_byQE(self,Etol,Qtol,L=False):
        b=range(len(self))
        for i in range(len(self)):
            if b[i] < 0 : continue
            self._filter_byQE(i,b,Etol,Qtol,L)
        print(b)
        return AllStr(x for i,x in enumerate(self) if b[i]>-1)

    def _filter_byQE(self,i,b,Etol,Qtol,L=False):
        for j in range(len(self)):
            if i>=j or b[j]<0: continue
            if (abs(self[i].energy -self[j].energy) < Etol
                and ( (L and abs(self[i].Q[3]-self[j].Q[3]) < Qtol*self[i].Q[3]
                and  abs(self[i].Q[4]-self[j].Q[4]) < Qtol*self[i].Q[4]
                and  abs(self[i].Q[5]-self[j].Q[5]) < Qtol*self[i].Q[5] ) or
                ( not L and abs(self[i].Q[0]-self[j].Q[0]) < Qtol*self[i].Q[0]
                and  abs(self[i].Q[1]-self[j].Q[1]) < Qtol*self[i].Q[1]
                and  abs(self[i].Q[2]-self[j].Q[2]) < Qtol*self[i].Q[2] )) ) : b[j]=-1

    def sort_by_force(self):
        '''sort Strs by force, return new allstr'''
        return AllStr(sorted(self, key=lambda x: x.maxF))
    
    def sort_by_natom(self):
        """sort structure by natom, return new allstr"""
        return AllStr(sorted(self, key=lambda x: x.natom))


    
        

    '''
    def get_all_smi_name(self,numproc=24,flag =2):
        if flag == 1:
            self.cal_all_bond_matrix(numproc=numproc)
        if flag == 2:
            self.cal_all_fake_bond_maxtrix(numproc=numproc)
        self.cal_all_seg_molecular (numproc=numproc)

        allgroup = []
        for str in self:
            substr = [[] for i in np.unique(str.group)]
            for id,atom in enumerate(str.atom):
                atom.id = id
                substr[str.group[id]-1].append(atom)
            substr = sorted(substr, key=lambda x:calmass(x), reverse=True)
            if flag == 1: allgroup.append((substr, str.bmx2D, str.Cell,[],1,[]))
            if flag == 2: allgroup.append((substr, str.bmx2D, str.Cell,str.bondneed,2,str.surfaceatom))
            #print str.bondneed

        pool = Pool(processes=numproc)
        result = pool.map_async(calAllName, allgroup)
        pool.close(); pool.join()

        for istr,(str,re) in enumerate(zip(self,result.get())):
            str.allmol = re
            str.id     = istr

        for str in self:
            str.sminame, strflag = glueSegStr(str.allmol)
    

    def get_all_pure_smi_name(self,numproc=24,flag =2):
        #self.calAllSegmolecular (numproc=numproc)
        if flag == 1:
            self.cal_all_bond_matrix(numproc=numproc)
        if flag == 2:
            self.cal_all_fake_bond_maxtrix(numproc=numproc)
        self.cal_all_seg_molecular (numproc=numproc)
        allgroup = []
        for str in self:
            substr = [[] for i in np.unique(str.group)]
            for id,atom in enumerate(str.atom):
                atom.id = id
                substr[str.group[id]-1].append(atom)
            substr = sorted(substr, key=lambda x:calmass(x), reverse=True)
            if flag == 1: allgroup.append((substr, str.bmx2D, str.Cell,[],1,[]))
            if flag == 2: allgroup.append((substr, str.bmx2D, str.Cell,str.bondneed,2,str.surfaceatom))
            #print str.bondneed

            #substr = [[]]
            #for id,atom in enumerate(str.atom):
            #    atom.id = id
            #    substr[0].append(atom)
            #if flag == 1: allgroup.append((substr, str.bmx2D, str.Cell,[],1,[]))
            #if flag == 2: allgroup.append((substr, str.bmx2D, str.Cell,str.bondneed,2,str.surfaceatom))

        pool = Pool(processes=numproc)
        result = pool.map_async(calAllpureName, allgroup)
        pool.close(); pool.join()

        for istr,(str,re) in enumerate(zip(self,result.get())):
            str.allmol = re
            str.id     = istr

        for str in self:
            str.puresminame, strflag = glueSegStr_pure(str.allmol)

    '''

# def ParaWrap_CoorPatterns(x):
#     return x.gen_coordination_patterns()


if __name__ == '__main__' :
    test= AllStr()
    test.arcinit([0,0],'input.arc')
    test.print_str([0],'test.arc')
    #test.read_arc('input.arc',allformat = 1)
    #test.GenGIN_GULP(0)
#    test= change()
#    test.read_arc('input.arc')
#    for str in test:
#        for iatom in range(str.natom):
#            if str.atom[iatom].ele_num == 6:
#                neighO = 0
#                atomneib =str.specialneighbour(iatom,8)
#                print atomneib
#                for dis in atomneib.values():
#                    if dis < 2.0:
#                        neighO= neighO +1
#                print neighO
#    test[0].outPOSCAR('fe')
