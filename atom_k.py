# noted by JamesBourbon in 20220401
import PeriodicTable as PT

class S_atom(object):
    '''Single Atom object for anything of an atom
    
    noted by JamesBourbon in 20220401
    '''                                                        
    def __init__(self,coord,element_num):                                        
        self.xyz= coord                                                      
        self.ele_num = element_num  
        self.ele_symbol = PT.Eletable[self.ele_num-1]                                                 
        self.force = []                                                      
        self.symf = []                                                       
        self.species = 0                                                     
        self.bondlist = []                                                   
        self.bondtype = []                                                   
        self.Ctye = 'empty'                                                  
        self.charge = 0                                                                   
                            
