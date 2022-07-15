# coordination pattern lib
# JamesBourbon update all 20220607
# NOT recommend to use decendence
# use set format to store all data
import json
import os
import time    
    
class CoordinationPatterns():
    """All Coordination Patterns Object list and Operations
    
    Example for init_patterns and CoordinationPatterns.patterns:
        {('O': (('Au', 2.0, 1), ('Pd', 2.0, 2))), 
        ('Pd': (('O', 2.0, 4), ('Au', 1.8, 1), ('Au', 1.9, 1))),
        ('Pd': (('O', 2.0, 3), ('O', 1.9: 1))), 
        ('Pd': (('O', 2.0, 5), ('O', 2.1, 1))),
        ('Au': (('O', 2.0, 3), ('O', 2.1, 1)))}
    type: a set contain name and all OneCoordinatonPattern tuple
    
    Args:
        init_patterns: a set contain all OneCoordiantionPattern.patterns
        name: Patterns name
        output: running print-out message file, default to None
        limit: update limit; <1 to rate; >1 to num; should >0
    """
    def __init__(self, init_patterns:set=set(), name="lattice", output="", limit=1):
        self.name = name
        self.patterns = set()
        self.output_file = output
        self.print_out = ""
        if limit<=0:
            print("Wrong update_limit setting !!!")
            raise ValueError
        self.update_limit = limit
        if bool(init_patterns):
            self.update_from_patterns(patterns=init_patterns)
        self.set_operation_log()
    
    def add_one_pattern(self, pattern:tuple):
        self.patterns.add(pattern) 
        
    def is_pattern_in(self, pattern:tuple):
        return pattern in self.patterns

    def is_patterns_allin(self, patterns:set):
        return patterns.issubset(self.patterns)
    
    def update_from_patterns(self, patterns: set):
        '''update CoordinationPatterns.patterns from set or another coordination patterns'''
        self.patterns.update(patterns)
    
    def update_patterns_from_structure(self, coor_patterns: set, name="a"):
        """update patterns from structure_coordination_set, set made by API, 
        program will judge the patterns should be update or not by limit
        
        Args:
            coor_patterns (list): list of OneCoordiantionPattern objects
            name (str): name of coor_patterns added. Defaults to "a".

        Returns:
            bool: CoordinationPatterns.patterns have updated or not
        """
        mes = f"read pattern from {name} coordination-pattern list\n"
        self.print_out += mes  
        # use set to get new patterns 
        new_patterns = coor_patterns.difference(self.patterns)
        # statistica 
        update_num = len(new_patterns)
        all_num = len(coor_patterns)
        update_rate = update_num / all_num
        limit = self.update_limit
        self.print_out += f"Num of Add Patterns is {all_num}\n"
        self.print_out += f"Num of New Patterns is {update_num}\n"
        self.print_out += f"Update rate is {update_rate:.4}\n"
        self.print_out += f"Update limit is {limit}\n"
        # limit can be >1 ref to new_num or <1 ref to update_rate
        if limit < 1:
            updating = (update_rate >= limit)
            self.print_out += "limit<1, judge by update_rate\n"
        else: 
            updating = (update_num >= limit)
            self.print_out += "limit>1, judge by updage_num\n"
        if updating:
            mes = "Coordination Patterns Updated!!\n"
            # update buffer to coordination_pattern
            self.patterns.update(new_patterns)
            # print message
            self.print_out += mes 
            return True
        else:
            mes = "rate BELOW limit! Not update Coordination Patterns\n "
            self.print_out += mes
            return False
        
    def bond_num(self):
        '''get number of coordination bonds in patterns'''
        count = 0
        for pattern in self.patterns:
            count += len(pattern[1])
        return count

    # print-out attribute
    def set_operation_log(self):
        '''set running log'''
        start_str = f"running coor-pattern module in {time.ctime()}\n"
        if bool(self.output_file):
            if os.path.exists(self.output_file):
                with open(self.output_file, "w") as fo:
                    fo.write(start_str)
        else:
            self.print_out += start_str
        
    
    def print_operation_log(self):
        '''print running log in self.print_out'''
        filename = self.output_file
        if bool(filename): 
            with open(filename, "a") as fo:
                fo.write(self.print_out)
                print(f"print_out write in {filename}")
        else:
            print(self.print_out)
        self.print = "" # reset print-out when print
        
    
    def coordinations_dict(self):
        """return the simple coordiantion pattern of all atom
        
        Returns:
            (dict): name : patterns_values
        """
        coordination_dict = {
            self.name: self.patterns
        }
        return coordination_dict
    
    def print_all_coordinations_simple(self):
        simple_str = f"{self.name}: \n"
        coordination_str_list = []
        for pattern in self.patterns:
            coordination_str_list.append(str(pattern))
        coordination_str = ",\n".join(coordination_str_list)
        return simple_str + coordination_str
        
        
    def print_all_coordinations_full(self):
        """give a json_format string for all coordination patterns
        
        Returns:
        (str): json format string for coordination
        """ 
        all_coordination_patterns = []
        # pattern example : (Au, ((O, 1.6, 2), (Au, 1.8. 1)))
        for pattern in self.patterns:
            central_atom: str = pattern[0]
            atom_coor_pattern: tuple = pattern[1]
            coordination_objs_list = []
            for atom_coor in atom_coor_pattern:
                coordination_obj = {
                    "coordination_atom": atom_coor[0],
                    "distance": float(atom_coor[1]),
                    "coordination_numbers": int(atom_coor[2])
                }
                coordination_objs_list.append(coordination_obj)
            coordination_pattern_obj = {
                "central_atom": central_atom,
                "coordination_pattern": coordination_objs_list
            }
            all_coordination_patterns.append(coordination_pattern_obj)
        all_coordination_patterns_json = {
            "patterns_name": self.name,
            "all_coordination_patterns": all_coordination_patterns
        }
        all_patterns_json_formstr = json.dumps(
                all_coordination_patterns_json, indent=4)
        return all_patterns_json_formstr
    
    def save_file(self, string: str, filename: str):
        '''save string to file'''
        with open(filename, 'w') as fo:
            fo.write(string)
        
    # coding in 2021-12-18, refine to set_format
    def read_coordination_json(self, filename):
        """set coordination pattern by existed json file
        """    
        with open(filename, 'r') as fo:
            json_string = fo.read()
            json_dict = json.loads(json_string)
        patterns_all = json_dict["all_coordination_patterns"]
        for coor_pattern_json in patterns_all:
            coor_list = []
            for coor_json in coor_pattern_json["coordination_pattern"]:
                coor_list.append((coor_json["coordination_atom"], 
                    coor_json["distance"], coor_json["coordination_numbers"]))
            coor_tuple = tuple(coor_list)
            one_coor_pattern = (coor_pattern_json["central_atom"], coor_tuple)
            self.patterns.add(one_coor_pattern)
        message = f"-------read existed coordination_patterns from {filename}--------\n\n"
        self.print_out += message 

if __name__ == "__main__":
    filename = 'test/Au8Pd8O24_patterns_full.json'
    arc_coor = CoordinationPatterns()
    arc_coor.read_coordination_json(filename)
    print(arc_coor.patterns)