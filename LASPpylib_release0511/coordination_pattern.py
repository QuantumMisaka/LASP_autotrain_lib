import json

class OneCoordinationPattern():
    """Coordination Pattern Object by one central atom
        Central Coordinate in Cartesian can be specified
    
    Example for OneCoordinationPattern.pattern:
        ('Pd', {('O', 2.0): 2, ('O', 2.1): 2})
    Or
        ('Au', {('O', 2.0): 3, ('O', 2.1): 1})
    type: str-dict-pair in tuple
    
    Args:
        pattern_tuple (tuple): pattern tuple
        central_coordination (float-tuple): default (0.0, 0.0, 0.0)
    """
    def __init__(self, pattern_tuple: tuple, name="", 
                 central_coordinate=[0.0,0.0,0.0]):
        self.pattern = pattern_tuple
        self.coordinate = central_coordinate
        self.name = name
    
    def central_atom(self):
        """return the central atom in Coordination"""
        return self.pattern[0]
    
    def coordination(self):
        """return the simple coordiantion pattern of central atom"""
        return self.pattern
    
    def coordination_full(self):
        """return the full coordination pattern of central atom 
        in json format"""
        full_json = {}
        full_list = []
        for atom_coor, atom_num in self.pattern[1]:
            coor_obj = {
                "coordination_atom": atom_coor[0],
                "coordination_distance": float(atom_coor[1]),
                "coordination_number": int(atom_num)
            }
            full_list.append(coor_obj)
        full_json['center'] = self.pattern[0]
        full_json['center_coordinate'] = self.coordinate
        full_json["coordination_pattern"] = full_list
        full_json_formstr = json.dumps(full_json, indent=4)
        return full_json_formstr
    
    
class CoordinationPatterns():
    """All Coordination Patterns Object and Operations
    
    Example for init_patterns and CoordinationPatterns.patterns:
        [('O', {('Au', 2.0): 1, ('Pd', 2.0): 2}), 
        ('Pd', {('O', 2.0): 4}), 
        ('Pd', {('O', 2.0): 3, ('O', 1.9): 1}), 
        ('Pd', {('O', 2.0): 5, ('O', 2.1): 1}),
        ('Au', {('O', 2.0): 3, ('O', 2.1): 1})]
    type: a list contain name and all OneCoordinatonPattern.pattern tuple
    
    Args:
        init_patterns: [] or list contain all OneCoordiantionPattern
    """
    def __init__(self, init_patterns: list, name="lattice"):
        self.name = name
        self.patterns = []
        self.pattern_pack = []
        self.print_out = ""
        if bool(init_patterns):
            self.add_patterns_from_list(init_patterns)
          
    def is_pattern_replica(self, coor_pattern: OneCoordinationPattern):
        return (coor_pattern.pattern in self.patterns)
    
    def add_pattern(self, coor_pattern: OneCoordinationPattern):
        """Add NEW coor_pattern to CoordinationPatterns.patterns
        
        Args:
            coor_pattern (OneCoordinationPattern): 
        Returns:
            [bool]: [Judge this coor_pattern is NEW or not]
        """
        if coor_pattern.pattern not in self.patterns:
            self.patterns.append(coor_pattern.pattern)
            self.pattern_pack.append(coor_pattern)
            return True
        else:
            return False
    
    def add_patterns_from_list(self, coor_patterns: list, name="a"):
        """add patterns from list or another CoordinationPattern.patterns
        
        Args:
            coor_patterns (list): list of OneCoordiantionPattern objects
            name (str): name of coor_patterns added. Defaults to "a".

        Returns:
            bool: CoordinationPatterns.patterns have updated or not
        """
        patterns_replica_degree = [] # bool-list
        mes = f"read pattern from {name} coordination-pattern list\n"
        self.print_out += mes
        for coor_pattern in coor_patterns:
            one_pattern_degree = self.add_pattern(coor_pattern=coor_pattern)
            if one_pattern_degree:
                update_notice = f"New Pattern: {coor_pattern.pattern} updated!\n"
                self.print_out += update_notice
            if not one_pattern_degree:
                replica_notice = f"Pattern: {coor_pattern.pattern} already have!\n"
                self.print_out += replica_notice
            patterns_replica_degree.append(one_pattern_degree)
        if True in patterns_replica_degree:
            mes = "Coordination Patterns Updated!\n"
            self.print_out += mes
            # bool can be add like 0/1
            new_patterns_num = sum(patterns_replica_degree)
            patterns_sum_num = len(patterns_replica_degree)
            replica_num = patterns_sum_num - new_patterns_num
            mes2 = f"Update {new_patterns_num} NEW coordination patterns!\n"
            mes3 = f"Find {replica_num} Replica patterns!\n"
            self.print_out += mes2
            self.print_out += mes3
            return True
        else:
            mes = "NO update! No new Coordination Patterns\n "
            self.print_out += mes
            return False
    
    
    def print_operation_log(self, filename=""):
        if bool(filename):
            with open(filename, "w") as fo:
                fo.write(self.print_out)
                print(f"print_out write in {filename}")
        else:
            print(self.print_out)
        self.print = "" # reset print-out when print
        
    
    def coordinations_out(self):
        """return the simple coordiantion pattern of all atom
        
        Returns:
            (dict): name : patterns_values
        """
        coordination_dict = {
            self.name: self.patterns
        }
        return coordination_dict
    
    def all_coordinations_simple(self):
        simple_str = f"{self.name}: \n"
        coordination_str_list = []
        for pattern in self.patterns:
            coordination_str_list.append(str(pattern))
        coordination_str = ",\n".join(coordination_str_list)
        return simple_str + coordination_str
        
        
    def all_coordination_full(self):
        """give a json_format string for all coordination patterns
        
        Returns:
        (str): json format string for coordination
        """ 
        all_coordination_patterns = []
        for pattern in self.patterns:
            central_atom = pattern[0]
            atom_coor_pattern: dict = pattern[1]
            coordination_objs_list = []
            for atom_coor, atom_num in atom_coor_pattern.items():
                coordination_obj = {
                    "coordination_atom": atom_coor[0],
                    "distance": float(atom_coor[1]),
                    "coordination_numbers": int(atom_num)
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
        
        
    # coding in 2021-12-18
    def coordination_json_read(self, filename):
        """set coordination pattern by existed json file
        """    
        with open(filename, 'r') as fo:
            json_string = fo.read()
            json_dict = json.loads(json_string)
        patterns_all = json_dict["all_coordination_patterns"]
        for coor_pattern_json in patterns_all:
            coor_dict = {}
            for coor_json in coor_pattern_json["coordination_pattern"]:
                coor_dict[(coor_json["coordination_atom"], 
                    coor_json["distance"])] = coor_json["coordination_numbers"]
            one_coor_pattern = (coor_pattern_json["central_atom"], coor_dict)
            self.patterns.append(one_coor_pattern)
        message = f"-------read existed coordination_patterns from {filename}--------\n\n"
        self.print_out += message 
