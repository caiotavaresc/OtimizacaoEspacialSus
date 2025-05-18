"""
Here we reproduce some repairing operators
"""
from collections import Counter
import random

def can_remove_unit_from_region(unit_index, region):
    return True

def evaluate_region_neighbours(unit_index, region, genotype, graph, max_units):
    neigh_regions = []
    
    unit = list(graph.keys())[unit_index]
    for i in graph[unit]["vizinhos"]:
        reg_viz = genotype[graph[i]["seq"]]
        if reg_viz != region and reg_viz not in neigh_regions:
            neigh_regions.append(reg_viz)
            
    if len(neigh_regions) == 0:
        return None
            
    freq_dict = get_ordered_freq_dict(genotype)
    freq_dict_keys = list(freq_dict.keys())
    for i in freq_dict_keys:
        if i in neigh_regions and freq_dict[i] < max_units:
            return i
        
    return None

def find_indexes(genotype, region):
    return [i for i, valor in enumerate(genotype) if valor == region]

def get_freq_dict(genotype):
    frequency = Counter(genotype)
    freq_dict = dict(frequency)
    return freq_dict

def get_ordered_freq_dict(genotype):
    freq_dict = get_freq_dict(genotype)
    ordered = dict(sorted(freq_dict.items(), key=lambda item: item[1]))
    return ordered

def get_ordered_freq_dict_desc(genotype):
    freq_dict = get_freq_dict(genotype)
    ordered = dict(sorted(freq_dict.items(), key=lambda item: item[1], reverse=True))
    return ordered

def get_a_region_with_less_units_then_minimum(genotype, min_size):
    freq_dict = get_freq_dict(genotype)
    
    for i in freq_dict:
        if freq_dict[i] < min_size:
            return i
    
    return None

def change_all_index(from_ind, to_ind, genotype):
    for i in range(len(genotype)):
        if genotype[i] == from_ind:
            genotype[i] = to_ind

def transform_region_indexes(genotype, num_of_regions):
    freq_dict = get_freq_dict(genotype)
    
    keys = list(freq_dict.keys())
    indexes = range(num_of_regions)
    
    non_used_indexes = []
    for i in indexes:
        if i not in keys:
            non_used_indexes.append(i)
    
    for i in keys:
        if i not in indexes:
            new_index = non_used_indexes.pop()
            change_all_index(i, new_index, genotype)

def merge_region_in_neighbours(genotype, reg_i, graph):
    units_list = list(graph.keys())
    
    while reg_i in genotype:
        for k in range(len(genotype)):
            if(genotype[k]==reg_i):
                neigh = list(graph[units_list[k]]["vizinhos"].keys())
                random.shuffle(neigh)
                
                for n in neigh:
                    if(genotype[graph[n]["seq"]] != reg_i):
                        genotype[k] = genotype[graph[n]["seq"]]
                        break
                
                continue

def repair_steiner(genotype, num_of_regions, min_size, max_size, graph):
    repaired = genotype.copy()
    
    #Step 1: Avoid embedded microregions -> Not implemented
    
    #Step 2: Maintaining the minimum size of a microregion
    i = get_a_region_with_less_units_then_minimum(repaired, min_size)
    while(i is not None):        
        #Merge in an adjacent microregion
        merge_region_in_neighbours(repaired, i, graph)
        
        i = get_a_region_with_less_units_then_minimum(repaired, min_size)
    
    #Step 3: Maintainin the maximum number of microregions
    freq_dict = get_ordered_freq_dict(repaired)
    while len(freq_dict) > num_of_regions:        
        #Get the smallest region
        reg = list(freq_dict.keys())[0]
        #Merge in an adjacent microregion
        merge_region_in_neighbours(repaired, reg, graph)
        
        freq_dict = get_ordered_freq_dict(repaired)
        
    #Step 4: Maintaining the maximum size of a microregion
    freq_dict = get_ordered_freq_dict_desc(repaired)
    reg = list(freq_dict.keys())[0]
    while(freq_dict[reg] > max_size):
        
        #Num units to expell
        units_to_expell = freq_dict[reg] - max_size 
        
        while(units_to_expell > 0):
            
            units_in_region = find_indexes(repaired, reg)
            random.shuffle(units_in_region)
            
            for i in units_in_region:
                if can_remove_unit_from_region(i, reg) and units_to_expell > 0:
                    #Remove unit i to another region
                    chosen_neigh_region = evaluate_region_neighbours(i, reg, repaired, graph, max_size)
                    
                    if chosen_neigh_region is not None:
                        repaired[i] = chosen_neigh_region
                        #One unit minus
                        units_to_expell = units_to_expell - 1 
        
        freq_dict = get_ordered_freq_dict_desc(repaired)
        reg = list(freq_dict.keys())[0]
        
    #Step 5: Transform region indexes
    transform_region_indexes(repaired, num_of_regions)
        
    return repaired