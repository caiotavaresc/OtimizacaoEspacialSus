# -*- coding: utf-8 -*-
"""
This code is responsible to generate initial solutions for heuristic approaches
"""

import utils
import random

def generate_initial_solution_random(graph, number_of_partitions, min_size, max_size):
    num_units = len(graph)
    solution = [-1] * num_units
    list_units = list(graph.keys())
    non_used_units = list(graph.keys())
    regions = {}
    expansible_regions = []
    
    #seeding phase
    for i in range(number_of_partitions):
        chosen = random.randrange(0, len(non_used_units))
        
        #pegar a unidade, e remover da lista
        unit = non_used_units.pop(chosen)
        
        #atribuir a raiz
        solution[list_units.index(unit)] = i
        
        regions[i] = {}
        regions[i]["units"] = [unit]
        regions[i]["expansible_units"] = [unit]
        expansible_regions.append(i)
        
    #expansion phase
    while len(non_used_units) > 0:
        
        #If there's no regions to expand, return an infeasible solution
        if len(expansible_regions) == 0:
            return [0] * num_units
        
        #select next region to be expanded
        next_region = expansible_regions[random.randrange(len(expansible_regions))]
        #select an unused neighbor
        region_units = regions[next_region]["expansible_units"]
        unit_to_expand = region_units[random.randrange(len(region_units))]
        neighbors = list(graph[unit_to_expand]['vizinhos'].keys())
        
        #choose a randon unused neighbor
        chosen_neighbor = None
        while len(neighbors) > 0:
            next_neighbor = neighbors.pop(random.randrange(len(neighbors)))
            if next_neighbor in non_used_units:
                chosen_neighbor = next_neighbor
                break
        
        if chosen_neighbor is not None:
            index = list_units.index(chosen_neighbor)
            solution[index] = next_region
            regions[next_region]["units"].append(chosen_neighbor)
            regions[next_region]["expansible_units"].append(chosen_neighbor)
            del non_used_units[non_used_units.index(chosen_neighbor)]
            
            #check if region is still expansible
            if len(regions[next_region]["units"]) >= max_size:
                expansible_regions.remove(next_region)
        else:
            regions[next_region]["expansible_units"].remove(unit_to_expand)
            if(len(regions[next_region]["expansible_units"]) == 0):
                expansible_regions.remove(next_region)
    
    return solution