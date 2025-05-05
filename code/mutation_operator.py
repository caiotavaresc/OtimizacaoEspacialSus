
"""
Here we reproduce some mutation operators
"""
import random
from collections import Counter

def mutation_steiner(genotype, graph, max_moving_units):
    mutant = genotype.copy()
    units = list(graph.keys())
    
    largest_region = Counter(mutant).most_common(1)[0][0]
    units_in_reg = [i for i, val in enumerate(mutant) if val == largest_region]
    
    #Select a boundary unit from the largest region to move
    num_units_to_move = random.randint(1, max_moving_units)
    moved_units = 0
    
    while(moved_units < num_units_to_move and len(units_in_reg) > 0):
        unit_to_move = units_in_reg[random.randrange(len(units_in_reg))]
        neighbors = list(graph[units[unit_to_move]]["vizinhos"].keys())
        
        chosen_neighbor = None
        while(chosen_neighbor is None and len(neighbors) > 0):
            candidate = neighbors[random.randrange(len(neighbors))]
            candidate_index = graph[candidate]["seq"]
            if(mutant[candidate_index] != largest_region):
                chosen_neighbor = candidate_index
            else:
                neighbors.remove(candidate)
        
        if(chosen_neighbor is None):
            units_in_reg.remove(unit_to_move)
        else:
            mutant[unit_to_move] = mutant[chosen_neighbor]
            moved_units = moved_units + 1
            units_in_reg.remove(unit_to_move)
    
    return mutant