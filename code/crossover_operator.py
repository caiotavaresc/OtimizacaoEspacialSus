
"""
Here we reproduce some crossover operators
"""
import random

def crossover_steiner(sol_a, sol_b, num_of_regions, graph, max_num_regions_to_include):
    offspring = []
    b_regions = list(range(num_of_regions))
    offs_gen = sol_a["atribuicao"].copy()
    
    num_regions_to_include = random.randrange(1,max_num_regions_to_include+1)
    
    for reg in range(1, num_regions_to_include+1):
        n_reg = reg * -1
        region_to_include = b_regions[random.randrange(len(b_regions))]
        b_regions.remove(region_to_include)
        
        for i in range(len(sol_b["atribuicao"])):
            if sol_b["atribuicao"][i] == region_to_include:
                offs_gen[i] = n_reg
    
    offspring.append(offs_gen)
    
    return offspring