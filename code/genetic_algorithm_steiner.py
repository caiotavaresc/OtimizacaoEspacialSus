# -*- coding: utf-8 -*-
"""
Created on Sat May  3 22:56:47 2025

@author: ctcca
"""
import random
from initial_solutions_generator import generate_initial_solution_random
from utils import workWithStateData, add_medical_procedures, isFeasible, checkFeasibility, computeDistanceMatrix, dominance_ranking, calculate_crowding_distance, crowded_binary_tournament, remove_equal_solutions
from crossover_operator import crossover_steiner
from mutation_operator import mutation_steiner
from repair_operator import repair_steiner
from collections import defaultdict
from objective_functions import f1_homogeneity_inhabitants, f2_variety_medical_procedures, f3_intra_regional_traveling_distance

"""
GA Parameters
"""
population_size = 300
mating_pool_size = population_size*2
random_seed_min = 0.01
random_seed_max = 0.30
crossover_prob = 0.9
mutation_prob_min = 0.0001
mutation_prob_max = 0.05
number_of_generations = 100
elitism = True

"""
GA Data
"""
file_state_data = "../data/State_Roraima/Municipios_RR.xlsx"
nrows_file_state_data = 16
file_state_map_json = "../data/State_Roraima/RR_Municipios_2022/RR_Municipios_2022.json"
object_name = "RR_Municipios_2022"
file_state_map_shapefile = "../data/State_Roraima/RR_Municipios_2022/RR_Municipios_2022.shp"
file_medical_procedures = "../data/State_Roraima/Procedimentos_RR.xlsx"
nrows_file_medical_procedures = 964
NUM_OF_REGIONS = 2

stateMap, municipalities, mun_list, stateArea, statePop, arcs = workWithStateData(file_state_data, nrows_file_state_data, file_state_map_json, object_name, file_state_map_shapefile)
add_medical_procedures(municipalities, file_medical_procedures, nrows_file_medical_procedures)
distances = computeDistanceMatrix(municipalities)

"""
GA Structures
"""
population = []
mating_pool = []

"""
GA Initial Solution Generator
"""
def gen_solution(graph, num_regions, totalPopulation):
    sol = generate_initial_solution_random(graph, num_regions, 1, len(graph)-(num_regions-1))
    
    solucao = {}
    solucao['atribuicao'] = sol
    solucao['fitness'] = {}
    solucao['fitness']['f1'] = f1_homogeneity_inhabitants(graph, sol, num_regions, totalPopulation)
    solucao['fitness']['f2'] = f2_variety_medical_procedures(graph, sol, num_regions)
    solucao['fitness']['f3'] = f3_intra_regional_traveling_distance(graph, sol, num_regions, distances)
    solucao["id"] = "".join(str(x) for x in sol)
    return solucao

def mount_solution_from_genotype(genotype, graph, num_regions, totalPopulation):
    solucao = {}
    solucao['atribuicao'] = genotype
    solucao['fitness'] = {}
    solucao['fitness']['f1'] = f1_homogeneity_inhabitants(graph, genotype, num_regions, totalPopulation)
    solucao['fitness']['f2'] = f2_variety_medical_procedures(graph, genotype, num_regions)
    solucao['fitness']['f3'] = f3_intra_regional_traveling_distance(graph, genotype, num_regions, distances)
    solucao["id"] = "".join(str(x) for x in genotype)
    
    return solucao

def choose_solutions_from_dominance_and_diversity_rank(new_solutions, new_pop_size):
    new_pop_ranking = dominance_ranking(new_solutions, fun_list)
    
    positions = defaultdict(list)
    
    for idx, val in enumerate(new_pop_ranking):
        positions[val].append(idx)

    # Ordena o dicionário pelas chaves
    ord_positions = dict(sorted(positions.items()))
    
    ranks_to_copy = []
    ranks_to_order_diversity = []
    
    ranks_keys = list(ord_positions.keys())
    
    control_size = 0
    num_to_add = 0
    while True:
        next_one = ranks_keys.pop(0)
        
        if(control_size + len(ord_positions[next_one]) <= new_pop_size):
            #Add the solutions
            ranks_to_copy.append(next_one)
            control_size = control_size + len(ord_positions[next_one])
        else:
            #order solutions, add by diversity ranking
            ranks_to_order_diversity.append(next_one)
            num_to_add = new_pop_size - control_size
            break
    
    finalPop = []
    
    #Copy first
    for ind in ranks_to_copy:
        for pos in ord_positions[ind]:
            finalPop.append(new_solutions[pos])
            
    #Mount the diversity rank
    if len(ranks_to_order_diversity) > 0:
        fitness_keys = list(fun_list.keys())
        new_pop_diversity = calculate_crowding_distance(new_solutions, new_pop_ranking, fitness_keys)
        residual = []
        
        for pos in ord_positions[ranks_to_order_diversity[0]]:
            new_solutions[pos]["diversity"] = new_pop_diversity[pos]
            residual.append(new_solutions[pos])
        
        #Order residual from diversity
        residual = sorted(residual, key=lambda x: x['diversity'], reverse=True)

        for i in range(num_to_add):
            sol = residual.pop(0)
            del sol["diversity"]
            finalPop.append(sol)
    
    return finalPop

"""
GA EXEC
"""

#Generate initial solutions
while(len(population) < population_size):
    solTemp = gen_solution(municipalities, NUM_OF_REGIONS, statePop)
    if(isFeasible(solTemp, NUM_OF_REGIONS, municipalities)):
        population.append(solTemp)
        
#List of functions and types
fun_list = {}
fun_list["f1"] = "min"
fun_list["f2"] = "max"
fun_list["f3"] = "min"
fitness_keys = list(fun_list.keys())

worst = {}
worst["f1"] = 0
worst["f2"] = float("inf")
worst["f3"] = 0
        
#Evolution
gen_atu = 1
while(gen_atu <= number_of_generations):
    print("gen", gen_atu)
    
    #Selection operator - Mount the mating pool
    #Get two random solutions and copy the best one, 
    #based on the convergence and diversity of the solutions
    
    #Convergence ranking
    sol_ranking = dominance_ranking(population, fun_list)
    
    #Diversity
    sol_diversity = calculate_crowding_distance(population, sol_ranking, fitness_keys)
    
    #Mount the mating pool
    mating_pool = []
    while(len(mating_pool) < mating_pool_size):
        sol_a = random.randrange(len(population))
        sol_b = random.randrange(len(population))
        winner = crowded_binary_tournament(sol_a, sol_b, sol_ranking, sol_diversity)
        mating_pool.append(winner)
    
    new_population = []
    feasible_ones = []
    infeasible_ones = []
    
    while(len(mating_pool) > 0):
        pai1 = population[mating_pool.pop()]
        pai2 = population[mating_pool.pop()]
        
        if random.random() < crossover_prob:
            new_items = crossover_steiner(pai1, pai2, NUM_OF_REGIONS, municipalities, 1)
        else:
            new_items = [pai1["atribuicao"], pai2["atribuicao"]]
            
        for i in range(len(new_items)):
            new_item = new_items[i]
            
            #Mutation probability
            mutation_prob = random.uniform(mutation_prob_min, mutation_prob_max)
            
            if random.random() < mutation_prob:
                new_item = mutation_steiner(new_item, municipalities, 3)
            
            new_item = repair_steiner(new_item, NUM_OF_REGIONS, 1, len(municipalities)-(NUM_OF_REGIONS-1), municipalities)
            
            new_solution = mount_solution_from_genotype(new_item, municipalities, NUM_OF_REGIONS, statePop)
            
            solutionIsFeasible, distance = checkFeasibility(new_solution, NUM_OF_REGIONS, municipalities)
            if solutionIsFeasible:
                feasible_ones.append(new_solution)
            else:
                new_solution["distance"] = distance
                infeasible_ones.append(new_solution)
                
    #Elite preservation
    if elitism:
        new_population.extend(population)
        new_population.extend(feasible_ones)
    else:
        new_population.extend(feasible_ones)
    
    #Penalty-parameter-less-constraint handling approach
    if(len(infeasible_ones) > 0):
        for i in new_population:
            for fun in fun_list:
                if(fun_list[fun]=="min"):
                    if(i["fitness"][fun] > worst[fun]):
                        worst[fun] = i["fitness"][fun]
                else:
                    if(i["fitness"][fun] < worst[fun]):
                        worst[fun] = i["fitness"][fun]
    
        for i in infeasible_ones:
            for fun in fun_list:
                if(fun_list[fun]=="min"):
                    i["fitness"][fun] = worst[fun] + i["distance"]
                else:
                    i["fitness"][fun] = worst[fun] - i["distance"]
                    
            del i["distance"]
            
        #Add the infeasible solutions in new population
        new_population.extend(infeasible_ones)
    
    #Make the new population
    fp = choose_solutions_from_dominance_and_diversity_rank(new_population, population_size)
    random.shuffle(fp)
    population = fp
    
    #Update generator indicator
    gen_atu = gen_atu + 1
    
#remove equal solutions
final_pop = remove_equal_solutions(population)
sol_ranking = dominance_ranking(final_pop, fun_list)