# -*- coding: utf-8 -*-
"""
Created on Sat May  3 22:56:47 2025

@author: ctcca
"""
import random
from initial_solutions_generator import generate_initial_solution_random
from utils import workWithStateData, add_medical_procedures, isFeasible, computeDistanceMatrix, dominance_ranking, calculate_crowding_distance, crowded_binary_tournament
from crossover_operator import crossover_steiner
from mutation_operator import mutation_steiner
from repair_operator import repair_steiner

"""
GA Parameters
"""
population_size = 100
mating_pool_size = population_size
random_seed_min = 0.01
random_seed_max = 0.30
crossover_prob = 0.9
mutation_prob_min = 0
mutation_prob_max = 0.05
number_of_generations = 500

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
GA Objective Functions
"""

def f1_homogeneity_inhabitants(graph, sol, num_regions, totalPop):
    units_list = list(graph.keys())
    #Get each region population
    pops = {}
    
    for i in range(num_regions):
        pops[i] = 0
    
    for i in range(len(sol)):
        pops[sol[i]] = pops[sol[i]] + graph[units_list[i]]['NUM_HABITANTES']
    
    mean_pop = totalPop/num_regions
    desvio_total = 0
    
    for i in pops:
        desvio = abs(pops[i] - mean_pop)
        desvio_total = desvio_total + desvio
        
    desvio_medio = desvio_total/num_regions
    
    return desvio_medio


def f2_variety_medical_procedures(graph, sol, num_regions):
    units_list = list(graph.keys())
    procedures = {}
    
    for i in range(num_regions):
        procedures[i] = set()
        
    for i in range(len(sol)):
        for j in graph[units_list[i]]["procedimentos"]:
            procedures[sol[i]].add(j)
    
    sum_procs_reg = 0
    
    for i in range(num_regions):
        sum_procs_reg = sum_procs_reg + len(procedures[i])
        
    valor_medio = sum_procs_reg / num_regions
    
    return valor_medio


def f3_intra_regional_traveling_distance(graph, sol, num_regions):
    min_dist_region = {}
    units_list = list(graph.keys())
    
    for i in range(num_regions):
        min_dist_region[i] = float("inf")
    
    for i in range(len(sol)):
        for j in range(len(sol)):
            if sol[i] != sol[j]:
                dist = distances[(units_list[i],units_list[j])]
                if dist < min_dist_region[sol[i]]:
                    min_dist_region[sol[i]] = dist
                if dist < min_dist_region[sol[j]]:
                    min_dist_region[sol[j]] = dist
    
    sum_min_dist = 0
    for k in min_dist_region:
        sum_min_dist = sum_min_dist + min_dist_region[k]
    
    return sum_min_dist

"""
GA Initial Solution Generator
"""
def gen_solution(graph, num_regions, totalPopulation):
    sol = generate_initial_solution_random(graph, num_regions, 2, len(graph)-((num_regions-1)*2))
    
    solucao = {}
    solucao['atribuicao'] = sol
    solucao['fitness'] = {}
    solucao['fitness']['f1'] = f1_homogeneity_inhabitants(graph, sol, num_regions, totalPopulation)
    solucao['fitness']['f2'] = f2_variety_medical_procedures(graph, sol, num_regions)
    solucao['fitness']['f3'] = f3_intra_regional_traveling_distance(graph, sol, num_regions)
    
    return solucao

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
        
#Evolution
gen_atu = 1
#while(gen_atu <= number_of_generations):
    
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
new_item = crossover_steiner(population[mating_pool[0]], population[mating_pool[1]], NUM_OF_REGIONS, municipalities, 1)

print("Pai 1", population[mating_pool[0]]["atribuicao"])
print("Pai 2", population[mating_pool[1]]["atribuicao"])
print("filho", new_item[0])

mutant = mutation_steiner(new_item[0], municipalities, 3)
print("mutante", mutant)

#Update generator indicator
gen_atu = gen_atu + 1