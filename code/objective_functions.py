# -*- coding: utf-8 -*-
"""
This script contains the objective functions
"""



"""
GA Objective Functions (Steiner)
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


def f3_intra_regional_traveling_distance(graph, sol, num_regions, distances):
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