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


def f3_intra_regional_traveling_distance_old(graph, sol, num_regions, distances):
    min_dist_region = {}
    units_list = list(graph.keys())
    
    for i in range(num_regions):
        min_dist_region[i] = float("inf")
    
    for i in range(len(sol)):
        for j in range(len(sol)):
            if sol[i] != sol[j]:
                dist = distances[(units_list[i],units_list[j])]
                comparing_i = (dist * graph[units_list[j]]["NUM_HABITANTES"])
                comparing_j = (dist * graph[units_list[i]]["NUM_HABITANTES"])
                if  comparing_i < min_dist_region[sol[i]]:
                    min_dist_region[sol[i]] = comparing_i
                if comparing_j < min_dist_region[sol[j]]:
                    min_dist_region[sol[j]] = comparing_j
    
    sum_min_dist = 0
    for k in min_dist_region:
        sum_min_dist = sum_min_dist + min_dist_region[k]
    
    return sum_min_dist

           
def f3_intra_regional_traveling_distance(graph, sol, num_regions, distances):
    """
    Calcula a função objetivo f3 para minimizar a distância de viagem inter-regional.

    Parâmetros:
    - graph: dicionário onde cada chave é o código da cidade e contém "NUM_HABITANTES"
    - sol: lista onde sol[i] é a região atribuída à cidade de índice i em list(graph.keys())
    - num_regions: número total de regiões
    - distances: dicionário com tuplas (cidade_i, cidade_j) como chave e distância como valor

    Retorna:
    - valor da função f3
    """
    city_codes = list(graph.keys())  # lista ordenada das cidades
    f3 = 0.0

    for region in range(num_regions):
        for i, city_i in enumerate(city_codes):
            if sol[i] != region:  # cidade i não pertence à região atual
                hi = graph[city_i]["NUM_HABITANTES"]
                min_distance = float('inf')

                for j, city_j in enumerate(city_codes):
                    if sol[j] == region:  # cidade j pertence à região atual
                        # Tenta encontrar a distância (em ambas as ordens)
                        dij = distances.get((city_i, city_j), distances.get((city_j, city_i), None))
                        if dij is not None:
                            weighted_distance = hi * dij
                            if weighted_distance < min_distance:
                                min_distance = weighted_distance

                if min_distance != float('inf'):
                    f3 += min_distance

    return f3
