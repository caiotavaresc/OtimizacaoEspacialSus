# -*- coding: utf-8 -*-
"""
This script generates all possible solutions and get the pareto front for the brazilian state of Roraima
"""

from utils import workWithStateData, add_medical_procedures, isFeasible, computeDistanceMatrix, dominance_ranking
from objective_functions import f1_homogeneity_inhabitants, f2_variety_medical_procedures, f3_intra_regional_traveling_distance
import json
"""
State Data
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

NUM_UNITS = len(municipalities)

#Generate all possible solutions
def genRec(ant, pos, solucoes):
    
    for i in range(NUM_OF_REGIONS):
        sol = ant.copy()
        sol.append(i)
        if(pos == NUM_UNITS - 1):
            finSol = {}
            finSol["atribuicao"] = sol
            finSol["contiguo"] = isFeasible(finSol, NUM_OF_REGIONS, municipalities)
            if(finSol["contiguo"] == True):
                finSol['fitness'] = {}
                finSol['fitness']['f1'] = f1_homogeneity_inhabitants(municipalities, sol, NUM_OF_REGIONS, statePop)
                finSol['fitness']['f2'] = f2_variety_medical_procedures(municipalities, sol, NUM_OF_REGIONS)
                finSol['fitness']['f3'] = f3_intra_regional_traveling_distance(municipalities, sol, NUM_OF_REGIONS, distances)
            solucoes.append(finSol)
        else:
            genRec(sol, pos+1, solucoes)
            
        if len(solucoes)%1000 == 0:
            print(len(solucoes), "solutions generated of", NUM_OF_REGIONS ** NUM_UNITS)

solucoes = []
ant = []
genRec(ant, 0, solucoes)

contiguous_solutions = [s for s in solucoes if s['contiguo']]

#Generate the dominance ranking
#List of functions and types
fun_list = {}
fun_list["f1"] = "min"
fun_list["f2"] = "max"
fun_list["f3"] = "min"

cont_sol_dic = {}
for ind in range(len(contiguous_solutions)):
    i = contiguous_solutions[ind]
    str_id = "".join(str(x) for x in i["atribuicao"])
    i["id"] = str_id
    cont_sol_dic[str_id] = i
    
rank = dominance_ranking(contiguous_solutions, fun_list)

#Generate the pareto front
pareto_front = []
for i in range(len(rank)):
    if rank[i] == 1:
        pareto_front.append(contiguous_solutions[i])
        
#Write a json file with the pareto front
with open("pareto_front_rr_steiner.json", "w", encoding="utf-8") as f:
    json.dump(pareto_front, f, indent=4, ensure_ascii=False)