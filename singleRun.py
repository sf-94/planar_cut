# This file gives an example for how to run all Cut algorithms provided by PlanarCuts
# and generates Test outputs and statistics to compare the different algorithms for a randomly generated graph

from PlanarCuts import PlGraph, RandGraphGenError
from tabulate import tabulate



#======== Parameters for graph generation
#number of vertices generated
numVertices = 100
#size of area for vertex positions
maxX = 1000
maxY = 500
#Pct of edges kept from triangulation: 0 = Tree graph, 1 = fully triangulated graph
pctEdges = 1
#possible values: "random", "distance"
calcWeight = "random"
#ignored for "distance"
randMin = 0
randMax = 1
#if True, use delauney triangulation, for False, use greedy triangulation
useDelauney = False

#Number of repetitions for averages of random/approximation algorithms
repetitions = 100



try:
	plg = PlGraph.random(numVertices, maxX, maxY, pctEdges, calcWeight, randMin, randMax, useDelauney)
except RandGraphGenError:
	print("ERROR in random graph generation.\nAvoid using large numbers of vertices on small areas, especially on greedy triangulations, as it can cause floating point precision problems.")
	quit()



plg.draw()


max_cut = plg.max_cut()
plg.draw_cut(*max_cut[:3],filename="output_maxCut.svg")

min_cut = plg.min_cut()
plg.draw_cut(*min_cut[:3],filename="output_minCut.svg")

rand_cut = plg.random_cut()
plg.draw_cut(*rand_cut[:3],filename="output_randomCut.svg")

apprx_max_cut = plg.approx_max_cut()
plg.draw_cut(*apprx_max_cut[:3],filename="output_apprxMaxCut.svg")

apprx_min_cut = plg.approx_min_cut()
plg.draw_cut(*apprx_min_cut[:3],filename="output_apprxMinCut.svg")


avg_rand = rand_max = rand_min = rand_cut[3]
avg_apprx_max = apprx_max_max = apprx_max_min = apprx_max_cut[3]
avg_apprx_min = apprx_min_max = apprx_min_min = apprx_min_cut[3]
for i in range(repetitions-1):
	rand_current = plg.random_cut()[3]
	avg_rand += rand_current
	rand_max = max(rand_max, rand_current)
	rand_min = min(rand_min, rand_current)

	apmax_current = plg.approx_max_cut()[3]
	avg_apprx_max += apmax_current
	apprx_max_max = max(apprx_max_max, apmax_current)
	apprx_max_min = min(apprx_max_min, apmax_current)

	apmin_current = plg.approx_min_cut()[3]
	avg_apprx_min += apmin_current
	apprx_min_max = max(apprx_min_max, apmin_current)
	apprx_min_min = min(apprx_min_min, apmin_current)

avg_rand = avg_rand/repetitions
avg_apprx_max = avg_apprx_max/repetitions
avg_apprx_min = avg_apprx_min/repetitions


total_edge_weight = plg.total_edge_weight()

print("Graph:")
print(plg.to_list())

print()

print("Results for Max-Cut:")
print("Cut Edges: ",max_cut[0],"\nPartition A: ",max_cut[1],"\nPartition B: ",max_cut[2])
print()






table = [
	["Maximum-Cut",max_cut[3]],
	["Minimum-Cut",min_cut[3]],
	["Random-Cut",rand_cut[3],avg_rand,rand_min,rand_max],
	["Apprx-Max-Cut",apprx_max_cut[3],avg_apprx_max,apprx_max_min,apprx_max_max],
	["Apprx-Min-Cut",apprx_min_cut[3],avg_apprx_min,apprx_min_min,apprx_min_max]
]




print("Total Edge Weight: "+str(total_edge_weight)+"\n")
print(tabulate(table, headers=["Algorithm","Total\nCut Weight'","Avg. Cut\n Weight*","Min. Cut\nWeight*","Max. Cut\nWeight*"]))
print("('): result for a single run of Random and Approximations, \n        for which graphic output is provided")
print("(*): over "+str(repetitions)+" repetitions\n")



table2 = [
	["Random","Avg",avg_rand,100*avg_rand/total_edge_weight],
	["~","Min",rand_min,100*rand_min/total_edge_weight],
	["~","Max",rand_max,100*rand_max/total_edge_weight],
	[],
	["Apprx-Max-Cut","Avg",avg_apprx_max,100*avg_apprx_max/total_edge_weight,(100*avg_apprx_max/max_cut[3] if max_cut[3] !=0 else None)],
	["~","Min",apprx_max_min,100*apprx_max_min/total_edge_weight,(100*apprx_max_min/max_cut[3] if max_cut[3] !=0 else None)],
	["~","Max",apprx_max_max,100*apprx_max_max/total_edge_weight,(100*apprx_max_max/max_cut[3] if max_cut[3] !=0 else None)]
]

print("Random and Approx. Results compared to correct cut weight and total edge weight:")
print(tabulate(table2,headers=["","\n(*)","Cut  \nWeight","%Total \nE-Weight","%Correct \nCut Weight"]))
print("(*): over "+str(repetitions)+" repetitions\n")