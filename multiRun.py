# This file gives an example for how to run max Cut algorithms provided by PlanarCuts.
# It logs statistics on how max cut approximations perform over large groups of random graphs
# with similar parameters.


from PlanarCuts import PlGraph, RandGraphGenError
from pathlib import Path

#======== Parameters for graph generation
#number of vertices generated
numVertices = 100
#size of area for vertex positions
maxX = 1000
maxY = 1000
#Pct of edges kept from triangulation: 0 = Tree graph, 1 = fully triangulated graph
pctEdges = 1
#possible values: "random", "distance"
calcWeight = "random"
#ignored for "distance"
randMin = 0
randMax = 1
#if True, use delauney triangulation, for False, use greedy triangulation
useDelauney = False


#Number of distinct graphs randomly generated
numGraphs = 100
#Number of repetitions for averages of random/approximation algorithms per graph
repetitions = 100


Path("logs").mkdir(exist_ok=True)
#Filename for detailed log for each graph in this specific multi run (will be overwritten if name exists)
#log_this = "tmp.log"
log_this = "logs/"+"-".join(("run",str(numVertices),str(pctEdges),str(calcWeight),"delauney" if useDelauney else "greedy"))+".log"
#Filename to log all multi runs (result summary will be appended)
log_all = "logs/overall.log"


#number of failed random graph generation attempts before program is stopped
allowedFails = int(numGraphs/10)


file_this = open(log_this,"w")

overall_avg_pct = 0
overall_max_pct = 0
overall_min_pct = 100
overall_avg_num_edges = 0
overall_avg_total_weight = 0

counter = 0
fail_counter = 0
while counter < numGraphs:

	try:
		plg = PlGraph.random(numVertices, maxX, maxY, pctEdges, calcWeight, randMin, randMax, useDelauney)
	except RandGraphGenError:
		if fail_counter >= allowedFails:
			print("Aborted for exceeding allowed number of graph generation failures.\nAvoid using large numbers of vertices on small areas, especially on greedy triangulations, as it can cause floating point precision problems.")
			break
		fail_counter += 1
		continue

	
	max_cut = plg.max_cut()

	total_weight = plg.total_edge_weight()
	cut_weight = max_cut[3]
	num_edges = len(plg.edges)

	overall_avg_total_weight += total_weight
	overall_avg_num_edges += num_edges

	avg_apprx = min_apprx = max_apprx = plg.approx_max_cut()[3]
	for i in range(repetitions-1):
		apmax_current = plg.approx_max_cut()[3]
		avg_apprx += apmax_current
		max_apprx = max(max_apprx, apmax_current)
		min_apprx = min(min_apprx, apmax_current)

	avg_apprx = avg_apprx/repetitions
	
	avg_apprx_pct = 100*avg_apprx/cut_weight
	max_apprx_pct = 100*max_apprx/cut_weight
	min_apprx_pct = 100*min_apprx/cut_weight


	overall_avg_pct += avg_apprx_pct
	overall_max_pct = max(overall_max_pct,max_apprx_pct)
	overall_min_pct = min(overall_min_pct,min_apprx_pct)

	file_this.write(",".join(map(str, (num_edges,total_weight,cut_weight,avg_apprx,avg_apprx_pct,min_apprx,min_apprx_pct,max_apprx,max_apprx_pct)))+"\n")

	counter += 1

	print(counter)

overall_avg_pct = overall_avg_pct/counter
overall_avg_total_weight = overall_avg_total_weight/counter
overall_avg_num_edges = overall_avg_num_edges/counter



file_all = open(log_all,"a")

file_all.write(",".join(map(str, (
	numVertices, maxX, maxY, pctEdges, calcWeight, randMin, randMax, useDelauney,
	counter, overall_avg_num_edges, overall_avg_total_weight, overall_avg_pct, overall_min_pct, overall_max_pct
)))+"\n")


file_this.close()
file_all.close()

print("Ran for "+str(counter)+" graphs with "+str(fail_counter)+" graph generation failures")