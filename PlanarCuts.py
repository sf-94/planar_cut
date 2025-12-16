# Module that provides a planar graph implementation,
#  including implementations of random and
#  local search min/max cut algorithms,
#  and an exact min/max cut implementation
#  based on the algorithm by Liers and Pardella
#
# all methods returning a cut will return the set of cut edge ids,
# as well as the set of vertex ids for each partition,
# and the total cut weight



#basic imports
from dataclasses import dataclass
from typing import Any
import random
import itertools
from operator import sub, lt, gt
import math

#External libraries for Delaunay and Matching
import numpy as np
from scipy.spatial import Delaunay
import networkx as nx

#External library to draw svg outputs
import drawsvg as draw

#from MyUnionFind import SimpleUF
from Utility import edges_intersect, SimpleUF


# exception raised in case the random graph generation was unsuccessful,
#  mostly due to floating point precision problems
class RandGraphGenError(Exception):
	pass


# Planar Graph - main class to be imported and used outside the module
@dataclass
class PlGraph: 
	vertices: dict
	edges: dict

	def __len__(self):
		return len(self.vertices)
	

	# draw graph into svg file
	# possible edgeAnnotations "fullweight", "weight" (rounded to 3 decimal digits), and "id"
	def draw(self, edgeAnnotation = "weight", filename = "output.svg"):
		maxX = max([v.pos[0] for v in self.vertices.values()])
		maxY = max([v.pos[1] for v in self.vertices.values()])

		d = draw.Drawing(maxX+20, maxY+20, origin = (-10,-10))

		for e in self.edges.values():
			d.append(draw.Line(*self.vertices[e.vertices[0]].pos,*self.vertices[e.vertices[1]].pos,stroke="lightgray", stroke_width=1))
			

			if edgeAnnotation is not None and edgeAnnotation != "":
				annotation = ""
				if edgeAnnotation == "fullweight":
					annotation = str(e.weight)
				elif edgeAnnotation == "weight":
					annotation = str(round(e.weight,3))
				elif edgeAnnotation == "id":
					annotation = str(e.id)
				d.append(draw.Text(annotation,8,*list(map(lambda a, b: (a+b)*0.5,self.vertices[e.vertices[0]].pos,self.vertices[e.vertices[1]].pos)),fill="lightgrey"))

		for v in self.vertices.values():
			d.append(draw.Circle(*v.pos,3,fill="gray"))
			d.append(draw.Text(str(v.id),8,*map(lambda x: x+3,v.pos)))

		d.set_pixel_scale(2)
		d.save_svg(filename)


	# draw cut into svg file, highlight given cut
	# possible edgeAnnotations "fullweight", "weight" (rounded to 3 decimal digits), and "id"
	def draw_cut(self, cut_edges, part_a, part_b,  edgeAnnotation = "weight", filename = "output_cut.svg"):
		maxX = max([v.pos[0] for v in self.vertices.values()])
		maxY = max([v.pos[1] for v in self.vertices.values()])

		d = draw.Drawing(maxX+20, maxY+20, origin = (-10,-10))

		for e in self.edges.values():
			if e.id in cut_edges: #edge in cut
				d.append(draw.Line(*self.vertices[e.vertices[0]].pos,*self.vertices[e.vertices[1]].pos,stroke="lightgray", stroke_dasharray="10,5", stroke_width=1))
			elif set(e.vertices) < part_a: #edge within partition a
				d.append(draw.Line(*self.vertices[e.vertices[0]].pos,*self.vertices[e.vertices[1]].pos,stroke="red", stroke_width=1))
			elif set(e.vertices) < part_b: #edge within partition b
				d.append(draw.Line(*self.vertices[e.vertices[0]].pos,*self.vertices[e.vertices[1]].pos,stroke="blue", stroke_width=1))
			else: #failsafe
				d.append(draw.Line(*self.vertices[e.vertices[0]].pos,*self.vertices[e.vertices[1]].pos,stroke="lightgray", stroke_width=1))
			
			if edgeAnnotation is not None and edgeAnnotation != "":
				annotation = ""
				if edgeAnnotation == "fullweight":
					annotation = str(e.weight)
				elif edgeAnnotation == "weight":
					annotation = str(round(e.weight,3))
				elif edgeAnnotation == "id":
					annotation = str(e.id)
				d.append(draw.Text(annotation,8,*list(map(lambda a, b: (a+b)*0.5,self.vertices[e.vertices[0]].pos,self.vertices[e.vertices[1]].pos)),fill="lightgrey"))
	
		for v in self.vertices.values():
			if v.id in part_a:
				d.append(draw.Circle(*v.pos,3,fill="darkred"))
			elif v.id in part_b:
				d.append(draw.Rectangle(*map(sub,v.pos,[3,3]),6,6,fill="darkblue"))
			else:
				d.append(draw.Circle(*v.pos,3,fill="gray"))
			d.append(draw.Text(str(v.id),8,*map(lambda x: x+3,v.pos)))

		d.set_pixel_scale(2)
		d.save_svg(filename)


	# provide list representation of random graphs - used by PlGraphs random method
	@classmethod
	def _random_list(cls, numVertices=100, maxX=None, maxY=None, pctEdges=1, calcWeight="random", randMin=0, randMax=1, useDelauney=False, seed=None):

		if maxX is None:
			maxX = numVertices * 10
		
		if maxY is None:
			maxY = numVertices * 5
		
		if seed is not None:
			random.seed(seed)

		#generate list of random points, ensure uniqueness by initially inserting into a set
		random_points = set()
		while len(random_points) < numVertices:
			random_points.add((random.randint(0,maxX),random.randint(0,maxY)))
		random_points = list(random_points)

		#generate random edges, either by delauney or greedy triangulation
		if useDelauney:
			triangulation = Delaunay(np.array(random_points))

			#convert list of triangles represented by points into list of edges 
			random_edges = list(set([tuple(sorted(map(int,e2))) for e2 in sum([[e[i:i+2] for i in range(len(e)-1)]+[[e[-1],e[0]]] for e in triangulation.simplices],[])]))
		else:
			#greedy method: gerate list of all possible potential edges
			all_edges = list(itertools.combinations(range(len(random_points)), 2))
			random_edges = []

			#attempt to include random edge one at a time
			cutoff = len(all_edges)
			while cutoff > 0:
				rand_id = random.randrange(0,cutoff)
				rand_edge = all_edges[rand_id]
				cutoff -= 1
				all_edges[rand_id], all_edges[cutoff] = all_edges[cutoff], all_edges[rand_id]

				#reject edges if they would intersect with any already selected edge
				rejected = False
				for edge in random_edges:
					if edges_intersect(random_points[rand_edge[0]], random_points[rand_edge[1]], random_points[edge[0]], random_points[edge[1]]):
						rejected = True
						break

				if not rejected:
					random_edges += [rand_edge]

		if pctEdges < 1:
			#Union-Find datastructure used to keep track of connected components
			components = SimpleUF(len(random_points))
			cutoff = len(random_edges)
			selected = 0
			#ensure connectiveness by selecting |nodes|-1 edges, ignoring edges between nodes that already are in the same component
			while selected < len(random_points)-1:
					
					# this is one of the locations where floating point precision problems in _random_list can be easily detected
					if selected >= cutoff:
						raise RandGraphGenError

					rand_id = random.randrange(selected,cutoff)
					rand_edge = random_edges[rand_id]	

					if not components.match(*rand_edge):
						random_edges[rand_id], random_edges[selected] = random_edges[selected], random_edges[rand_id]
						selected += 1
						components.union(*rand_edge)
					else:
						cutoff -= 1
						random_edges[rand_id], random_edges[cutoff] = random_edges[cutoff], random_edges[rand_id]

			#select x percent of remaining edges (this means 0 pct will result in a tree that contains only the edges of the last step)
			if pctEdges > 0:
				target = int((len(random_edges)-selected)*pctEdges)+selected
				while selected < target:
					rand_id = random.randrange(selected,len(random_edges))
					random_edges[rand_id], random_edges[selected] = random_edges[selected], random_edges[rand_id]
					selected += 1

			random_edge_subset = random_edges[:selected]


			vlist, elist = random_points, random_edge_subset
		else:
			vlist, elist = random_points, random_edges

		#calculate and append edge weights to list
		for i in range(len(elist)):
			if calcWeight == "distance":
				weight = math.dist(vlist[elist[i][0]],vlist[elist[i][1]])
			else: 
				weight = random.uniform(randMin, randMax)

			elist[i] = tuple(list(elist[i])+[weight])

		return vlist, elist



	# create random graph, by generating a list representation of random graph, and let it be parsed by PlGraphs from_list method
	@classmethod
	def random(cls, numVertices=100, maxX=None, maxY=None, pctEdges=1, calcWeight="random", randMin=0, randMax=1, useDelauney=False, seed=None):
		vlist, elist = cls._random_list(numVertices, maxX, maxY, pctEdges, calcWeight, randMin, randMax, useDelauney, seed)
		return cls.from_list(vlist, elist)



	# provide small test graphs with all positive or all negative edge weights
	@classmethod
	def small_test_graph(cls):
		return cls.from_list([(10,10),(90,10),(100,50),(20,60),(40,70)],[(0,1,80), (1,2,42), (2,0,98), (2,4,63),(3,4,22),(4,0,67)])
	@classmethod
	def small_negative_test_graph(cls):
		return cls.from_list([(10,10),(90,10),(100,50),(20,60),(40,70)],[(0,1,-80), (1,2,-42), (2,0,-98), (2,4,-63),(3,4,-22),(4,0,-67)])


	# create planar graph from list representation (list of x, y Position tuples for vertices, list of u,v, weight triples for edges)
	@classmethod
	def from_list(cls, vlist, elist):
		plg = cls({},{})

		plg.source = vlist, elist

		vdict = {i: PlVertex(i, [], vlist[i],graph=plg) for i in range(len(vlist))}
		edict = {i: PlEdge(i, list(elist[i][:2]), elist[i][2],graph=plg) for i in range(len(elist))}

		# add very edge to their respective vertices
		for e in edict.values():
			vdict[e.vertices[0]].edges += [e.id]
			vdict[e.vertices[1]].edges += [e.id]


		# for every vertex, sort edge list counter-clockwise
		for v in vdict.values():
			# reversed sort necessary, cause pos for drawing has positive y pointing down instead of up
			# v.neighbours.sort(reverse=True,key=lambda neighbour: math.atan2(*reversed(list(map(sub,vdict[neighbour].pos,v.pos)))))

			v.edges.sort(reverse=True,key=lambda edge: math.atan2(*reversed(list(
				map(sub,vdict[edict[edge].vertices[0]].pos, v.pos) if vdict[edict[edge].vertices[0]].pos != v.pos else map(sub,vdict[edict[edge].vertices[1]].pos, v.pos)
				))))
			
			# this is one of the locations where floating point precision problems in _random_list can be easily detected
			if len(v.edges) == 0:
				raise RandGraphGenError

			# prime nextEdge cache for each vertex
			_ = v.nextEdge(v.edges[0])

		plg.vertices, plg.edges = vdict, edict

		return plg
	

	# return list representationn of given graph
	def to_list(self):
		return self.source


	# generate dual graph for given primal planar graph
	def get_dual(self):

		face_id_counter = 0

		# append face information to each edge, counting each face traversed
		for e in self.edges.values():
			if not hasattr(e,"_faces") or e._faces[0] is None:
				e.setFace(e.vertices[0],face_id_counter)
				face_id_counter += 1
			if not hasattr(e,"_faces") or e._faces[1] is None:
				e.setFace(e.vertices[1],face_id_counter)
				face_id_counter += 1

		# create empty graph
		dual = DualGraph([],[])

		# create vertex for each edge found
		for i in range(face_id_counter):
			dual.vertices.append(DualVertex(i, [], dual))
		
		# for each edge, add associated edge to dual, using primal face ids as vertex ids in dual
		for e in self.edges.values():
			dual.edges.append(DualEdge(len(dual.edges), e.id, e._faces, e.weight, dual))

			for v_id in e._faces:
				dual.vertices[v_id].edges.append(len(dual.edges)-1)
		
		# clear face information to avoid conflicts in case of later calls of method
		for e in self.edges.values():
			del e._faces

		return dual


	# return random cut by assigning each vertex to a partition by coin flip
	def random_cut(self, seed=None):
		if seed is not None:
			random.seed(seed)
		
		cut_edges = set()
		part_a = set()
		part_b = set()

		for v in self.vertices.values():
			if random.randrange(2) == 0:
				part_a.add(v.id)
			else:
				part_b.add(v.id)
		
		for e in self.edges.values():
			if e.vertices[0] in part_a and e.vertices[1] in part_b or e.vertices[1] in part_a and e.vertices[0] in part_b:
				cut_edges.add(e.id)

		return cut_edges, part_a, part_b, self._total_cut_weight(cut_edges)

	#approximate min or max cut by local search
	def approx_max_cut(self, seed=None):
		return self._approx_cut(gt,seed)

	def approx_min_cut(self, seed=None):
		return self._approx_cut(lt,seed)

	# helper method, not to be called from outside,
	#  depending on wheather lower than or greater than is given as priority function,
	#  min or max cut is calculated
	def _approx_cut(self,priority,seed=None):
		cut_edges, part_a, part_b, cut_weight = self.random_cut(seed=seed)
		
		if seed is not None:
			random.seed(seed)

		queue = set.union(part_a,part_b)
		waiting = set()

		while len(queue) != 0:
			current = random.choice(tuple(queue))

			new_cut_edges = cut_edges.copy()

			for e in self.vertices[current].edges:
				if e in new_cut_edges:
					new_cut_edges.remove(e)
				else:
					new_cut_edges.add(e)

			new_cut_weight = self._total_cut_weight(new_cut_edges)

			if priority(new_cut_weight, cut_weight):
				cut_edges = new_cut_edges
				cut_weight = new_cut_weight

				if current in part_a:
					part_a.remove(current)
					part_b.add(current)
				else:
					part_b.remove(current)
					part_a.add(current)

				queue.update(waiting)
				waiting.clear()

			queue.remove(current)
			waiting.add(current)

		return cut_edges, part_a, part_b, cut_weight


	# helper method, not to be called from outside, calculates total cut weight
	def _total_cut_weight(self, cut_edges):
		cut_value = 0
		for e in cut_edges:
			cut_value += self.edges[e].weight
		
		return cut_value
	

	# calculates total weight over all edges
	def total_edge_weight(self):
		weight = 0
		for e in self.edges.values():
			weight += e.weight
		
		return weight


	# exact max/min cut algorithm, calculated depending on which matching and priority function
	#  is given to the universal _cut helper function
	def max_cut(self):
		return self._cut(lambda nxg : nx.max_weight_matching(nxg, maxcardinality=True), gt)
	
	def min_cut(self):
		return self._cut(lambda nxg : nx.min_weight_matching(nxg), lt)

	# helper method, not to be called from outside the module
	def _cut(self, match, priority):
		#Step 1: calculate dual graph
		graph = self.get_dual()
		#Step 2: split all nodes of degree > 4
		#  transform into chain of nodes with 0 edges inbetween
		#iterate over indices to leave out vertices added in the process
		for i in range(len(graph.vertices)):
			v = graph.vertices[i]
			if len(v.edges) > 4:
				
				#keep first 3 edges on original node, store rest for moving to new nodes
				moving_edges = v.edges[3:]
				v.edges = v.edges[:3]
				#partition edges into sets to assign to new nodes
				edge_sets = [moving_edges[i:i+2] for i in range(0, len(moving_edges),2)]
				#last node in chain may have 2 or three original edges assigned to
				if len(edge_sets[-1]) == 1:
					edge_sets[-2] = edge_sets[-2]+edge_sets[-1]
					edge_sets = edge_sets[:-1]
				
				#create new vertices and update old edges
				chain_ids = [v.id]
				for eset in edge_sets:
					#create new node, recieve id and add to chain
					new_id = graph.addVertex(eset)
					#track associated vertices for debug purpose
					graph.vertices[new_id]._assocVertex = v.id
					chain_ids += [new_id]
					for e in eset:
						graph.edges[e].updateVertex(v.id, new_id)

				#add zero weight edges inbetween chain nodes
				for i in range(len(chain_ids)-1):
					new_id = graph.addZeroEdge(chain_ids[i:i+2])
					graph.vertices[chain_ids[i]].edges.append(new_id)
					graph.vertices[chain_ids[i+1]].edges.append(new_id)

		#Step 3: turn each node into a K_4
		#iterate over indices to leave out vertices added in the process
		for i in range(len(graph.vertices)):
			v = graph.vertices[i]
			kc_ids = [v.id]
			
			#v will be either degree 3 or 4, so loop will run 2 or 3 times
			for e in v.edges[1:]:
				new_id = graph.addVertex([e])
				graph.vertices[new_id]._assocVertex = v.id
				kc_ids += [new_id]
				graph.edges[e].updateVertex(v.id, new_id)
			
			v.edges = v.edges[:1]
			
			#add one more node to complete K_4 if deg(v) = 3
			if len(kc_ids) == 3:
				new_id = graph.addVertex([])
				graph.vertices[new_id]._assocVertex = v.id
				kc_ids += [new_id]
			
			#fully connect K_4 with zero edges
			for i in range(len(kc_ids)-1):
				for j in range(i+1,len(kc_ids)):
					new_id = graph.addZeroEdge([kc_ids[i],kc_ids[j]])
					graph.vertices[kc_ids[i]].edges.append(new_id)
					graph.vertices[kc_ids[j]].edges.append(new_id)

		lst=[]
		for v in graph.vertices:
			if hasattr(v,"_assocVertex"):
				lst += [(v.id, v._assocVertex)]
			else: lst += [v.id]
		
		#Step 4: Compute min-weight perfect matching
		# build networkX graph data structure, store original edge id in networkX edges
		# to extract cut edges after matching
		nx_graph = nx.Graph()
		for e in graph.edges:
			# as networkX does not support multigraph matching, select single edge from multiedges by priority helper function
			# note: multiedges can only appear within K_4, if dual graph contains self loops
			if nx_graph.has_edge(*e.vertices):
				if priority(e.weight, nx_graph.get_edge_data(*e.vertices)["weight"]):
					nx_graph.add_edge(*e.vertices, weight=e.weight, original_id=e.original_id)
			else:
				nx_graph.add_edge(*e.vertices, weight=e.weight, original_id=e.original_id)
		
		# use given (min/max) matching function to calculate matching
		matching = match(nx_graph)

		# get cut edge id set by extracting original edge id from each edge
		cut_edges = set([nx_graph.get_edge_data(*e)["original_id"] for e in matching])
		# discard None, which will be included due to edges not accociated to original graph edges having original id set to None
		cut_edges.discard(None)


		#include single vertex to prime partitions for DFS traversal
		part_a = set([self.vertices[0].id])
		part_b = set()

		#sorting vertices into partitions by DFS traversal
		stack = [self.vertices[0].id]
		while stack != []:
			current = stack.pop()

			for e in self.vertices[current].edges:
				neighbour = self.edges[e].otherVertex(current)
				if neighbour not in part_a and neighbour not in part_b:
					if current in part_a:
						if e in cut_edges:
							part_b.add(neighbour)
						else:
							part_a.add(neighbour)
					else: #current in part_b
						if e in cut_edges:
							part_a.add(neighbour)
						else:
							part_b.add(neighbour)
					stack.append(neighbour)
		
		return cut_edges, part_a, part_b, self._total_cut_weight(cut_edges)



@dataclass
class PlVertex:
	id: int
	edges: list #sorted counter-clockwise
	pos: tuple = None #for the purpose of drawing, x right, y down
	name: str = ""
	graph: PlGraph = None
	
	def nextEdge(self, edge):
		if not hasattr(self, "_nextEdge"):
			self._nextEdge = {self.edges[i-1]: self.edges[i] for i in range(len(self.edges))}
		return self._nextEdge[edge]



@dataclass
class PlEdge:
	id: int
	vertices: list
	weight: Any = 1
	graph: PlGraph = None

	def setFace(self, source, face_id):
		#track adjacent faces for each edge, with
		#[0] ^= face left of edge from v[0] to v[1]
		#[1] ^= face left of edge from v[1] to v[0]
		if not hasattr(self, "_faces"):
			self._faces = [None,None]
		
		if source == self.vertices[0] and self._faces[0] != face_id:
			self._faces[0] = face_id
			next_source = self.vertices[1]
		elif source == self.vertices[1] and self._faces[1] != face_id:
			self._faces[1] = face_id
			next_source = self.vertices[0]
		else:
			return
		
		next_edge = self.graph.vertices[next_source].nextEdge(self.id)
		self.graph.edges[next_edge].setFace(next_source, face_id)

	def otherVertex(self, source):
		if self.vertices[0] == source:
			return self.vertices[1]
		if self.vertices[1] == source:
			return self.vertices[0]


#===================================================================================================


# helper class, creating a dual graph to PlGraph as part of exact max cut algorithms
#  double functions as working graph for exact max cut algorithms
@dataclass
class DualGraph:
	vertices: list
	edges: list

	def addVertex(self, edges):
		new_id = len(self.vertices)
		self.vertices.append(DualVertex(new_id,edges,self))
		return new_id

	# adding Edge with zero weight
	def addZeroEdge(self, vertices):
		new_id = len(self.edges)
		self.edges.append(DualEdge(new_id,None,vertices,0,self))
		return new_id



@dataclass
class DualVertex:
	id: int
	edges: list
	graph: DualGraph = None



@dataclass
class DualEdge:
	id: int
	original_id: int
	vertices: int
	weight: Any
	graph: DualGraph = None

	def updateVertex(self, old, new):
		if self.vertices[0] == old:
			self.vertices[0] = new
		elif self.vertices[1] == old:
			self.vertices[1] = new