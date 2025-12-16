# module of utility functions not directly related to the planar max cut problem


# test a pair of edges defined by respective start and end points for intersection
# based on segment intersection check, as described 
#    in "Introduction to Algorithms" by Cormen, Leiserson, Rivest & Stein
# adapted to edges by allowing intersections in starting/ending points
def edges_intersect(p1,p2,p3,p4):
	d1 = direction(p3,p4,p1)
	d2 = direction(p3,p4,p2)
	d3 = direction(p1,p2,p3)
	d4 = direction(p1,p2,p4)
	if ((d1 > 0 and d2 < 0) or (d1 < 0 and d2 > 0 )) and ((d3 > 0 and d4 < 0) or (d3 < 0 and d4 > 0 )):
		return True
	elif d1 == 0 and on_segment(p3,p4,p1):
		return True
	elif d2 == 0 and on_segment(p3,p4,p2):
		return True
	elif d3 == 0 and on_segment(p1,p2,p3):
		return True
	elif d4 == 0 and on_segment(p1,p2,p4):
		return True
	return False

def direction(pi,pj,pk):
	xi,yi = pi
	xj,yj = pj
	xk,yk = pk
	return (xk-xi)*(yj-yi) - (xj-xi)*(yk-yi)

def on_segment(pi,pj,pk):
	xi,yi = pi
	xj,yj = pj
	xk,yk = pk
	if min(xi,xj) <= xk <= max(xi, xj) and min(yi, yj) <= yk <= max(yi,yj):
		#allow edges to intersect on their endpoints
		if pi == pk or pj == pk:
			return False
		return True
	return False



# the following union-find implementation was taken from an assignment
# I personally coded and submitted for an Advanced Algorithm homework
#
#
#========== Union-Find =======================================================
# Simple implementation of a union-find datastructure (simple as it will only
# handle the partitioning of integers from 0 to n-1). It will use 
# union-by-size and path compression to improve efficiency as decribed in
# Wikipedia's article.
class SimpleUF:
	def __init__(self, n):
		self._parent = [ None ] * n
		self._size = [ 1 ] * n
	
	def __len__(self):
		return len(self._parent)

	# If the elements are not already in the same part, append the smaller
	# part to the representative of the larger part.
	def union(self, x, y):
		repr_x, repr_y = self.find(x), self.find(y)
		if repr_x == repr_y:
			pass
		elif self._size[repr_x] >= self._size[repr_y]:
			self._parent[repr_y] = repr_x
			self._size[repr_x] += self._size[repr_y]
		else: # <
			self._parent[repr_x] = repr_y
			self._size[repr_y] += self._size[repr_x]

	# Find representative for a certain element. Flatten the path from
	# the element to the representative on each call.
	def find(self, x):
		if self._parent[x] is None:
			return x
		self._parent[x] = self.find(self._parent[x])
		return self._parent[x]
	
	# Check whether 2 elements are in the same part.
	def match(self, x, y):
		return (self.find(x) == self.find(y))
	
	# Give size of part that contains an element.
	def partSize(self, x):
		return self._size[self.find(x)]
#=============================================================================
