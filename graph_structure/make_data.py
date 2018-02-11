import numpy as np 
from scipy.stats import norm 
import sys 
sys.path.append('../') 
from core.algorithm import find_leaves 

# Create Graph.
def generate_graph(nums_by_depth, kernel_size, assign_order = 'bottom-up', seed = 0):
	"""
	This function generates an adjacency matrix of a DAG. 

	Inputs: 
	nums_by_depth: an array whose ith element indicating 
	the number of nodes at the (i+1)th depth.
	
	kernel_size: Indicating how many parents a child is randomly connected to 
	if assign_order is bottom-up; how many children a parent is randomly
	connected to if assign_order is top-down. 

	assign_order: Whether one assigns edges top-down or bottom-up.

	seed: random seed.

	Output:
	adjacency matrix of the graph. 
	"""
	np.random.seed(0)
	num_nodes = sum(nums_by_depth)  

	D = len(nums_by_depth)

	# The ith entry is the last node's index at depth i. 
	ind_last_nodes = np.concatenate(([0],np.cumsum(nums_by_depth)))

	adj_matrix = np.tile(0,(num_nodes,num_nodes))

	if assign_order == 'top-down':
		# Top down seems not a good approach to specify 
		# the number of nodes at each depth.
		for d_ in range(D-1):
			# At depth d_+1. 
			for node in xrange(ind_last_nodes[d_], ind_last_nodes[d_+1]):
				connected_children = np.random.choice(np.arange(ind_last_nodes[d_+1],
					ind_last_nodes[d_+2]), kernel_size, replace = False)
				adj_matrix[node, connected_children] = 1

			# Avoid a node at depth d from having no parents at depth d-1. 
			for node in xrange(ind_last_nodes[d_+1], ind_last_nodes[d_+2]): 
				if sum(adj_matrix[ind_last_nodes[d_]:ind_last_nodes[d_+1],node]) == 0:
					connected_parents = np.random.choice(np.arange(ind_last_nodes[d_],
						ind_last_nodes[d_+1]),
					max(1, int(kernel_size * nums_by_depth[d_]/nums_by_depth[d_+1])),
					replace = False)
					adj_matrix[connected_parents, node] = 1

	elif assign_order == 'bottom-up':
		if type(kernel_size) != list:
			for d_ in xrange(D-1,0,-1):
				# At depth d_+1. 
				for node in xrange(ind_last_nodes[d_], ind_last_nodes[d_+1]):
					connected_parents = np.random.choice(np.arange(ind_last_nodes[d_-1],
						ind_last_nodes[d_]), kernel_size, replace = False)
					adj_matrix[connected_parents,node] = 1
		else:
			for d_ in xrange(D-1,0,-1):
				# At depth d_+1. 
				for node in xrange(ind_last_nodes[d_], ind_last_nodes[d_+1]):
					connected_parents = np.random.choice(np.arange(ind_last_nodes[d_-1],
						ind_last_nodes[d_]), kernel_size[d_-1], replace = False)
					adj_matrix[connected_parents,node] = 1

				# The below is not necessary.
				# Avoid a node at depth d from having no children at depth d+1. 
				# for node in xrange(ind_last_nodes[d_-1], ind_last_nodes[d_]): 
				# 	if sum(adj_matrix[node, ind_last_nodes[d_]:ind_last_nodes[d_+1]]) == 0:
				# 		connected_children = np.random.choice(np.arange(ind_last_nodes[d_],
				# 			ind_last_nodes[d_+1]),
				# 		max(1, int(kernel_size[d_-1] * nums_by_depth[d_]/nums_by_depth[d_-1])),
				# 		replace = False)
				# 		adj_matrix[node, connected_children] = 1
	elif assign_order == 'hourglass': 
		k = num_nodes / 5
		for node in xrange(ind_last_nodes[0], ind_last_nodes[1]):
			adj_matrix[node, 2*k+node/2] = 1
		for node in xrange(ind_last_nodes[2],ind_last_nodes[3]):
			adj_matrix[(node-3*k) / 2 + 2*k, node] = 1
	elif assign_order == 'diamond':
		k = num_nodes / 4
		for node in xrange(ind_last_nodes[1], ind_last_nodes[2]):
			adj_matrix[(node - k)/2, node] = 1
			adj_matrix[node, (node-k)/2 + 3*k] = 1


	return adj_matrix

# Proportion of nulls on leaves.
def assign_hypothesis(children_list, pi0, seed):
	"""
	This function assigns each node a null or alternative hypothesis.

	Inputs:
	children_list: A list of arrays. The ith item of the list is 
	an array recording the ith node's children.

	pi0: proportion of nulls on the leaves. 

	Output:
	hypos: a 0-1 array indicating whether the hypothesis 
	at node i is null (0) or alternative (1)
	""" 

	num_nodes = len(children_list)

	hypos = np.tile(False,num_nodes) 

	leaves = np.where(find_leaves(children_list)==True)[0] 
	np.random.seed(seed)
	alter_leaves = np.random.choice(leaves, 
		size = int(round((1-pi0) * len(leaves),0)), 
		replace = False)  
	
	hypos[alter_leaves] = True

	for node in range(num_nodes-1,-1,-1):
		if len(children_list[node]) != 0:
			hypos[node] = any(hypos[children_list[node]])

	return hypos


# Compute p-values on all nodes.
def compute_pvalues(depths, hypos, signal_strength, signal_vs_depth, seed = 0):
	"""
	This function computes the p-value at each node. 

	Inputs:
	depths: A list recording the depth of each node.  

	hypos: a 0-1 array indicating whether the hypothesis 
	at node i is null (0) or alternative (1).

	signal_strength: Base signal strength mu0 of alternatives.

	signal_vs_depth: How signal strength varies with depth. 
	choices = ['reciprocal', 'arithmetic', an array that assigns value to each depth]

	seed: random seed.

	Output:
	pvalues
	"""

	np.random.seed(seed)

	num_nodes = len(hypos)

	D = max(depths) 
	if type(signal_vs_depth) == str:
		if signal_vs_depth == 'reciprocal':
			mus = float(signal_strength) / depths 
		elif signal_vs_depth == 'arithmetic':
			mus = float(signal_strength) * (D+1-depths) 
	else:
		signals_by_depths = signal_vs_depth * signal_strength  
		mus = signals_by_depths[depths-1]

	zs = np.random.randn(num_nodes) 

	# random variables at each node.
	rvs = np.where(hypos, zs + mus, zs) 

	pvalues = 1 - norm.cdf(rvs) 

	return pvalues

