from make_data import *
import sys 
sys.path.append('../') 
from core.utils import estimate_fdr_and_power_with_arrays
from core.algorithm import find_children_parents, find_depths,DAGGER_topo 
import statsmodels.sandbox.stats.multicomp 
import numpy as np

def run_single_experiment(nums_by_depth, 
	kernel_size = 3, 
	assign_order = 'bottom-up', 
	pi0 = 0.7, 
	signal_strength = 1, 
	signal_vs_depth = 'reciprocal', 
	alpha = 0.05,   
	n_replications = 10,
	fix_graph = True,
	if_seq = False):
	"""
	This function runs a single experiment for n_replications times. 

	Output:

	output: a [3,4] array, each row records the mean, 
	std of fdrs, the mean, std of powers respectively. 
	Three rows are for our algorithm, BY, BH respectively. 
	"""
	fdrs = np.zeros((n_replications,2))
	powers = np.zeros((n_replications,2))

	
	if fix_graph:
		adj_matrix = generate_graph(nums_by_depth, 
			kernel_size, 
			assign_order = assign_order, 
			seed = np.random.randint(100)) 
		children_list, parents_list = find_children_parents(adj_matrix) 

		depths = find_depths(parents_list) 
		hypos = assign_hypothesis(children_list, pi0, 
				seed = np.random.randint(100))

		pi0_ = 1 - sum(hypos)/float(len(hypos))
		print 'The proportion of nulls is {}'.format(pi0_)
	else:
		avg_prop_H1 = 0.0

	for r in range(n_replications):  
		pvalues = compute_pvalues(depths, hypos, signal_strength, 
			signal_vs_depth, seed = r) 
		rejections, _ = DAGGER_topo(adj_matrix, pvalues, alpha, 
			reshaping_func = 'ID') 


		bh_rejections = statsmodels.sandbox.stats.multicomp.\
		multipletests(pvalues,alpha,'fdr_bh')[0]

		# Estimate fdr and power. 
		fdrs[r,0], powers[r,0] = estimate_fdr_and_power_with_arrays(rejections,
		 hypos)  
		fdrs[r,1], powers[r,1] = estimate_fdr_and_power_with_arrays(bh_rejections,
		 hypos)  

	if not fix_graph:
		print 'The proportion of H1 is {}.'.format(avg_prop_H1 / float(n_replications)) 

	output = np.zeros((2,4))
	output[:,0] = np.mean(fdrs,axis = 0)
	output[:,1] = np.std(fdrs,axis = 0)
	output[:,2] = np.mean(powers,axis = 0)
	output[:,3] = np.std(powers,axis = 0) 
	if not fix_graph:
		return output, 0
	else:
		return output, pi0_ 

def shallow_vs_deep(var_name = 'pi0', 
	var_range = np.arange(0.05,1,0.05),
	num_nodes = 100,
	kernel_size = 3, 
	assign_order = 'bottom-up', 
	pi0 = 0.7, 
	signal_strength = 2,  
	alpha = 0.2,  
	n_replications = 10):
	nums_by_depth_shallow = np.array([num_nodes/2,num_nodes/2])
	nums_by_depth_deep = np.array([num_nodes/4] * 4)
	# np.array([num_nodes/15 * i for i in range(1,6)]) 

	# Shallow, deep, bh respectively. 
	#[l,3,4]
	output = np.zeros((len(var_range), 2, 4))
	#[l,2]
	pi0s = np.zeros((len(var_range),2))

	for i, var in enumerate(var_range):
		if var_name == 'pi0':
			pi0 = var
		elif var_name == 'kernelsize':
			print 'not implemented for true pi0'
			return 
			kernel_size = var 

		shallow_output, pi0s[i, 0] = run_single_experiment(nums_by_depth_shallow, 
			kernel_size = kernel_size, 
			assign_order = assign_order, 
			pi0 = pi0, 
			signal_strength = signal_strength, 
			signal_vs_depth = np.array([2, 2]), 
			alpha = alpha,   
			n_replications = n_replications) 

		output[i, 0, :] = shallow_output[0, :]

		deep_output, pi0s[i, 1] = run_single_experiment(nums_by_depth_deep, 
			kernel_size = kernel_size, 
			assign_order = assign_order, 
			pi0 = pi0, 
			signal_strength = signal_strength, 
			signal_vs_depth = np.array([2, 2, 2, 2]), 
			alpha = alpha,   
			n_replications = n_replications)	

		output[i, 1, :] = deep_output[0, :]  

	output = np.concatenate((output, pi0s.reshape(len(var_range),2,1)), -1) 

	name = 'results/shallow_vs_deep/{}-k-{}-order-{}-pi0-{}-signal-{}-alpha-{}-num-{}-rep-{}'.format(var_name, kernel_size, assign_order, 
		pi0, signal_strength, alpha, num_nodes, n_replications)

	np.save(name + '_output.save',output)
	return output

def diamond_vs_hourglass(var_name = 'pi0', 
	var_range = np.arange(0.1,1,0.1),
	num_nodes = 100,
	kernel_size = 3, 
	assign_order = 'bottom-up', 
	pi0 = 0.7, 
	signal_strength = 2,  
	alpha = 0.2,  
	n_replications = 10):
	nums_by_depth_diamond = np.array([num_nodes/4,num_nodes/4*2,num_nodes/4])
	nums_by_depth_hourglass = np.array([num_nodes/5*2,num_nodes/5,num_nodes/5*2])

	# Shallow, deep, bh respectively. 
	output = np.zeros((len(var_range), 2, 4))

	pi0s = np.zeros((len(var_range),2))

	for i, var in enumerate(var_range):
		if var_name == 'pi0':
			pi0 = var 
		diamond_output, pi0s[i, 0] = run_single_experiment(nums_by_depth_diamond, 
			kernel_size = [1, 2], 
			assign_order = assign_order,#'diamond', 
			pi0 = pi0, 
			signal_strength = signal_strength, 
			signal_vs_depth = np.array([2,2,2]), 
			alpha = alpha,   
			n_replications = n_replications)
		output[i, 0, :] = diamond_output[0, :]

		hourglass_output, pi0s[i, 1] = run_single_experiment(nums_by_depth_hourglass, 
			kernel_size = [2, 1], 
			assign_order = assign_order,#'hourglass', 
			pi0 = pi0, 
			signal_strength = signal_strength, 
			signal_vs_depth = np.array([2,2,2]), 
			alpha = alpha,   
			n_replications = n_replications)	 
		output[i, 1, :] = hourglass_output[0,:]

	output = np.concatenate((output, pi0s.reshape(len(var_range),2,1)),-1)

	name = 'results/diamond_vs_hourglass/{}-k-{}-order-{}-pi0-{}-signal-{}-alpha-{}-num-{}-rep-{}'.format(var_name,
		kernel_size, assign_order, pi0, signal_strength, alpha, num_nodes, n_replications)

	np.save(name + '_output.save',output)

def mountain_vs_valley(var_name = 'pi0', 
	var_range = np.arange(0.1,1,0.1),
	num_nodes = 100,
	kernel_size = 3, 
	assign_order = 'bottom-up', 
	pi0 = 0.7, 
	signal_strength = 2,  
	alpha = 0.2,  
	n_replications = 10):
	nums_by_depth_mountain = np.array([num_nodes/6,num_nodes/6*2,num_nodes/6*3])
	nums_by_depth_valley = np.array([num_nodes/6*3,num_nodes/6*2,num_nodes/6])

	# Shallow, deep, bh respectively. 
	output = np.zeros((len(var_range), 2, 4))

	pi0s = np.zeros((len(var_range),2))

	for i, var in enumerate(var_range):
		if var_name == 'pi0':
			pi0 = var 
		mountain_output, pi0s[i, 0] = run_single_experiment(nums_by_depth_mountain, 
			kernel_size = 1, 
			assign_order = 'bottom-up',#'diamond', 
			pi0 = pi0, 
			signal_strength = signal_strength, 
			signal_vs_depth = np.array([2,2,2]), 
			alpha = alpha,   
			n_replications = n_replications)
		output[i, 0, :] = mountain_output[0, :]

		valley_output, pi0s[i, 1] = run_single_experiment(nums_by_depth_valley, 
			kernel_size = 2, 
			assign_order = 'bottom-up',#'hourglass', 
			pi0 = pi0, 
			signal_strength = signal_strength, 
			signal_vs_depth = np.array([2,2,2]), 
			alpha = alpha,   
			n_replications = n_replications)	 
		output[i, 1, :] = valley_output[0,:]

	output = np.concatenate((output, pi0s.reshape(len(var_range),2,1)),-1)

	name = 'results/mountain_vs_valley/{}-k-{}-order-{}-pi0-{}-signal-{}-alpha-{}-num-{}-rep-{}'.format(var_name,
		kernel_size, assign_order, pi0, signal_strength, alpha, num_nodes, n_replications)

	np.save(name + '_output.save',output)

	return output

def BH_vs_DAGGER(var_range = np.arange(0,1,0.1),
	num_nodes = 100,
	kernel_size = 3, 
	assign_order = 'top-down', 
	ratio = 1, 
	signal_strength = 2,  
	alpha = 0.2,  
	n_replications = 10):
	
	nums_by_depth = np.array([num_nodes/2,num_nodes/2])

	# Shallow, deep, bh respectively. 
	#[l,3,4]
	output = np.zeros((len(var_range), 2, 4))
	#[l,2] 

	for i, var in enumerate(var_range): 
		signal_vs_depth = np.array([signal_strength * ratio, signal_strength])

		shallow_output, _ = run_single_experiment(nums_by_depth,
			kernel_size = kernel_size, 
			assign_order = assign_order, 
			pi0 = var, 
			signal_strength = signal_strength, 
			signal_vs_depth = signal_vs_depth, 
			alpha = alpha,   
			n_replications = n_replications) 

		output[i, :, :] = shallow_output 

	output = np.concatenate((output, np.tile(var_range, (1,2,1)).T),2)
	name = 'results/BH_vs_DAGGER/k-{}-order-{}-signal-{}-alpha-{}-num-{}-rep-{}-ratio-{}'.format(kernel_size, assign_order, 
			signal_strength, alpha, 
			num_nodes, n_replications, ratio)

	np.save(name + '_output.save',output)
	return output

# 	















# 	