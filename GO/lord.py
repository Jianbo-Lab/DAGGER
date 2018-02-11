import numpy as np  
import sys 
sys.path.append('../')
from core.algorithm import find_children_parents
import cPickle as pkl   
import argparse
import time
import pandas as pd
import time
def find_descendants_each_node(children_list):
	num_nodes = len(children_list)
	descendants_list = [[] for i in range(num_nodes)]
	for i in range(num_nodes-1,-1,-1):
		if len(children_list[i]) == 0:
			descendants_list[i] = [i]
		else:
			descendants_list[i] = [i]
			for child in children_list[i]:
				descendants_list[i] += descendants_list[child]

			descendants_list[i] = list(set(descendants_list[i]))

	return descendants_list

def lord(pvalues, adj_matrix, alpha, wealth=None, gammas=None, pvalues_type = 'ind'):
	if wealth == None:
		wealth = 0.5 * alpha
		if pvalues_type == 'simes': 
			wealth_at_previous_rejection = 0.5 * alpha

	if gammas == None:
		gammas = [0.0722*np.log(max(i+1,2))/((i+1)*np.e**(np.sqrt(np.log(i+1)))) for i in range(len(pvalues))]
		if pvalues_type == 'simes':
			gammas = [gamma / (1+np.log(i+1)) for i,gamma in enumerate(gammas)]

	children_list, parents_list = find_children_parents(adj_matrix) 
	descendants_list = find_descendants_each_node(children_list)

	rejections = [False for i in range(len(pvalues))]
	counter = 0
	deleted_nodes = []
	rejected_nodes_indices = []

	for i,pvalue in enumerate(pvalues):
		if i not in deleted_nodes:
			# B = alpha - wealth if len(rejected_nodes_indices) == 0 else alpha
			if pvalues_type == 'ind':
				alphat = gammas[counter] * wealth + sum([(alpha - wealth) * gammas[counter-s] if j == 0 else alpha * gammas[counter-s] for j,s in enumerate(rejected_nodes_indices)]) 
			elif pvalues_type == 'simes':
				alphat = gammas[counter] * wealth_at_previous_rejection
				wealth = wealth - alphat

			if pvalues[i] < alphat:
				rejections[i] = True
				rejected_nodes_indices.append(counter)
				if pvalues_type == 'simes':
					wealth += alpha
					wealth_at_previous_rejection = wealth
			else:
				deleted_nodes += descendants_list[i]
				deleted_nodes = list(set(deleted_nodes))

			counter += 1

	return np.array(rejections)

def extract(seed, pi0, set_name = 'cellprolif',pvalues_type='ind'):   
	pvalues = pd.read_table('data/pvalues_{}_{}_{}_{}.txt'.format(set_name, pvalues_type, seed, pi0), sep=' ',header = None).values  


	hypo_list = pd.read_table('data/hypo_{}_{}_{}_{}.txt'.format(set_name, pvalues_type, seed, pi0), sep=' ',header=None).values

	return pvalues, hypo_list 

def compute_fdp_power(rejections, hypo_list):
	fdp = sum(np.logical_and(hypo_list,rejections)) / float(sum(rejections)) if float(sum(rejections)) != 0.0 else 0.0

	power = sum(np.logical_and(rejections, np.logical_not(hypo_list))) / float(sum(np.logical_not(hypo_list))) if float(sum(np.logical_not(hypo_list))) != 0.0 else 1.0

	return fdp, power
def run_single_experiment(adj_matrix, alpha, pi0, seed, set_name, pvalues_type):
	print "Running algorithm for the {}th replication with alpha = {} and pi0 = {} and set_name = {} and pvalues_type= {}".format(seed, alpha, pi0, set_name, pvalues_type)
	pvalues,hypo = \
			extract(seed,pi0, set_name, pvalues_type)
	pvalues, hypo = pvalues.flatten(), hypo.flatten()
	st = time.time()
	rejections = lord(pvalues,adj_matrix, alpha, 
		pvalues_type=pvalues_type) 
	duration = time.time() - st
	fdp, power = compute_fdp_power(rejections, hypo) 
	return fdp, power, duration

def main(set_name, pvalues_type, seed, pi0):
	num = 100
	alpha = 0.2 
	# alphas = np.arange(0.01,0.33,0.02) 

	adj_matrix = pd.read_table('data/adjmatrix_{}.txt'.format(set_name),
		sep=' ',header=None).values

	pi0 = '{:0.2f}'.format(pi0)

	fdp, power, duration = run_single_experiment(adj_matrix, alpha, pi0, seed, set_name, pvalues_type)
				
	with open('output/{}_{}_{}_{}_{}.txt'.format(set_name, pi0, pvalues_type, 'LORD', seed), 'wb') as f:
		f.write('{}\n{}\n{}'.format(power,fdp,duration)) 

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument("--set_name", type = str, default = 'cellprolif',help="which set of genes are being considered.")
	parser.add_argument("--pvalues_type", type = str, default = 'ind')

	parser.add_argument('--seed', type=int, 
	  default=1) 

	parser.add_argument('--pi0', type=float, 
	  default=0.95)

	a = parser.parse_args()
	main(a.set_name, a.pvalues_type,a.seed, a.pi0) 
