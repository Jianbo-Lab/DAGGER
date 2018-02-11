import numpy as np 
import sys 
sys.path.append('../')
from core.algorithm import DAGGER_topo, find_leaves_adj
import cPickle as pkl 
from statsmodels.sandbox.stats.multicomp import multipletests 
import argparse
import time
import pandas as pd
import time

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
	rejections_FDR, rejections_on_leaves_FDR = DAGGER_topo(adj_matrix, pvalues, alpha) 
	duration = time.time() - st
	fdp, power = compute_fdp_power(rejections_FDR, hypo) 
	return fdp, power, duration

def main(set_name, pvalues_type, seed, pi0):
	num = 100
	alpha = 0.2 
	# alphas = np.arange(0.01,0.33,0.02)
	pi0s = np.arange(0.15,1.0,0.05) 

	adj_matrix = pd.read_table(set_name+'/adjmatrix.txt',
		sep=' ',header=None).values

	pi0 = '{:0.2f}'.format(pi0)

	fdp, power, duration = run_single_experiment(adj_matrix, alpha, pi0, seed, set_name, pvalues_type)
				
	with open('output/{}_{}_{}_{}_{}.txt'.format(set_name, pi0, pvalues_type, 'DAGGER', seed), 'wb') as f:
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

# if __name__ == '__main__':
# 	set_name = 'full'
# 	main(set_name)