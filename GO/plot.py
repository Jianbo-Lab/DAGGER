import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import cPickle as pkl
import argparse
import os
import csv
import pandas as pd
font = {'family' : 'normal',
		'weight' : 'bold',
		'size'   : 30}

matplotlib.rc('font', **font)

def map_to_latex(variable, simplified = False): 
	if simplified:
		if variable == 'pi0':
			return r'$\pi_0$' 
	else: 
		if variable == 'pi0':
			variable_in_word = r'Proportion $\pi_0$ of $H_0$'  

def plot(x,y,title, 
	xlabel, ylabel, ylabels, markers, name, 
	graphtype = 'power', errorbar = None, yerr=None):
	"""
	x: a one-dim array. 
	y: a two-dim array. The second dim is the same as the dim of x.
	yerr: a two-dim array. The second dim is the same as the dim of x.
	ylabels: a one-dim array. The length is the same as the first dim of y.
	""" 
	plt.rc('text', usetex = True)
	plt.rc('font', family = 'serif')   

	fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize = (10,8))  
	for i, y_single in enumerate(y): 
		if 'independent' in title or (ylabels[i] != 'BH.DAG' and ylabels[i] != 'SCR.DAG'):
			if errorbar: 
				ax.errorbar(x, y_single, marker=markers[i], yerr = yerr[i], 
					label = ylabels[i]) 
			else:   
				ax.errorbar(x, y_single, marker=markers[i], markersize=12,
					label = ylabels[i])   

	if graphtype == 'fdr': 
		plt.plot(np.array([x[0],x[-1]]),np.array([0.2,0.2]), 
				linestyle='dashed',
				label = 'Targeted FDR',
				color='black')
		
	ax.legend(loc = 'best',fontsize = 16)
	ax.set_title(title,fontsize = 22)
	ax.set_xlabel(xlabel,fontsize = 28)
	ax.set_ylabel(ylabel,fontsize = 28)
	if graphtype == 'power' or graphtype == 'fdr':
		ax.set_ylim(bottom = 0,top=1)
	elif graphtype == 'time':
		ax.set_yscale('log')
	plt.savefig(name, format='eps', dpi=200) 	 

def extract_float_matrix(filename):
	with open(filename,'r') as f:
		lines = f.read().split('\n')[:-1] 

	mat = [line.split(' ') for line in lines]
	for i in range(len(mat)):
		for j in range(len(mat[0])):
			mat[i][j] = float(mat[i][j])
	mat = np.array(mat) 
	return mat

def extract_output(method, pvalue_type, num, pi0s, set_name):
	result_file = 'output/{}-{}-{}.pkl'.format(pvalue_type,
			 method, set_name)
	if os.path.isfile(result_file):
		with open(result_file,'rb') as f:
			output = pkl.load(f)	
	else:
		output = np.zeros((len(pi0s),num,3))
		for j,pi0 in enumerate(pi0s): 
			if pvalue_type != 'simes' or (method != 'BH.DAG' and method != 'SCR.DAG'): 
				pi0 = '%0.2f'%pi0

				for seed in range(1, num+1):
					filename = '{}_{}_{}_{}_{}.txt'.format(set_name,pi0, pvalue_type, method, seed)
					if filename in os.listdir('output'):  
						output[j, seed-1] = pd.read_table(os.path.join('output',filename), sep='\n', header=None).values.flatten()  
						output[j,seed-1,-1] = max(output[j,seed-1,-1],0.001) 
					else:
						print filename	
		output = np.mean(output, axis=1)
		with open(result_file,'wb') as f:
			pkl.dump(output, f)

	return output	




def main(methods, methods_name, set_name = 'cellprolif', pvalue_type = 'ind'):

	pi0s = np.arange(0.15,1.0,0.05) 
	ys = np.zeros((2,len(methods), len(pi0s), 3))

	# for (pvalue_type in c('independent','simes')){
	
	for s,pvalue_type in enumerate([pvalue_type]):
		for i,method in enumerate(methods):  
			num = 100 if method == 'LORD' or method == 'DAGGER' else 10#or method == 'LORD' else 10
			ys[s,i] = extract_output(method, pvalue_type, num, pi0s, set_name)
			print 'method {} for {} pvalue extracted'.format(method, pvalue_type)
	graph_type = 'subgraph' if set_name == 'cellprolif' else 'full'
	markers = [".","o","v","^","<",">","p","s"]
	for ylabel in ['Power','FDR','Time']:
		index = ['Power','FDR','Time'].index(ylabel)
		for s,pvalue_type in enumerate([pvalue_type]): 
			title = r'{} vs. $\pi_0^L$ on {} GO with {} p-values'.format(ylabel, graph_type, 'Simes' if pvalue_type == 'simes' else 'independent')
			xlabel = r'Proportion of nulls on leaves $\pi_0^L$'

			ylabels = methods_name 
			name = 'results/{}-{}-{}.eps'.format(pvalue_type, graph_type, ylabel) 

			yvalues = ys[s,:,:,index]
			#np.mean(y[pvalue_type][:,:,:,index], 
				#axis = -1)

			plot(pi0s,yvalues,title, xlabel, ylabel, 
				ylabels, markers, name, 
				graphtype = ylabel.lower()) 

# line styles: ['-', '--', '-.', ':']
# markers: ".",".","o","v","^","<",">"
if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument("--set_name", type = str, default = 'cellprolif',
		help="which set of genes are being considered.")
	parser.add_argument("--pvalues_type", type = str, default = 'ind',
		help="which set of genes are being considered.")
	
	a = parser.parse_args()
	methods =['Holm','BH','all-goeman','any-goeman','DAGGER']
	methods_name = ['Structured Holm','BH','MG-b1','MG-b2','DAGGER']
	methods += ['LORD']
	methods_name += ['LORD']
	# methods = []
	# methods_name = []
	if a.set_name == 'cellprolif':
		methods += ['SCR.DAG','BH.DAG']
		methods_name += ['SCR-DAG','BH-DAG']

	main(methods, methods_name, a.set_name, a.pvalues_type) 



