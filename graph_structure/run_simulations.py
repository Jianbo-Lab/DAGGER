from experiment import *
import numpy as np
from plot import *
import argparse 
import os
def main():
	parser = argparse.ArgumentParser() 
	# Number of replications for the same experiment.
	parser.add_argument('--n_replications',type = int, default = 10) 
	# Which experiment to run.
	parser.add_argument('--experiment',type = str, 
		choices = ['shallow_vs_deep','diamond_vs_hourglass','BH_vs_DAGGER',
		'mountain_vs_valley']) 

	args = parser.parse_args()     

	# Mkdir to put results.
	if 'results' not in os.listdir('./'):
		os.mkdir('results')
	for folder in ['figs','shallow_vs_deep','diamond_vs_hourglass','BH_vs_DAGGER',
	'mountain_vs_valley']:
		if folder not in os.listdir('./results'):
			os.mkdir('results/{}'.format(folder)) 

	if args.experiment == 'shallow_vs_deep': 
		output = shallow_vs_deep(var_name = 'pi0', 
			var_range = np.concatenate((np.array([0.05]), np.arange(0.1,0.9,0.05),
				np.arange(0.90,0.98,0.02), np.array([0.97,0.98]))),
			num_nodes = 500,
			kernel_size = 2, 
			assign_order = 'bottom-up', 
			pi0 = None, 
			signal_strength = 1,  
			alpha = 0.2,  
			n_replications = args.n_replications) 
		print('Data have been saved to results/')

		path = "results/shallow_vs_deep/"
		for file in os.listdir(path):
			if file.endswith(".save.npy"):
				plot_shallow_vs_deep(path, file)
		print('Figures have been saved to results/')

	elif args.experiment == 'diamond_vs_hourglass':
		output = diamond_vs_hourglass(var_name = 'pi0', 
			var_range = np.concatenate((np.arange(0.1,0.95,0.05),
				np.arange(0.91,0.98,0.02))),
			num_nodes = 500,
			kernel_size = 2, 
			assign_order = 'bottom-up', 
			pi0 = None, 
			signal_strength = 1,  
			alpha = 0.2,  
			n_replications = args.n_replications) 
		print('Data have been saved to results/'  )
		path = "results/diamond_vs_hourglass/"
		for file in os.listdir(path):
			if file.endswith(".save.npy"):
				plot_diamond_vs_hourglass(path, file)
		print('Figures have been saved to results/')	

	elif args.experiment == 'BH_vs_DAGGER':
		output = BH_vs_DAGGER(var_range = np.arange(0,1,0.1),
			num_nodes = 200,
			kernel_size = 2, 
			assign_order = 'bottom-up', 
			ratio = 5, 
			signal_strength = 1,  
			alpha = 0.2,  
			n_replications = 10)
		print('Data have been saved to results/')
		path = "results/BH_vs_DAGGER/"
		for file in os.listdir(path):
			if file.endswith(".save.npy"): 
				plot_BH_vs_DAGGER(path, file) 
		print('Figures have been saved to results/')

	elif args.experiment == 'mountain_vs_valley':
		output = mountain_vs_valley(var_name = 'pi0', 
			var_range = np.concatenate((np.arange(0.1,0.95,0.05),
				np.arange(0.91,0.98,0.02))),
			num_nodes = 498,
			kernel_size = 2, 
			assign_order = 'bottom-up', 
			pi0 = None, 
			signal_strength = 1,  
			alpha = 0.2,  
			n_replications = args.n_replications) 

		print('Data have been saved to results/' )
		path = "results/mountain_vs_valley/"
		for file in os.listdir(path):
			if file.endswith(".save.npy"):
				plot_mountain_vs_valley(path, file)
		print('Figures have been saved to results/')

					
if __name__ == '__main__':
	main()

