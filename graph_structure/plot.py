import matplotlib.pyplot as plt
import matplotlib
import numpy as np
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
		elif variable == 'kernelsize':
			variable_in_word = r'Kernel size $k$'
		elif variable == 'mu':
			variable_in_word = r'$\mu$'

		return variable_in_word 

def plot_shallow_vs_deep_single(x,y,yerr, title, 
	xlabel, ylabel, name, graphtype='power'):
	"""
	This function plots the figure with error bar 

	and save it to local.
	"""
	plt.rc('text', usetex = True)
	plt.rc('font', family = 'serif') 

	fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize = (10,8)) 
	ax.errorbar(x[:-1,0],y[:-1,0],marker='o',#yerr = yerr[:-1,0], 
		label = 'DAGGER on shallow DAGs',
		color='red')
	ax.errorbar(x[1:,1],y[1:,1],marker='o',#yerr = yerr[1:,1], 
		label = 'DAGGER on deep DAGs',color='blue') 

	if graphtype == 'FDR': 
		if len(x.shape) == 2: 
			plt.plot(np.array([0,1]),np.array([0.2,0.2]), 
				linestyle='dashed',
				label = 'Targeted FDR',
				color='black')
		else:
			plt.plot(np.array([x[0],x[-1]]),np.array([0.2,0.2]), 
				linestyle='dashed',
				label = 'Targeted FDR',
				color='black')
		
	ax.legend(loc = 'upper right',fontsize = 22)
	ax.set_title(title)
	ax.set_xlabel(xlabel,fontsize = 28)
	ax.set_ylabel(ylabel,fontsize = 28)
	ax.set_ylim(bottom = 0,top=1)
	plt.savefig(name, format='eps', dpi=1000) 	


def plot_shallow_vs_deep(path, saved_file):
	output = np.load(path + saved_file)
	subtitle = saved_file.split('_')[0] 

	variable = subtitle.split('-')[0]
	rep = subtitle.split('-')[-1]

	y_fdr, yerr_fdr,y_pow,yerr_pow, var_range = \
	output[:,:,0],output[:,:,1], output[:,:,2], output[:,:,3], output[:,:,4]  

	variable_in_word = map_to_latex(variable) 

	plot_shallow_vs_deep_single(x=var_range,y=y_fdr,yerr=yerr_fdr,
		title = r'FDR vs. {}'.format(variable_in_word), 
		xlabel = variable_in_word,
		ylabel = r'FDR',#r'Achieved FDR',
		name = 'results/figs/shallow-deep-fdr.eps',  
		graphtype='FDR') 

	plot_shallow_vs_deep_single(x=var_range,y=y_pow,yerr=yerr_pow,
		title = r'Power vs. {}'.format(variable_in_word),
		xlabel = variable_in_word,
		ylabel = r'Power',
		name = 'results/figs/shallow-deep-power.eps',  
		graphtype='power') 


def plot_diamond_vs_hourglass_single(x,y,yerr, title, xlabel, ylabel, name,
		  graphtype='power'):
	"""
	This function plots the figure with error bar 

	and save it to local.
	"""
	plt.rc('text', usetex = True)
	plt.rc('font', family = 'serif') 

	fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize = (10,8)) 
	ax.errorbar(x[:,0],y[:,0],marker='o',#yerr = yerr[:,0], 
		label = 'DAGGER on diamond DAGs',
		color='red')
	ax.errorbar(x[:,1],y[:,1],marker='o',#yerr = yerr[:,1], 
		label = 'DAGGER on hourglass DAGs',color='blue') 

	if graphtype == 'FDR': 
		plt.plot(np.array([0, 1]),np.array([0.2,0.2]), 
			 linestyle='dashed',
			 label = 'Targeted FDR',
			color='black')
		
	ax.legend(loc = 'upper right',fontsize = 22)
	ax.set_title(title)
	ax.set_xlabel(xlabel,fontsize = 28)
	ax.set_ylabel(ylabel,fontsize = 28)
	ax.set_ylim(bottom = 0,top=1)
	plt.savefig(name, format='eps', dpi=1000) 	


def plot_diamond_vs_hourglass(path, saved_file):
	output = np.load(path + saved_file)
	subtitle = saved_file.split('_')[0] 

	variable = subtitle.split('-')[0]
	rep = subtitle.split('-')[-1]

	y_fdr, yerr_fdr,y_pow,yerr_pow, var_range = \
	output[:,:,0],output[:,:,1], output[:,:,2], output[:,:,3], output[:,:,4]  

	variable_in_word = map_to_latex(variable) 

	plot_diamond_vs_hourglass_single(x=var_range,y=y_fdr,yerr=yerr_fdr,
		title = r'FDR vs. {}'.format(variable_in_word), 
		xlabel = variable_in_word,
		ylabel = r'FDR',#r'Achieved FDR',
		name = 'results/figs/diamond-hourglass-fdr.eps',  
		 graphtype='FDR') 

	plot_diamond_vs_hourglass_single(x=var_range,y=y_pow,yerr=yerr_pow,
		title = r'Power vs. {}'.format(variable_in_word),
		xlabel = variable_in_word,
		ylabel = r'Power',
		name = 'results/figs/diamond-hourglass-power.eps',  
		 graphtype='power') 


def plot_mountain_vs_valley_single(x,y,yerr, title, xlabel, ylabel, name,
		  graphtype='power'):
	"""
	This function plots the figure with error bar 

	and save it to local.
	"""
	plt.rc('text', usetex = True)
	plt.rc('font', family = 'serif') 

	fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize = (10,8)) 
	ax.errorbar(x[:,0],y[:,0],marker='o',#yerr = yerr[:,0], 
		label = 'DAGGER on mountain DAGs',
		color='red')
	ax.errorbar(x[:-1,1],y[:-1,1],marker='o',#yerr = yerr[:-1,1], 
		label = 'DAGGER on valley DAGs',color='blue') 

	if graphtype == 'FDR': 
		plt.plot(np.array([0, 1]),np.array([0.2,0.2]), 
			 linestyle='dashed',
			 label = 'Targeted FDR',
			color='black')
		
	ax.legend(loc = 'upper right',fontsize = 22)
	ax.set_title(title)
	ax.set_xlabel(xlabel,fontsize = 28)
	ax.set_ylabel(ylabel,fontsize = 28)
	ax.set_ylim(bottom = 0,top=1)
	plt.savefig(name, format='eps', dpi=1000) 	


def plot_mountain_vs_valley(path, saved_file):
	output = np.load(path + saved_file)
	subtitle = saved_file.split('_')[0] 

	variable = subtitle.split('-')[0]
	rep = subtitle.split('-')[-1]

	y_fdr, yerr_fdr,y_pow,yerr_pow, var_range = \
	output[:,:,0],output[:,:,1], output[:,:,2], output[:,:,3], output[:,:,4]  

	variable_in_word = map_to_latex(variable) 

	plot_mountain_vs_valley_single(x=var_range,y=y_fdr,yerr=yerr_fdr,
		title = r'FDR vs. {}'.format(variable_in_word), 
		xlabel = variable_in_word,
		ylabel = r'FDR',#r'Achieved FDR',
		name = 'results/figs/mountain-valley-fdr.eps',  
		 graphtype='FDR') 

	plot_mountain_vs_valley_single(x=var_range,y=y_pow,yerr=yerr_pow,
		title = r'Power vs. {}'.format(variable_in_word),
		xlabel = variable_in_word,
		ylabel = r'Power',
		name = 'results/figs/mountain-valley-power.eps',  
		 graphtype='power') 

def plot_BH_vs_DAGGER_single(x,y,yerr, title, xlabel, ylabel, name,
		  graphtype='power'):
	"""
	This function plots the figure with error bar 

	and save it to local.
	"""
	plt.rc('text', usetex = True)
	plt.rc('font', family = 'serif') 

	fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize = (10,8)) 
	ax.errorbar(x,y[:,0],marker='o',yerr = yerr[:,0], label = 'DAGGER',
		color='red')
	ax.errorbar(x,y[:,1],marker='o',yerr = yerr[:,1], label = 'BH',color='green')

	if graphtype == 'FDR': 
		plt.plot(np.array([x[0],x[-1]]),np.array([0.2,0.2]), 
			 linestyle='dashed',
			 label = 'Targeted FDR',
			color='black')
		
	ax.legend(loc = 'upper right',fontsize = 22)
	ax.set_title(title)
	ax.set_xlabel(xlabel,fontsize = 28)
	ax.set_ylabel(ylabel,fontsize = 28)
	ax.set_ylim(bottom = 0,top=1)
	plt.savefig(name, format='eps', dpi=1000) 	

def plot_BH_vs_DAGGER(path, saved_file):
	output = np.load(path + saved_file)
	subtitle = saved_file.split('_')[0]  

	y_fdr, yerr_fdr,y_pow,yerr_pow, var_range = \
	output[:,:,0],output[:,:,1], output[:,:,2], output[:,:,3], output[:,0,4]  

	variable_in_word = map_to_latex('pi0',simplified = False)
	
	ratio = subtitle[-1]
	plot_BH_vs_DAGGER_single(x=var_range,y=y_fdr,yerr=yerr_fdr,
		title = r'FDR vs. {}, with ratio ={}'.format(variable_in_word, ratio), 
		xlabel = variable_in_word,
		ylabel = r'FDR',#r'Achieved FDR',
		name = 'results/figs/BH-DAGGER-fdr.eps',  
		 graphtype='FDR') 

	plot_BH_vs_DAGGER_single(x=var_range,y=y_pow,yerr=yerr_pow,
		title = r'Power vs. {}, with ratio ={}'.format(variable_in_word, ratio),
		xlabel = variable_in_word,
		ylabel = r'Power',
		name = 'results/figs/BH-DAGGER-power.eps',  
		 graphtype='power') 