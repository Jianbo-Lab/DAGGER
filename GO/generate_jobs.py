import itertools
import numpy as np 
import os
import random
set_names = ['full','cellprolif']#['cellprolif','full']
seeds = range(1,11)
pvalues_types = ['ind','simes']
algorithms = ['DAGGER','lord','BH','all-goeman','any-goeman','SCR.DAG','BH.DAG','Holm']
pi0s = np.arange(0.15,1.00,0.05)

paras_choices = list(itertools.product(set_names, seeds, 
	pvalues_types, algorithms, pi0s))

jobs_list = []
for paras in paras_choices:
	set_name, seed, pvalues_type, algorithm, pi0 = paras 
	key_words = ['python','run_simulations.py',
	'--set_name', set_name, 
	'--seed', str(seed),
	'--pvalues_type', pvalues_type,
	'--algorithm',algorithm,
	'--pi0',str(pi0),
	'\n'
	]
	job = ' '.join(key_words)
	jobs_list.append(job)

random.shuffle(jobs_list)

with open('jobs.sh','wb') as f:
	for job in jobs_list:
		f.write(job)

