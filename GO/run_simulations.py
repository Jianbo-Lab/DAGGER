import subprocess
import argparse
def main(algorithm, set_name, pvalues_type,seed, pi0):
	program = 'Rscript' if algorithm in ['BH','all-goeman','any-goeman','SCR.DAG','BH.DAG','Holm'] else 'python'
	if program == 'Rscript':
		script = 'competitors.R'
		cmd = [program, script, '--method', algorithm,'--set_name',set_name, '--pvalues_type', pvalues_type, '--seed', seed, '--pi0', pi0]
	elif program == 'python':
		script = 'lord.py' if algorithm == 'lord' else 'DAGGER.py'
		cmd = [program, script, '--set_name',set_name, '--pvalues_type', pvalues_type, '--seed', seed, '--pi0', pi0]

	cmd = ' '.join([str(i) for i in cmd]) 
	subprocess.call(cmd, shell=True)

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument("--algorithm", type = str, default = 'DAGGER',
		choices = ['DAGGER','lord','BH','all-goeman','any-goeman','SCR.DAG','BH.DAG','Holm'])
	parser.add_argument("--set_name", type = str, default = 'cellprolif',
		choices = ['cellprolif','full'],
		help="which set of genes are being considered.")
	parser.add_argument("--pvalues_type", type = str, choices = ['ind','simes'],
		default = 'ind', help = 'whether pvalues are independent or simes.')

	parser.add_argument('--seed', type=int, 
	  default=1, help = 'random seed') 

	parser.add_argument('--pi0', type=float, 
	  default=0.95, help = 'proportion of nulls.')

	a = parser.parse_args()

	main(a.algorithm, a.set_name, a.pvalues_type, a.seed, a.pi0) 