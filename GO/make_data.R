library(cherry)

load("sets.RData")
library(cherry)
library(igraph)
library(argparse)
source('utils.R')

write_to_table <- function(dag, set_name, pi0s, mu, seed = 0, 
	incre = 0.3, pvalues_type = 'both'){
	print(sprintf('Generate Data for the %i th experiment...', seed))
	set.seed(seed)  

	for (i in 1:length(pi0s)){
		pi0 <- pi0s[i]
		# print(length(compute_depth_vary_pvalues(dag, hypo_list, mu, incre=incre)))

		hypo_list <- assign_hypo_all(dag, pi0)
		pi0 <- sprintf('%0.2f',pi0) 
		if (pvalues_type == 'both'|pvalues_type == 'ind'){
			pvalues_ind <- compute_depth_vary_pvalues(dag, hypo_list, mu, incre=incre)
			write.table(pvalues_ind, file=sprintf("data/pvalues_%s_ind_%i_%s.txt",set_name, seed, pi0), 
	      row.names=FALSE, col.names=FALSE)
			write.table(hypo_list, file=sprintf("data/hypo_%s_ind_%i_%s.txt",set_name, seed, pi0), 
	      row.names=FALSE, col.names=FALSE)

		}
		if (pvalues_type == 'both'|pvalues_type == 'simes'){
			pvalues_simes <- compute_simes_pvalues(dag, hypo_list, mu+1)		
			write.table(pvalues_simes, file=sprintf("data/pvalues_%s_simes_%i_%s.txt",set_name, seed, pi0), 
	      row.names=FALSE, col.names=FALSE)
			write.table(hypo_list, file=sprintf("data/hypo_%s_simes_%i_%s.txt",set_name, seed, pi0), 
	      row.names=FALSE, col.names=FALSE)

		}

		
	}
	# print(pvalues_ind)
}

create_data <- function(set_name = 'cellprolif', mu = 2, incre = 0.3, 
	pvalues_type = 'both', n_replications = 10){
	if (sprintf("dag_%s.RData", set_name) %in% list.files('./data')){
		dag <- readRDS(sprintf("data/dag_%s.RData", set_name))
		print('DAG loaded from previous storage') 
		} 
	else {
		set <- if (set_name == 'cellprolif') cellprolifgenessets else allsets
		dag <- construct_sorted_dag(set)
		adj_matrix <- construct_adj_matrix(dag)

		print('Constructing its igraph correspondence...')
		DAG <- graph_from_adj_list(dag@children, mode='out')

		saveRDS(DAG, sprintf("data/igraph_%s.RData", set_name)) 
		saveRDS(dag, sprintf("data/dag_%s.RData", set_name))
		write.table(adj_matrix, file = sprintf("data/adjmatrix_%s.txt", set_name), 
	      row.names=FALSE, col.names = FALSE)
		print('New DAG created')
		}

	for (i in 1:n_replications){ 
		write_to_table(dag, set_name, seq(0.15,0.95,by = 0.05), mu, i, incre, pvalues_type = pvalues_type) 
	}
}


parser <- ArgumentParser()
parser$add_argument('--pvalues_type', type="character", default='ind') 
parser$add_argument('--set_name', type="character", default='full') 
parser$add_argument('--mu', type="double", default=1.0)#1.0)# 1.0)
parser$add_argument('--incre', type="double", default=0.3)#0.25)#0.3)
parser$add_argument('--n_replications', type="integer", default=10)#0.25)#0.3)

args <- parser$parse_args() 

if (!('data' %in% list.files('./'))) dir.create('data')
if (!('output' %in% list.files('./'))) dir.create('output')
if (!('results' %in% list.files('./'))) dir.create('results')

create_data(set_name = args$set_name, mu = args$mu, 
	incre = args$incre,pvalues_type = args$pvalues_type,
	n_replications = args$n_replications) 


