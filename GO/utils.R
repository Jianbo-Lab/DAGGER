library(cherry)

load("sets.RData")
library(cherry)
library(igraph)

# dag <- construct(cellprolifgenessets)
# construct a topologically sorted dag from sets:
construct_sorted_dag <- function(set){
	dag = construct(set)
	g = graph_from_adj_list(dag@children, mode='out')
	indices = topo_sort(g,mode='out')
	set = dag@sets[indices]
	dag = construct(set)
	return(dag)
}

find_depth_of_node <- function(lst, depths){
	# This function finds the depth of a certain node given 
	# its parent list for a DAG. 
	# Roots of a DAG has depth 1. 

	# Input:
	# lst: The parent list of the node.

	# depths: A list/array recording the depth of each node. 
	# (With undetermined depth set to 0.)

	# Output:
	# an integer denoting the depth of the node. 
	depth <- if (length(lst)==0) 1 else max(depths[lst]) + 1
}

find_depths <- function(parents_list){
	# Input: 
	# parents_list: a list of lists of parents for each node. 
	# (Ordered in topological order)

	# Output: 
	# An array indicating the depth of each node.
	num <- length(parents_list)
	depths <- rep(0, num)
	for (i in 1:num){
		depths[[i]] <- find_depth_of_node(parents_list[[i]], depths)
	}
	return(depths)
} 


find_leaves <- function(dag){
	# This function finds the leaves of a DAG.
	lst = lapply(dag@children, is.null)
	return(which(as.logical(lst)))
}
# Assign hypotheses to leaves.
# Return a list of booleans with nulls being true,
# alternatives being false and non-leaves being NA. 
assign_hypo_leaves <- function(dag, pi0){
	lst_hypos = rep(NA, length(dag@sets))
	names(lst_hypos) <- names(dag@sets)

	leave_ids <- find_leaves(dag)
	# leaves <- names(lst_hypos)[leave_ids]
	l <- length(leave_ids)
	lst_hypos[leave_ids] <- FALSE

	null_leaves <- sample(leave_ids, size = as.integer(pi0 * l),replace = FALSE)
	lst_hypos[null_leaves] <- TRUE 

	return(lst_hypos)
}
# Assign hypotheses to all nodes.
# Return a list of booleans with nulls being true and alternatives being false. 
assign_hypo_all <- function(dag, pi0){
	parents <- dag@parents
	lst_hypos <- assign_hypo_leaves(dag, pi0)
	l <- length(lst_hypos)
	for (i in l:1){  
		if (length(dag@children[i][[1]])!=0){
			lst_hypos[[i]] <- all(lst_hypos[dag@children[i][[1]]])
		}
	}
	return(lst_hypos)
}

# Compute p-values for each node. 
compute_pvalues <- function(hypo_list, mu){
	x <- rnorm(length(hypo_list), mean = 0, sd = 1)
	means = ifelse(hypo_list, 0, mu)
	x <- x + means 
	pvalues <- 1 - pnorm(x)
	names(pvalues) <- names(hypo_list)
	return(pvalues) 
}
# Compute simes for a single node given its children's pvalues.
compute_simes_single <- function(pvalues){
	pvalues <- sort(pvalues, decreasing = FALSE)
	scales <- length(pvalues) / 1:length(pvalues)
	return(min(pvalues*scales))
}

compute_depth_vary_pvalues <- function(dag, hypo_list, mu,incre=0.3){
	parents <- dag@parents
	depths <- find_depths(parents)

	x <- rnorm(length(hypo_list), mean = 0, sd = 1)
	means = ifelse(hypo_list, rep(0, length(hypo_list)), 
		mu + incre * (max(depths) - depths)) 
	x <- x + means 
	pvalues <- 1 - pnorm(x)
	names(pvalues) <- names(hypo_list)
	return(pvalues) 
}

# Compute simes' p-values for each node. 
compute_simes_pvalues <- function(dag, hypo_list, mu){
	# lst_hypos = rep(NA, length(dag@sets))
	# names(lst_hypos) <- names(dag@sets)

	leave_ids <- find_leaves(dag)
	# leaves <- names(dag@sets)[leave_ids] 

	pvalues <- rep(NA, length(dag@sets))

	names(pvalues) <- names(dag@sets)

	pvalues[leave_ids] <- compute_pvalues(hypo_list[leave_ids], mu)

	l <- length(pvalues)
	for (i in l:1){  
		if (length(dag@children[i][[1]])!=0){
			pvalues[[i]] <- compute_simes_single(pvalues[dag@children[i][[1]]])
		}
	}
	return(pvalues)
}
# Find leaves for each node.
find_leaves_each_node <- function(dag){
	leaves_each_node <- rep(NA, length(dag@sets))
	l <- length(dag@sets)
	for (i in l:1){
		if (length(dag@children[i][[1]])==0){
			leaves_each_node[i] <- list(i)
		} else {
			children <- dag@children[i][[1]]  
			leaves_each_node[i] <- list(unique(unlist(leaves_each_node[children])))
		}
	}
	return(leaves_each_node)
}

compute_simes_pvalues <- function(dag, hypo_list, mu){
	# lst_hypos = rep(NA, length(dag@sets))
	# names(lst_hypos) <- names(dag@sets)
	leave_ids <- find_leaves(dag)
	leaves_each_node <- find_leaves_each_node(dag)
	# leaves <- names(dag@sets)[leave_ids] 

	pvalues <- rep(NA, length(dag@sets))

	names(pvalues) <- names(dag@sets)

	pvalues[leave_ids] <- compute_pvalues(hypo_list[leave_ids], mu)

	l <- length(pvalues)
	for (i in l:1){  
		if (length(dag@children[i][[1]])!=0){
			pvalues[[i]] <- compute_simes_single(pvalues[unlist(leaves_each_node[i])])
		}
	}
	return(pvalues)
}

construct_adj_matrix <- function(dag){
	g <- graph_from_adj_list(dag@children, mode='out')
	adj_matrix <- as_adjacency_matrix(g, 
		names = FALSE, sparse = FALSE) 
	return(adj_matrix)
}
# sets = list(a=c('2','3'),b=c('4','5'),c=c('6','7'),d=c(),e=c(),f=c(),g=c())
# names(sets)=c(1,2,3,4,5,6,7)



