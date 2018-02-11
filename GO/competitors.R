library(cherry) # version 05-10
library(globaltest) # version 5.18.0 
load("sets.RData")
library(argparse)
source('utils.R')
print('Packages loaded')
source('Lynch.R') 
#returns rejected vector
Holm <- function(pvalues,alpha)
{
  p.adjust(pvalues,method="holm")<=alpha
}

#returns rejected vector
BH <- function(pvalues,alpha)
{
  p.adjust(pvalues,method="BH")<=alpha
}

#function that carries TRUE upwards (finds true number of false hypotheses, or makes more rejections)
propagate <- function(parents,children,rejected)
{
  n <- length(rejected)
  
  #initialize queue
  queue <- rep(0, n)
  head <- tail <- 1
  
  # reject parents
  # fill queue with nodes that were rejected
  for(i in 1:n)
  {
    if(rejected[i])
    {
      #put node on queue
      queue[tail] <- i
      tail <- tail + 1
    }
  }
  
  #put parents on queue
  while(head < tail)
  {
    current <- queue[head]
    head <- head + 1
    
    #put unrejected parents on queue
    for(j in parents[[current]])
    {
      if(!rejected[j])
      {
        queue[tail] <- j
        tail <- tail + 1
        
        #reject -> only on queue once
        rejected[j] <- TRUE
      }
    }
  }

  return(rejected)
  
}

#use cell proliferation sets with added genes







compute_fdp<- function(rejections,correct_rejections, hypo_list){
  FDP = ifelse(rejections == 0,0,(rejections - correct_rejections) / rejections)
  fdp = mean(FDP) 
  return(fdp)
}

compute_power<- function(rejections,correct_rejections, hypo_list){ 
  POWER = correct_rejections / sum(!hypo_list) 
  power = mean(POWER)
  return(power)
}


run_and_record_single_experiment <- function(set_name, alpha, 
  mu, pi0, pvalues_type, method, incre, seed){
  set.seed(seed) 
  dag <- readRDS(sprintf("data/dag_%s.RData", set_name))

  hypo_list <- read.table(sprintf("data/hypo_%s_%s_%i_%.2f.txt",set_name,pvalues_type, 
  seed, pi0))[,1]#[as.integer(pi0*20),] 

  pvalues <- read.table(sprintf("data/pvalues_%s_%s_%i_%.2f.txt", 
    set_name, pvalues_type, seed, pi0), sep = ' ')[,1]



  #[as.integer(pi0*20),]

  hypo_list <- as.integer(hypo_list)
  pvalues <- as.numeric(pvalues)
  # names(pvalues) <- as.character(dag@sets)

  parents <- dag@parents
  children <- dag@children
  toReject <- !hypo_list
  if (method == 'SCR.DAG' | method == 'BH.DAG'){
    DAG <- readRDS(sprintf("data/igraph_%s.RData", set_name))
  }
  start_time = proc.time()
  if (method == 'Holm'){
      rejected <- Holm(pvalues,alpha) #nodes rejected by Holm     
    } else if (method == 'BH'){
      rejected <- BH(pvalues,alpha)
    } else if (method == 'all-goeman'){
      all <- DAGmethod(dag, alpha_max = alpha, method = "all", 
      pvalues=pvalues)
      rejected <- all@rejected
    } else if (method == 'any-goeman'){
      any = DAGmethod(dag, alpha_max = alpha, 
        method = "any", pvalues=pvalues)
      rejected <- any@rejected 
    } else if (method == 'SCR.DAG'){ 
      rejected <- SCR.DAG(DAG, pvalues,
                    alpha.list = c(alpha))$'rej'[,1] 
      # write.table(rejected, 
      #   file='tmp/test_rejections', 
      #   row.names=FALSE, col.names=FALSE)
      # write.table(toReject, 
      #   file='tmp/test_alters', 
      #   row.names=FALSE, col.names=FALSE)

    } else if (method == 'BH.DAG'){ 
      rejected <- BH.DAG(DAG, pvalues, alpha.list = c(alpha))$'rej'[,1] 
    } else if (method == 'focus_level'){
      rejected <- focus_level(dag, alpha, pvalues)
    }
  # print(proc.time())
  # print(start_time)
  time = as.list(proc.time() - start_time)$'user.self'  
  rej <- sum(rejected)
  correctrej = sum(rejected & toReject) 

  fdp <- compute_fdp(rej, correctrej, hypo_list)
  power <- compute_power(rej, correctrej, hypo_list)
  
  write.table(c(power, fdp, time), 
    file=sprintf("output/%s_%.2f_%s_%s_%i.txt", set_name, pi0, 
      pvalues_type, method, seed), 
    row.names=FALSE, col.names=FALSE)
}



parser <- ArgumentParser()
parser$add_argument('--set_name', type="character", 
  default='cellprolif')
parser$add_argument('--pvalues_type', type="character", 
  default='ind')
parser$add_argument('--method', type="character", 
  default='SCR.DAG')
parser$add_argument('--mu', type="double", 
  default=1.0)
parser$add_argument('--seed', type="integer", 
  default=1)
parser$add_argument('--incre', type="double", 
  default=0.3)

parser$add_argument('--pi0', type="double", 
  default=0.95)

args <- parser$parse_args() 
# methods <-  c('Holm','BH','all-goeman','any-goeman','SCR.DAG','BH.DAG')
# methods <-  c('Holm','BH','all-goeman','any-goeman')
# methods <-  c('SCR.DAG','BH.DAG')
# methods <-  c('SCR.DAG')
# methods <- c('focus_level')
# methods <-  c('BH.DAG') 
# main(set_name = args$set_name, mu = args$mu, incre = args$incre,methods = methods, num = 3) 

run_and_record_single_experiment(args$set_name, 0.2, 
  args$mu, args$pi0, args$pvalues_type, args$method, 
  args$incre, args$seed)

