#### Tree
stepup <- function(pvals, alpha.mat){
    ## alpha.mat: row corresponds to r, column corresponds to i.
    n <- length(pvals)
    psi <- apply(alpha.mat, 1, function(alpha){
        sum(pvals < alpha)
    })
    inds <- which(psi == 1:n)
    if (length(inds) == 0){
        rej <- rep(FALSE, n)
    } else {
        max.ind <- max(inds)
        thresh <- sort(pvals)[max.ind]
        rej <- (pvals <= thresh)
    }
    return(rej)
}

subtree.info <- function(tree, node, max.depth){
    subtree.nodes <- neighborhood(tree, max.depth, node,
                                  mode = "out")[[1]]
    nnodes <- length(subtree.nodes)
    subtree <- induced_subgraph(tree, subtree.nodes)
    internal.nodes <- bfs(subtree, 1, father = TRUE)$father[-1]
    internal.nodes <- unique(internal.nodes)
    nleaves <- nnodes - length(internal.nodes)
    return(list(nnodes = nnodes, nleaves = nleaves))
}

lynch.guo <- function(pvals, tree, root = 1,
                      alpha.list = seq(0.01, 0.3, 0.01),
                      type = 1){
    n <- length(pvals)
    parse.tree <- bfs(tree, root = root, dist = TRUE, father = TRUE)
    depth <- parse.tree$dist
    max.depth <- max(depth)
    parents.list <- parse.tree$father
    info <- lapply(1: n, function(node){
        subtree.info(tree, node, max.depth)
    })
    nleaves.list <- sapply(info, function(x){x$nleaves})
    nnodes.list <- sapply(info, function(x){x$nnodes})    
    total.nleaves <- nleaves.list[1]
    rejs <- sapply(alpha.list, function(alpha){
        rej <- rep(FALSE, n)
        if (pvals[root] > alpha){
            return(rej)
        }
        R <- 1
        rej[1] <- TRUE
        for (d in 1:max.depth){
            candids <- which(depth == d)
            parents <- parents.list[candids]
            if (all(!rej[parents])){
                break
            }
            candids <- candids[rej[parents]]
            m <- length(candids)
            alpha.mat <- sapply(candids, function(candid){
                nleaves <- nleaves.list[candid]
                nnodes <- nnodes.list[candid]
                if (type == 1){
                    alpha.ir <- nleaves / total.nleaves * alpha *
                        (nnodes + R + 0: (m - 1)) / nnodes
                } else if (type == 2){
                    children <- neighborhood(tree, 1, candid,
                                             "out", 1)[[1]]
                    if (length(children) == 0){
                        alpha.ir <- (R + 1:m) * alpha /
                            total.nleaves
                    } else {
                        new.r <- R + 1:m
                        alpha.ir <- nleaves * new.r * alpha /
                            (total.nleaves +
                                 nleaves * (new.r - 1) * alpha)
                    }
                }
                return(alpha.ir)
            })
            tmp.pvals <- pvals[candids]
            tmp.rej <- stepup(tmp.pvals, alpha.mat)
            rej[candids] <- tmp.rej
            R <- R + sum(tmp.rej)
        }
        return(rej)
    })
    return(rejs)
}

#### DAG
DAG.hierarchy <- function(DAG){
    n <- length(V(DAG))
    nodes <- as.numeric(which(degree(DAG, mode = "out") < 1))
    new.nodes <- nodes
    hierarchy <- list()
    k <- 1
    hierarchy[[k]] <- nodes
    while (length(nodes)<n){
        # print(class(new.nodes))
        new.nodes <- neighborhood(DAG, 1, new.nodes, "in",
                                  mindist = 1)
        new.nodes <- as.numeric(unique(unlist(new.nodes)))

        valid <- unlist(sapply(new.nodes, function(j){
            children <- neighborhood(DAG, 1, j, "out",
                                    mindist = 1)
            children <- as.numeric(unlist(children)) 
            all(children %in% nodes)
        }))
        # print(class(valid))
        new.nodes <- new.nodes[valid] 
        if (length(new.nodes)==0) break
        k <- k + 1
        hierarchy[[k]] <- new.nodes
        nodes <- new.nodes
    }
    return(rev(hierarchy))
}



#     nodes <- as.numeric(which(degree(DAG, mode = "in") < 1))
#     new.nodes <- nodes
#     hierarchy <- list()
#     k <- 1
#     hierarchy[[k]] <- nodes
#     while (length(nodes) < n){
#         new.nodes <- neighborhood(DAG, 1, new.nodes, "out",
#                                   mindist = 1)
#         new.nodes <- as.numeric(unique(unlist(new.nodes)))
#         valid <- sapply(new.nodes, function(j){
#             parents <- neighborhood(DAG, 1, j, "in",
#                                     mindist = 1)
#             parents <- as.numeric(unlist(parents))
#             all(parents %in% nodes)
#         })
#         new.nodes <- new.nodes[valid]
#         k <- k + 1
#         hierarchy[[k]] <- new.nodes
#         nodes <- c(nodes, new.nodes)
#     }
#     return(hierarchy)
# }
#dag to depth. 
# effecitive num of leave.

eff.nleaves <- function(DAG){
    hierarchy <- DAG.hierarchy(DAG)
    m <- length(hierarchy)
    n <- length(V(DAG))
    s.mat <- diag(rep(1, n))
    total.ell <- length(hierarchy[[m]])
    if (m == 1){
        ell <- rep(1, n)
        return(list(ell = ell, total.ell = total.ell))
    }
    for (k in 1:(m-1)){
        nodes <- hierarchy[[k+1]]
        upper.nodes <- unlist(hierarchy[1:k])
        for (node in nodes){
            parents <- neighborhood(DAG, 1, node, "in",
                                    mindist = 1)
            parents <- as.numeric(unlist(parents))
            s.mat[upper.nodes, node] <-
                apply(s.mat[upper.nodes, parents, drop = FALSE],
                      1, mean)
        }
    } 
    leaves <- which(degree(DAG, mode = "out") < 1)
    ell <- apply(s.mat[, leaves], 1, sum)
    return(list(ell = ell, total.ell = total.ell,
                hierarchy = hierarchy))    
}
# thresh.
DAG.test <- function(DAG, pvals, thresh){
    n <- length(V(DAG))
    hierarchy <- DAG.hierarchy(DAG)
    m <- length(hierarchy)
    rej <- rep(FALSE, n)
    for (k in 1:m){
        nodes <- hierarchy[[k]] 
        rej.by.thresh <- (pvals[nodes] <= thresh[nodes]) 
        rej.nodes <- nodes[rej.by.thresh]
        rej[rej.nodes] <- TRUE
        if (!all(rej.by.thresh) && k < m){
            accept.nodes <- nodes[!rej.by.thresh]
            accept.children <-
                neighborhood(DAG, m, accept.nodes, "out")
            accept.children <- as.numeric(
                unique(unlist(accept.children)))
            for (r in (k + 1): m){
                hierarchy[[r]] <- setdiff(hierarchy[[r]],
                                          accept.children)
            }
        }
    }
    return(rej)
}
# V(DAG) = 1:n rej TRue:rej.

## SCR <- function(DAG, pvals, base.fun, alpha = 0.05, ...){
##     n <- length(V(DAG))
##     for (r in seq(n, 0, -1)){
##         if (r == 0){
##             WR.inv <- 1
##             break
##         }
##         beta <- alpha * r
##         rej <- base.fun(DAG = DAG, pvals = pvals, beta.list = beta)
##         R <- sum(rej)
##         if (r <= R){
##             WR.inv <- R
##             break
##         }
##     }
##     beta <- alpha * WR.inv
##     rej <- base.fun(DAG = DAG, pvals = pvals, beta = beta)
##     return(rej)
## }


## SCR <- function(DAG, pvals, base.fun,
##                 alpha.list = seq(0.01, 0.3, 0.01), ...){
##     n <- length(V(DAG))
##     m <- length(alpha.list)
##     WR.inv.list <- rep(1, m)
##     beta.list <- as.numeric((1:n) %*% t(alpha.list))
##     rej.list <- base.fun(DAG = DAG, pvals = pvals,
##                          beta.list = beta.list)
##     R.mat <- matrix(apply(rej.list, 2, sum), nrow = n, ncol = m)
##     WR.inv.list <- apply(R.mat, 2, function(R.list){
##         tmp <- which(R.list <= 1:n)
##         WR.inv <- ifelse(length(tmp) == 0, 1, max(tmp))
##         return(WR.inv)
##     })
##     final.beta.list <- alpha.list * WR.inv.list
##     rej <- base.fun(DAG = DAG, pvals = pvals,
##                     beta.list = final.beta.list)
##     return(list(rej = rej, WR.inv = WR.inv))
## }

## SCR <- function(DAG, pvals, base.fun,
##                 alpha.list = seq(0.01, 0.3, 0.01),
##                 print.quiet = FALSE,
##                 ...){
##     n <- length(V(DAG))
##     m <- length(alpha.list)
##     WR.inv.list <- rep(1, m)
##     R.max <- n
##     for (ind in seq(m, 1, -1)){
##         alpha <- alpha.list[ind]
##         beta.list <- alpha * (1:R.max)
##         rej.list <- base.fun(DAG = DAG, pvals = pvals,
##                              beta.list = beta.list,
##                              alpha = alpha, ...)
##         R.list <- apply(rej.list, 2, sum)
##         tmp <- which(R.list >= 1:n)
##         if (length(tmp) == 0){
##             break
##         } else {
##             WR.inv <- max(tmp)
##             WR.inv.list[ind] <- WR.inv
##             R.max <- WR.inv
##         }
##         if (!print.quiet){
##             print(paste0("alpha: ", alpha, ", R: ", WR.inv))
##         }
##     }
##     rej <- sapply(1:m, function(i){
##         alpha <- alpha.list[i]
##         beta <- alpha * WR.inv.list[i]
##         base.fun(DAG = DAG, pvals = pvals,
##                  beta.list = beta,
##                  alpha, ...)
##     })
##     return(list(rej = rej, WR.inv = WR.inv.list))
## }

SCR <- function(DAG, pvals, base.fun,
                alpha.list = seq(0.01, 0.3, 0.01),
                print.quiet = FALSE,
                tolerate = 0.5,
                ...){
    n <- length(V(DAG))
    m <- length(alpha.list)
    WR.inv.list <- rep(1, m)
    R.max <- n
    alpha <- alpha.list[m]
    while (TRUE){
        rej <- base.fun(DAG = DAG, pvals = pvals,
                        beta.list = alpha * R.max,
                        alpha = alpha, ...)
        R <- sum(rej)
        if (R >= R.max * tolerate){
            break
        }
        R.max <- ceiling(R.max / 2)
    }
    for (ind in seq(m, 1, -1)){
        alpha <- alpha.list[ind]
        tmp.rej <- base.fun(DAG = DAG, pvals = pvals,
                            beta.list = alpha * R.max,
                            alpha = alpha, ...)
        tmp.R <- sum(tmp.rej)
        if (tmp.R >= R.max){
            WR.inv <- R.max
            WR.inv.list[ind] <- WR.inv
            if (!print.quiet){
                print(paste0("alpha: ", alpha, ", R: ", WR.inv))
            }
            next
        }
        R.max <- R.max - 1
        beta.list <- alpha * (1:R.max)
        rej.list <- base.fun(DAG = DAG, pvals = pvals,
                             beta.list = beta.list,
                             alpha = alpha, ...)
        R.list <- apply(rej.list, 2, sum)
        tmp <- which(R.list >= 1:R.max)
        if (length(tmp) == 0){
            break
        } else {
            WR.inv <- max(tmp)
            WR.inv.list[ind] <- WR.inv
            R.max <- WR.inv
        }
        if (!print.quiet){
            print(paste0("alpha: ", alpha, ", R: ", WR.inv))
        }
    }
    rej <- sapply(1:m, function(i){
        alpha <- alpha.list[i]
        beta <- alpha * WR.inv.list[i]
        base.fun(DAG = DAG, pvals = pvals,
                 beta.list = beta,
                 alpha, ...)
    })
    return(list(rej = rej, WR.inv = WR.inv.list))
}

SCR.DAG.base <- function(DAG, pvals, beta.list, alpha,
                         lambda.c = 2,
                         info = NULL){
    if (is.null(info)){ 
        info <- eff.nleaves(DAG) 
    }
    ell <- info$ell
    L <- info$total.ell 
    lambda <- lambda.c * alpha
    alpha.base <- ell / (1 + ell * lambda)
    alpha.factor <- pmin(beta.list / L, lambda)
    thresh.list <- alpha.base %*% t(alpha.factor)

    rej.list <- apply(thresh.list, 2, function(thresh){
        DAG.test(DAG, pvals, thresh) 

    })
    return(rej.list)
}

# lambda.c=2alpha
SCR.DAG <- function(DAG, pvals, lambda.c = 2,
                    alpha.list = seq(0.01, 0.3, 0.01)){
    info <- eff.nleaves(DAG)     
    base.fun <- function(DAG, pvals, beta.list, alpha){
        SCR.DAG.base(DAG, pvals, beta.list, alpha, lambda.c, info)
    }
    result <- SCR(DAG, pvals, base.fun, alpha.list)
    return(result)
}

BH.DAG.base <- function(DAG, pvals, beta.list, alpha){
    n <- length(V(DAG))
    rej.list <- sapply(beta.list, function(beta){
        thresh <- rep(beta / n, n)
        DAG.test(DAG, pvals, thresh)
    })
    return(rej.list)
}

BH.DAG <- function(DAG, pvals, alpha.list = 0.05){
    result <- SCR(DAG, pvals, BH.DAG.base, alpha.list)
    return(result)
}

