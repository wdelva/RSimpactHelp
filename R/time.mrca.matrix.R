#' Matrix of time to the most recent common ancestor (MRCA)
#'
#' Coumpte time point to the  most recent common ancestor (MRCA) bewteen individuals
#'
#' @param tree A tree (phylo object) in nexus format (a phylogeny)
#' @return A symmetrix matrix where each entry is time where two sequences have the same most recent ancestor
#' @example
#' time.mrca <- time.mrca.matrix(tree = phylo.tree)
#' @import ape

time.mrca.matrix <- function(tree = phylo.tree){

    Mytree <- tree


    # Use of mrca() and branching.time() functions from ape package

    # Symmetric matrix with tips which where entries are the internal nodes which represents the MRCA between two tips
    # from the phylogenetic tree

    Mytree.mrca.tips <- as.data.frame(mrca(Mytree, full = FALSE)) # MRCA for tips only >> full = FALSE

    # The assignment is to use the internal nodes and replace them by branching time


    # tips which represent individuals (sequences)
    inds <- Mytree$tip.label

    # branching time
    branch.time <- branching.times(Mytree) # each internal_node is associated with a branching time

    # Preparing the matrix of time to the MRCA
    ##########################################
    # Change the Mytree.mrca.tips matrix a bit to be easily handle for as igraph or network obejct
    # we have to make the diagonal 0
    n = as.matrix(Mytree.mrca.tips)
    mat.fun.tree <- function(n){
        diag(n) <- 0
        m <- n
        return(m)
    }

    m <- as.data.frame(mat.fun.tree(n))

    # Get internal nodes of the phyloegentic tree

    i.nodes <- NULL
    int.nodes <- function(m){
        for(i in 1:nrow(m)){
            for(j in 1:nrow(m)){
                i.nodes <- c(i.nodes,m[i,j])
            }
        }
        return(unique(i.nodes))
    }

    h = int.nodes(m)
    # remove the first element which is the 0 on diagonal
    g <- h[(2:length(h))] # these are really the internal nodes which represents the MRCA in the phylogenetic tree

    # MRCA times
    # each internal node has its branching time
    # thus, we have to replace the internal nodes in the matrix with their branching times respectively

    mrca.time.val <- as.data.frame(matrix(,nrow(m),nrow(m))) # initiate an n*n empty matrix (n tips)

    mrca.time.matrix <- function(g,m,branch.time){ # write a function which take as parms the internal nodes (g)
                                                   # the matrix of internal nodes (MRCA) which is m and
                                                   # the branching time which correspond with time where the two seq eveoled separelty
        k <- as.data.frame(branch.time)
        p <- k$branch.time

        for(l in 1:length(g)){
            for(i in 1:nrow(m)){
                for(j in 1:nrow(m)){
                    if(m[i,j] == g[l]){
                        mrca.time.val[i,j] <- p[l]
                    }
                }
            }

        }
        return(mrca.time.val)
    }

    v = mrca.time.matrix(g,m,branch.time) # the diagonal is empty with entries NA

    # remove the NA in the diagonal and put 0 since no MRCA of a given sequence comparing at itself
    matrix.mrca.time.func <- function(v){
        diag(v) <- 0
        times.mrca.fine <- v
        return(times.mrca.fine)
    }

    names.inds <- names(m) # get the names of the tips
    mrca.times.final <- matrix.mrca.time.func(v)
    names(mrca.times.final) <- names.inds # put the names of the tips on the columns

    mrca.times.done <- mrca.times.final

    return(mrca.times.done)
}



