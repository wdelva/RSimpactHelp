#' Compte reconstructed transmission network optimisation criteria

# we choose a network with

# v.link.dens ~ 1
# v.clust << 1
# v.r2 \in 80 - 90
# v.links \in [tips, 2*tips]

optimEpsilon <- function(phylo.tree = epi.tree){

  # look on branching time
  branch.time <- branching.times(phylo.tree)

  # Branch length of the tips
  tree <- phylo.tree
  # tree$edge.length<-round(tree$edge.length,3)
  tree$edge.length<-tree$edge.length # remove rounding
  n<-length(tree$tip.label)
  ee<-setNames(tree$edge.length[sapply(1:n,function(x,y)
    which(y==x),y=tree$edge[,2])],tree$tip.label)

  d <- as.data.frame(ee)
  branc.leng.tips <- as.data.frame(d$ee) # yes we can find them in phylo.tree$edge.length but here they are not rounded whereas
  # with ee are rounded with 3 digit
  tips.id <- as.data.frame(as.character(tree$tip.label))

  tips.id <- tree$tip.label
  tips.branch.leng.raw <- cbind(tips.id,branc.leng.tips)

  names(tips.branch.leng.raw) <- c("Ind", "branch.len")

  tips.branch.leng <- tips.branch.leng.raw

  # return(tips.branch.leng)

  # define a symetric matrix of difference between pairs of tips
  tips.diff.matrix <- matrix(nrow = nrow(tips.branch.leng), ncol = nrow(tips.branch.leng))
  for(i in 1:nrow(tips.branch.leng)){
    for(j in 1:nrow(tips.branch.leng)){
      tips.diff.matrix[i,j] <- abs(tips.branch.leng[,2][i] - tips.branch.leng[,2][j])
    }
  }

  # retrieve elements of the matrix in a vector removing 0

  ele.matrix.raw <- vector()
  for(i in 1:nrow(tips.branch.leng)){
      ele.matrix.raw <- c(ele.matrix.raw,tips.diff.matrix[,i])
  }

  a <- sort(unique(ele.matrix.raw))
  b <- a[-1]

  ele.matrix <- b[1:nrow(tips.branch.leng)]

  source("/home/david/RSimpactHelp/R/Projects_2017/Portugal_David/properties_network.R")
  source("/home/david/RSimpactHelp/R/ConnectNearBy.R")

  v.r2 <- vector()
  v.eps <- vector()
  v.links <- vector()
  v.link.dens <- vector()
  v.clust <- vector()
  for(i in 1: length(ele.matrix)){
    net.elt <- ConnectNearBy(phylo.tree = epi.tree, epsilon = ele.matrix[i])
    properties.net.elt <- properties_network(graph = net.elt)
    Rsquare <- properties.net.elt[[6]]
    links <- properties.net.elt[[2]]
    links.dens <- properties.net.elt[[3]]
    clust.coef <- properties.net.elt[[7]]

    v.eps <- c(v.eps, ele.matrix[[i]])
    v.r2 <- c(v.r2, Rsquare)
    v.links <- c(v.links,links)
    v.link.dens <- c(v.link.dens, links.dens)
    v.clust <- c(v.clust,clust.coef)

  }

  opt.crit <- as.data.frame(cbind(v.eps,v.r2,v.links,v.link.dens,v.clust))

  return(opt.crit)
}
