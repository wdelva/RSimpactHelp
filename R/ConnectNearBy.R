#' Construct a transmission network from a phylogenetic tree (phylo object)
#' Using the branch lengths of leaves (tips), join leaves (tips) with approximately equal branch lengths
#'
#'
#' @param phylo.tree a phylo object
#' @param epsilon cut-off value as difference between two branch lengths
#' @return a transmission matrix

#' @examples
#' transm.mat <- ConnectNearBy(phylo.tree = tree)

#' @import igraph
#' @import ape
#' @import phytools

ConnectNearBy <- function(phylo.tree = tree, epsilon=0.1/2){

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

  ## Algorithm: SEARCH YOUR LINKS


  Ind.g <- vector()
  Ind.g.1 <- vector()
  Ind.g.2 <- vector()
  branch.len.g <- vector()
  branch.len.g.1 <- vector()
  branch.len.g.2 <- vector()
  # | tips.branch.leng$branch.len[i] - tips.branch.leng$branch.len[j] <= epsilon
  for(i in 1:nrow(tips.branch.leng)){
    for(j in 1:nrow(tips.branch.leng)){
      if(tips.branch.leng$branch.len[i]==tips.branch.leng$branch.len[j] | (abs(tips.branch.leng$branch.len[i] - tips.branch.leng$branch.len[j]) <= epsilon) ){
        Ind.g <- c(Ind.g,tips.id[i],tips.id[j]) # there are
        branch.len.g <- c(branch.len.g,tips.branch.leng$branch.len[i],tips.branch.leng$branch.len[j])
        reconst.grap.raw <- as.data.frame(cbind(Ind.g,branch.len.g))

        Ind.g.1 <- c(Ind.g.1, tips.id[i])
        Ind.g.2 <- c(Ind.g.2, tips.id[j])
        branch.len.g.1 <- c(branch.len.g.1, tips.branch.leng$branch.len[i])
        branch.len.g.2 <- c(branch.len.g.2, tips.branch.leng$branch.len[j])

        graph.trans <- as.data.frame(cbind(Ind.g.1,Ind.g.2,branch.len.g.1,branch.len.g.2))

      }
    }
  }

  gaga = as.matrix(graph.trans)
  f = graph.edgelist(gaga[,1:2])

  g <- simplify(f,remove.loops = TRUE) # remove loops

  l <- as.undirected(g, mode = "collapse") # remove bi-directions
  #
  # V(g)
  # E(g)
  return(l)
}
