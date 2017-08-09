# phylogenetic tree characteristics

properties_tree <- function(tree=epi.tree){

  library(phangorn)
  library(ape)
  library(apTreeshape)

  epi.tree.shape <- as.treeshape(tree)

  ntips <- length(tree$tip.label)

  colless.indice <- colless(epi.tree.shape) # shape statistic
  sackin.indice <- sackin(epi.tree.shape)

  tree.prop <- as.data.frame(cbind(ntips, colless.indice, sackin.indice))
  return(tree.prop)

}
