#' Compute phylogenetic error
#'

phylo_error <- function(tree=treephylo){
  a <- length(tree$edge.length)-length(tree$tip.label)
  return(a)

}
