#' Construct phylogenetic tree with Phangorn and calibrate same tree with treedater from sequence alignment
#'
#' @param simseqfile File containing sequences
#' @return A time-stamped phylogenetic  tree with annotation of internal nodes dates
#' @import ape
#' @import phangorn
#' @export


phylogenetic.tree.phangorn.par <- function(simseqfile = sequ.dna){
  # Construct phylogenetic tree
  # Many tools to build the trees: within R like ape, phangorn, and outside compiled tools like iq-tree and FastTree

  simseq.df <- as.phyDat(simseqfile)
  dm <- dist.ml(simseq.df)
  treeNJ <- NJ(dm)
  # A first naive fit
  fit <- pml(treeNJ, data = simseq.df)
  fitGTRGI <- update(fit, k = 4, inv = 0.35)
  # Now estimating the parameters
  fitGTRGI <- optim.pml(fitGTRGI, model="GTR", optInv=TRUE, optGamma=TRUE,
                        rearrangement = "NNI", control = pml.control(trace = 0))

  fitGTRGI.top <- optim.pml(fitGTRGI, model="GTR", optNni=TRUE, optEdge=TRUE,
                        rearrangement = "NNI", control = pml.control(trace = 0))

  return(fitGTRGI.top)
}
