#' Simulate the sequences for individual in the transmissiion network
#'
#' From the \code{\link{transnetwork2tree()}} we have the object of the transmission tree.
#' We then can simulate the consesus sequences of individuals using GTR+Gamma model
#' @param transtree Tree (or phylo) object from  \code{\link{trans.network2tree()}}
#' which respresents the transmisssion tree built from transmission networks
#' \code{\link{transmNetworkBuilder.baseline()}} or \code{\link{transmNetworkBuilder.diff()}}
#' @param seedSeq Input sequence of the seed, it might be a string
#' @param alpha Shape parameter of Gamma distribution for rates variation
#' @param r A vector of bases rate of substitutions
#' @param f A list of bases frequences
#' @return Consensus sequences of recipients for each seed
#' @example
#' tree0 <- trans.network2tree(transm.ls[[1]]) # tree from the transmission network of seed 1
#' freq <- c(0.3353293,0.2035928,0.2628077,0.1982701) # bases frequencies
#' rate <- list("a"=0.2, "b"=0.6, "c"=0.12,"d"=0.001, "e"=0.25, "f"=0.24) # substitution rates
#' sim <- sequence.simulation(transtree = tree0, seedSeq = hivSeq, alpha = 0.90,
#' rate.list = rate, base.freq = freq)
#' saveAlignment.PhyloSim(sim, file = paste("file_name",sep="")) alignment saved with internal nodes
#' saveAlignment.PhyloSim(sim, file = paste("file_name",sep=""), skip.internal = TRUE, paranoid = TRUE) alignment saved without internal nodes

#' @import phylosim
sequence.simulation <- function(transtree = tree0, seedSeq = hivSeq,
                                base.freq = freq){

  # define the substitution processes >> package phylosim
  # proc <- GTR(rate.params = rate.list,
  #             base.freqs = base.freq)

  proc <- F81(base.freqs = base.freq)

  # attach process to the nucleotides sequence >> package phylosim
  nucleproc <- NucleotideSequence(string = seedSeq, processes = list(list(proc)))

  # plusGamma(nucleproc,proc,alpha)

  # simulate the sequences >> package phylosim
  simsequence <- Simulate(PhyloSim(root.seq = nucleproc,
                                   phylo = transtree))

  # return the alignement
  return(simsequence)
}
