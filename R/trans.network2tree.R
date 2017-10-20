#' Build a transmission tree
#'
#' From the output of the \code{\link{transmNetworkBuilder.baseline()}} or \code{\link{transmNetworkBuilder.diff()}},
#' we can build the tranmission tree
#' @param transnetwork The phylo object produced by \code{\link{transmNetworkBuilder.baseline}} or \code{\link{transmNetworkBuilder.diff()}}.
#' @return Transmission tree object (a phylo object)
#' @examples
#'
#' transm.ls <- transmNetworkBuilder.baseline(datalist = datalist, endpoint = 40)
#' transnetwork <- transm.ls[[16]]
#' tree0 <- trans.network2tree(transnetwork = transnetwork) # transmission tree of seeds 16
#'
#' @importFrom expoTree epi2tree
#'

trans.network2tree <- function(transnetwork = transnetwork){

  source("/home/david/RSimpactHelp/R/epi2tree2.R")

  transm.phylo <- epi2tree2(transnetwork)

  return(transm.phylo)

}
