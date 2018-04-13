#' Build a transmission tree
#'
#' @param transnetwork The phylo object produced by \code{\link{transmission.network.builder}}.
#' @return Transmission tree object (a phylo object)
#' @examples
#' transm.ls <- transmission.network.builder(datalist = datalist, endpoint = 40)
#' transnetwork <- transm.ls[[16]]
#' tree0 <- trans.network2tree(transnetwork = transnetwork) # transmission tree of seeds 16
#'
#' @importFrom expoTree epi2tree
#' @export
#'
trans.network2tree <- function(transnetwork = transnetwork){

  transm.phylo <- epi2tree2(transnetwork)

  return(transm.phylo)

}
