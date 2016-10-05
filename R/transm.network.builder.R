#' Build a list with HIV transmission network data.
#'
#' Calculate the HIV prevalence at a point in time, for specific age groups and gender strata.
#'
#' @param datalist The datalist that is produced by \code{\link{readthedata()}}
#' @param endtime Only transmission events that took place before this point in simulation time, are captured in the output.
#' @return a list with the transmission network data, as required by the epi2tree function.
#' @examples
#' transm.ls <- transm.network.builder(datalist = datalist, endpoint = 30)

transm.network.builder <- function(datalist = datalist,
                                   endpoint = 30){
  TE = datalist$etable[eventname=="transmission"] # TE = transmission events
  Seed <- datalist$ptable[InfectTime == as.numeric(datalist$itable$hivseed.time)]
  gender_prefix <- c("man_", "woman_")
  SeedID <- paste0(gender_prefix[Seed$Gender + 1], Seed$ID)

  transm.ls <-list()
  seed.time <- as.numeric(datalist$itable$hivseed.time)
  all.itimes <- endpoint - c(seed.time, TE$eventtime)
  transm.ls$itimes <- all.itimes[all.itimes > 0]
  transm.ls$dtimes <- rep(0, length(transm.ls$itimes)) # A placeholder for now.

  all.id.original <- c(SeedID, TE$p2name)
  all.parent.original <- c("root", TE$p1name)
  all.id.integer <- as.integer(1:length(all.id.original) - 1)

  parent.indices <- match(all.parent.original, all.id.original)
  all.parent.integer <- all.id.integer[parent.indices]
  all.parent.integer[1] <- as.integer(-1) # Setting the parent of the seed ID as -1

  transm.ls$id <- all.id.integer[all.itimes > 0]
  transm.ls$parent <- all.parent.integer[all.itimes > 0]
  return(transm.ls)
}
