#' Create a csv file with parameter combinations for the first wave of Simpact simulations
#'
#' A short description here...
#'
#' @param file.name Where the csv file will be stored, and its name
#' @param n.experiments Number of parameter combinations to be created
#' @param limits list of prior distributions, one per parameter
#' @return a csv file, written to disk
#' @importFrom randtoolbox sobol
#' @export

wave1.creator <- function(file.name = "/user/data/gent/vsc400/vsc40070/agemixing/scripts/wave1.csv",
                          n.experiments,
                          limits){ # e.g.: list(prior.x.1 = c(0,1), prior.x.2 = c(2,4), prior.x.3 = c(-1,1))
  lls <- as.numeric(as.data.frame(limits)[1, ])
  uls <- as.numeric(as.data.frame(limits)[2, ])
  range.width <- uls - lls
  ll.mat <- matrix(rep(lls, n.experiments), nrow = n.experiments, byrow = TRUE)
  range.width.mat <- matrix(rep(range.width, n.experiments), nrow = n.experiments, byrow = TRUE)
  sobol.seq.0.1 <- sobol(n = n.experiments, dim = length(lls), init = TRUE, scrambling = 1, seed = 1, normal = FALSE)
  experiments <- cbind(1:n.experiments, ll.mat + sobol.seq.0.1 * range.width.mat)
  experiments.col.names <- c("seed", paste0("x", seq(1:length(lls))))
  write.table(experiments,
              sep = ",",
              dec = ".",
              quote = FALSE,
            row.names = FALSE,
            col.names = experiments.col.names,
            file = file.name)
}
