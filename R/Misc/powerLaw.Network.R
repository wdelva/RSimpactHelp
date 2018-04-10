#' Power law fitting of network
#'
#' Fitting degree distribution of a network (transmission or relationship) by the power law  distribution
#' The power law fitting is done by linear regression between nodes and their degrees
#'
#'@param network A transmission or relationship network (igraph object)
#'@return Power law fitting coefficents and the R2
#'@example fit.net <- powerLaw.tree(network = trans.network)
#'@import igraph

powerLaw.tree <- function(network = trans.network) {

  graph <- network

  # calculate degree
  d = degree(graph, mode = "all") # degrees of nodes
  dd = degree.distribution(graph, mode = "all", cumulative = FALSE) # frequencies of nodes with degree 0 1 2 ... max(d)
  degree = 1:max(d)
  probability = dd[-1]
  # delete blank values
  nonzero.position = which(probability != 0)
  probability = probability[nonzero.position]
  degree = degree[nonzero.position]

  reg <- lm(log(degree) ~ log(probability))
  cozf = coef(reg)

  # we expect the sizes of subtrees to decrease follwoing power-law hypothetically
  power.law.fit = function(x) exp(cozf[[1]] + cozf[[2]] * log(x))
  alpha = -cozf[[2]]
  beta = -cozf[[1]]
  R.square = summary(reg)$r.squared

  coef.r2 <- list()
  coef.r2$intercept <- round(alpha, 4)
  coef.r2$slope <- round(beta, 4)
  coef.r2$RSquare <- round(R.square, 4)

  return(coef.r2)
}
