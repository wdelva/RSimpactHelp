#' Power law fitting of subtrees and their sizes by linear regression
#'
#' Fitting a tree spectrum (number of subtrees and their sizes from one tree) with power law
#' The power law fitting is done by linear regression between number of subtrees and their sizes
#' which can also be computed by \code{\link{phylogenetic.trend}}
#'
#'@param tree A phylogenetic tree object
#'@return Power law fitting coefficents and the R2
#'@example fit.coef.r2 <- powerLaw.tree(tree = tree)
#' @importFrom apTreeshape as.treeshape

powerLaw.tree <- function(tree = tree) {

  # Subtrees

  tree1 <- as.treeshape(tree)
  s <- spectrum.treeshape(tree1) # sequence of number of subtrees

  # Sizes of subtrees
  a <- as.phylo(tree)
  tips.labels <- a$tip.label
  numb.tips <- length(tips.labels)
  d <- rev(seq(from = 2, to = numb.tips, by = 1)) # sequence of sizes of subtrees (number of tips)



  # delete blank values: subtrees which don't have given zises
  nonzero.position = which(s != 0)
  s = s[nonzero.position]
  d = d[nonzero.position]

  reg <- lm(log(s) ~ log(d))
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
