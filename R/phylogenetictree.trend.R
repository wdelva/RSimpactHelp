#' Phylogenetic tree spectum
#'
#' Compute two vectors, the number of tips for a given subtree (or size of the tree)
#' and the number of trees (subtrees ) with corresponding tips for a given phylogenetic tree
#' This tree decomposition informthe general trend of the transmission and can be sometimes folloging
#' the pwer law \code{\link{powerLaw.tree}}
#'
#' @param tree A phylogenetic tree object
#' @return A list with two vector with same length, which is used for regression analysis
#' between the number of subtrees and their sizes.
#'
#' @examples tree <- rtree(23)
#'
#'examp1 <- phylogenetictree.trend(tree = tree)
#' x = examp1$num.tree
#' y = examp1$size.tree
#' reg <- lm(log(y) ~ log(x))
#' cozf = coef(reg)
#'
#' we expect the sizes of subtrees to decrease follwoing power-law hypothetically
#' power.law.fit = function(x) exp(cozf[[1]] + cozf[[2]] * log(x))
#' alpha = -cozf[[2]]
#' R.square = summary(reg)$r.squared
#' print(paste("Alpha =", round(alpha, 3)))
#' print(paste("R square =", round(R.square, 3)))
#'
#' plot(y ~ x, log = "xy", xlab = "Number of subtree (log)", ylab = "Subtree size (log)",
#'      col = 1, main = "Degree Distribution")
#' curve(power.law.fit, col = "red", add = T, n = length(y))
#'
#'@importFrom ape as.phylo
#'@importFrom apTreeshape as.treeshape
#'@importFrom apTreeshape spectrum.treeshape
#'
phylogenetictree.trend <- function(tree=tree){

  a <- as.phylo(tree)
  tips.labels <- a$tip.label
  numb.tips <- length(tips.labels)
  d <- rev(seq(from = 2, to = numb.tips, by = 1)) # sequence of sizes of subtrees (number of tips)

  tree1 <- as.treeshape(tree) # the tre must be a treeshape object
  s <- spectrum.treeshape(tree1) # sequence of number of subtrees

  # in d we read the number of tips for a given tree (or size of the tree)
  # in s we read the number of trees (subtrees ) with corresponding (in d) tips

  # delete blank values
  nonzero.position = which(s != 0)
  s = s[nonzero.position]
  d = d[nonzero.position]

  #num.tree.size <- vector("list",2)

  num.tree.size <- list()

  num.tree.size$num.tree <- s
  num.tree.size$size.tree <- d

  return(num.tree.size)

}

