# in Igraph

# 1. Clsuter and degree distribution

cluster.distribution(graph)

degree.distribution(graph)

# Histogram of the degree
hist(degree(ga.graph))



# 2. Cherries and tree balance: tree spectrum

spectrum.treeshape(tree) #  sequence containing the number of subtrees of size n, n-1, ..., 3, 2
# where n is the size of the tree. The 'k'th element of the sequence
# is the number of subtrees of size n-k+1 in the tree, where n is the number of tips of the tree.

# fit.stat.slope
# library(apTreeshape)
fit.stat.slope <- function(tree=tree11){

  a <- as.phylo(tree)
  tips.labels <- a$tip.label
  numb.tips <- length(tips.labels)
  d <- rev(seq(from = 2, to = numb.tips, by = 1)) # sequence of sizes of subtrees (number of tips)

  tree1 <- tree# as.treeshape(tree)
  s <- spectrum.treeshape(tree1) # sequence of number of subtrees

  # delete blank values
  nonzero.position = which(s != 0)
  s = s[nonzero.position]
  d = d[nonzero.position]

  reg <- lm(log(s) ~ log(d))
  cozf = coef(reg)

  # we expect the sizes of subtrees to decrease follwoing power-law hypothetically
  power.law.fit = function(x) exp(cozf[[1]] + cozf[[2]] * log(x))
  alpha = -cozf[[2]]
  R.square = summary(reg)$r.squared
  print(paste("Alpha =", round(alpha, 3)))
  print(paste("R square =", round(R.square, 3)))

  # plot
  plot(s ~ d, log = "xy", xlab = "Subtree size (log)", ylab = "Number of subtree (log)",
       col = 1, main = "Degree Distribution")
  curve(power.law.fit, col = "red", add = T, n = length(d))

}

fit.stat.slope(tree=tree11)



# 3. Degree distribution of a graph


# write a function to plot the degree distribution
plot_degree_distribution = function(graph) {
  # calculate degree
  d = degree(graph, mode = "all")
  dd = degree.distribution(graph, mode = "all", cumulative = FALSE)
  degree = 1:max(d)
  probability = dd[-1]
  # delete blank values
  nonzero.position = which(probability != 0)
  probability = probability[nonzero.position]
  degree = degree[nonzero.position]
  # plot
  plot(probability ~ degree, log = "xy", xlab = "Degree (log)", ylab = "Probability (log)",
       col = 1, main = "Degree Distribution")
}

plot_degree_distribution(graph = ga.graph)


# 4. Fitting degree distribution with power-law distribution

# plot and fit the power law distribution
fit_power_law = function(graph) {
  # calculate degree
  d = degree(graph, mode = "all")
  dd = degree.distribution(graph, mode = "all", cumulative = FALSE)
  degree = 1:max(d)
  probability = dd[-1]
  # delete blank values
  nonzero.position = which(probability != 0)
  probability = probability[nonzero.position]
  degree = degree[nonzero.position]
  reg = lm(log(probability) ~ log(degree))
  cozf = coef(reg)
  power.law.fit = function(x) exp(cozf[[1]] + cozf[[2]] * log(x))
  alpha = -cozf[[2]]
  R.square = summary(reg)$r.squared
  print(paste("Alpha =", round(alpha, 3)))
  print(paste("R square =", round(R.square, 3)))
  # plot
  plot(probability ~ degree, log = "xy", xlab = "Degree (log)", ylab = "Probability (log)",
       col = 1, main = "Degree Distribution")
  curve(power.law.fit, col = "red", add = T, n = length(d))
}

fit_power_law(graph = ga.graph)


# Branch length labels for all bracnhes

library(phylotools)

# tree<-pbtree(n=10)
tree <- phylo.tree
# tree$edge.length<-round(tree$edge.length,3)
tree$edge.length<-tree$edge.length # remove rounding
n<-length(tree$tip.label)
ee<-setNames(tree$edge.length[sapply(1:n,function(x,y)
  which(y==x),y=tree$edge[,2])],tree$tip.label)

# plotTree(tree)

plot(tree)
edgelabels(round(tree$edge.length,3),cex=0.6) # add edges labels


# Branch length of the tips

d <- as.data.frame(ee)
branc.leng.tips <- as.data.frame(d$ee) # yes we can find them in phylo.tree$edge.length but here they are not rounded whereas
# with ee are rounded with 3 digit
tips.id <- as.data.frame(as.character(tree$tip.label))

tips.id <- tree$tip.label
tips.branch.leng.raw <- cbind(tips.id,branc.leng.tips)

names(tips.branch.leng.raw) <- c("Ind", "branch.len")

tips.branch.leng <- tips.branch.leng.raw

plot(tree)
edgelabels(round(tips.branch.leng$branch.len,3),cex=0.6) # add edges labels

