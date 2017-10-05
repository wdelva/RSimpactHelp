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

library(phytools)

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




## Tree imbalance: different kinds of transmission networks


library(expoTree)
library(ape)
library(apTreeshape)
library(phangorn)


# Homogeneous transmission: one > two or three

# one infection two infections after
id1 = c(seq(from=0,to=22, by=1)) # 39
par1 = c(-1,0,0,1,1,2,2,3,3,4,4,5,5,6,6,
         7,7,8,8,9,9,10,10) #,11,11,12,12,
# 13,13,14,14,15,15,16,16,17,17,18,18,19)
t1 = rev(c(seq(from=1,to=23, by=1)))
epi1 <- list()
epi1$itimes <- t1
epi1$dtimes <- rep(0, length(t1))
epi1$id <- id1
epi1$parent <- par1
epitree1 <- epi2tree(epi1)
source("~/RSimpactHelp/R/phylogenetictree.trend.R")

source("C:/Users/niyukuri/Documents/New folder/RSimpactHelp/R/phylogenetictree.trend.R")

examp1 <- phylogenetictree.trend(tree = epitree1)
x1 = examp1$num.tree
y1 = examp1$size.tree
reg1 <- lm(log(y1) ~ log(x1))
cozf1 = coef(reg1)

col.index.1 = colless(as.treeshape(epitree1))
sak.index.1 = sackin(as.treeshape(epitree1))

# one infection one infection after: Chain of transmission > one to one

id2 = c(seq(from=0,to=22, by=1))
par2 = c(seq(from=-1, to=21, by=1))
t2 = rev(c(seq(from=1,to=23, by=1)))
epi2 <- list()
epi2$itimes <- t2
epi2$dtimes <- rep(0, length(t2))
epi2$id <- id2
epi2$parent <- par2
epitree2 <- epi2tree(epi2)
examp2 <- phylogenetictree.trend(tree = epitree2)
x2 = examp2$num.tree
y2 = examp2$size.tree
reg2 <- lm(log(y2) ~ log(x2))
cozf2 = coef(reg2)

col.index.2 = colless(as.treeshape(epitree2))
sak.index.2 = sackin(as.treeshape(epitree2))

# One infection many infections after more than 3: Super-spreader transmission > one to five or more
id3 = c(seq(from=0,to=22, by=1))
par3 = c(-1,0,0,0,1,1,1,1,1,4,4,1,5,5,6,
         7,7,8,8,9,9,5,0)#,5,5, 10,10,
#3,3,7,14,7,15,16,3,17,17,3,18,5)
t3 = rev(c(seq(from=1,to=23, by=1)))
epi3 <- list()
epi3$itimes <- t3
epi3$dtimes <- rep(0, length(t3))
epi3$id <- id3
epi3$parent <- par3
epitree3 <- epi2tree(epi3)
examp3 <- phylogenetictree.trend(tree = epitree3)
x3 = examp3$num.tree
y3 = examp3$size.tree
reg3 <- lm(log(y3) ~ log(x3))
cozf3 = coef(reg3)


## Types of transmission networks

##########################################
# 1. Chain of transmission > one to one ##
##########################################

# 23 individuals

# one infection > one infection
a.id1 = c(seq(from=0,to=22, by=1))
a.par1 = c(seq(from=-1, to=21, by=1))
a.t1 = rev(c(seq(from=1,to=23, by=1)))
a.epi1 <- list()
a.epi1$itimes <- a.t1
a.epi1$dtimes <- rep(0, length(a.t1))
a.epi1$id <- a.id1
a.epi1$parent <- a.par1
a.epitree1 <- epi2tree(a.epi1)

a.examp1 <- phylogenetictree.trend(tree = a.epitree1)
a.x1 = a.examp1$num.tree
a.y1 = a.examp1$size.tree
a.reg1 <- lm(log(a.y1) ~ log(a.x1))
a.cozf1 = coef(a.reg1)

a.col.index.1 = colless(as.treeshape(a.epitree1))
a.sak.index.1 = sackin(as.treeshape(a.epitree1))

# 51 individuals

# one infection > one infection
a.id2 = c(seq(from=0,to=50, by=1))
a.par2 = c(seq(from=-1, to=49, by=1))
a.t2 = rev(c(seq(from=1,to=51, by=1)))
a.epi2 <- list()
a.epi2$itimes <- a.t2
a.epi2$dtimes <- rep(0, length(a.t2))
a.epi2$id <- a.id2
a.epi2$parent <- a.par2
a.epitree2 <- epi2tree(a.epi2)

a.examp2 <- phylogenetictree.trend(tree = a.epitree2)
a.x2 = a.examp2$num.tree
a.y2 = a.examp2$size.tree
a.reg2 <- lm(log(a.y2) ~ log(a.x2))
a.cozf2 = coef(a.reg2)

a.col.index.2 = colless(as.treeshape(a.epitree2))
a.sak.index.2 = sackin(as.treeshape(a.epitree2))

# 101 individuals

# one infection > one infection
a.id3 = c(seq(from=0,to=100, by=1))
a.par3 = c(seq(from=-1, to=99, by=1))
a.t3 = rev(c(seq(from=1,to=101, by=1)))
a.epi3 <- list()
a.epi3$itimes <- a.t3
a.epi3$dtimes <- rep(0, length(a.t3))
a.epi3$id <- a.id3
a.epi3$parent <- a.par3
a.epitree3 <- epi2tree(a.epi3)

a.examp3 <- phylogenetictree.trend(tree = a.epitree3)
a.x3 = a.examp3$num.tree
a.y3 = a.examp3$size.tree
a.reg3 <- lm(log(a.y3) ~ log(a.x3))
a.cozf3 = coef(a.reg3)

a.col.index.3 = colless(as.treeshape(a.epitree3))
a.sak.index.3 = sackin(as.treeshape(a.epitree3))



# 151 individuals

# one infection > one infection
a.id4 = c(seq(from=0,to=150, by=1))
a.par4 = c(seq(from=-1, to=149, by=1))
a.t4 = rev(c(seq(from=1,to=151, by=1)))
a.epi4 <- list()
a.epi4$itimes <- a.t4
a.epi4$dtimes <- rep(0, length(a.t4))
a.epi4$id <- a.id4
a.epi4$parent <- a.par4
a.epitree4 <- epi2tree(a.epi4)


a.examp4 <- phylogenetictree.trend(tree = a.epitree4)
a.x4 = a.examp4$num.tree
a.y4 = a.examp4$size.tree
a.reg4 <- lm(log(a.y4) ~ log(a.x4))
a.cozf4 = coef(a.reg4)

a.col.index.4 = colless(as.treeshape(a.epitree4))
a.sak.index.4 = sackin(as.treeshape(a.epitree4))


# 201 individuals

# one infection > one infection
a.id5 = c(seq(from=0,to=200, by=1))
a.par5 = c(seq(from=-1, to=199, by=1))
a.t5 = rev(c(seq(from=1,to=201, by=1)))
a.epi5 <- list()
a.epi5$itimes <- a.t5
a.epi5$dtimes <- rep(0, length(a.t5))
a.epi5$id <- a.id5
a.epi5$parent <- a.par5
a.epitree5 <- epi2tree(a.epi5)

a.examp5 <- phylogenetictree.trend(tree = a.epitree5)
a.x5 = a.examp5$num.tree
a.y5 = a.examp5$size.tree
a.reg5 <- lm(log(a.y5) ~ log(a.x5))
a.cozf5 = coef(a.reg5)

a.col.index.5 = colless(as.treeshape(a.epitree5))
a.sak.index.5 = sackin(as.treeshape(a.epitree5))

####################################################
# 2. Homogeneous transmission: one > two or three ##
####################################################

# 23 individuals

# one infection > two or three infections after
b.id1 = c(seq(from=0,to=22, by=1)) # 39
b.par1 = c(-1,0,0,1,1,1,2,3,3,4,4,5,5,5,6,
         7,7,8,8,9,9,10,10) #,11,11,12,12,
# 13,13,14,14,15,15,16,16,17,17,18,18,19)
b.t1 = rev(c(seq(from=1,to=23, by=1)))
b.epi1 <- list()
b.epi1$itimes <- b.t1
b.epi1$dtimes <- rep(0, length(b.t1))
b.epi1$id <- b.id1
b.epi1$parent <- b.par1
b.epitree1 <- epi2tree(b.epi1)


b.examp1 <- phylogenetictree.trend(tree = b.epitree1)
b.x1 = b.examp1$num.tree
b.y1 = b.examp1$size.tree
b.reg1 <- lm(log(b.y1) ~ log(b.x1))
b.cozf1 = coef(b.reg1)

b.col.index.1 = colless(as.treeshape(b.epitree1))
b.sak.index.1 = sackin(as.treeshape(b.epitree1))


# 51 individuals

# one infection > two  or three infections after
b.id2 = c(seq(from=0,to=50, by=1)) # 39
b.par2 = c(-1,0,0,1,1,1,2,3,3,4,4,5,5,5,6,
         7,7,8,8,9,9,10,10,11,11,11,12,12,
         13,13,14,15,15,15,12,2,16,16,17,
         17,18,18,18,19,19,20,20,21,21,22,22) #,11,11,12,12,
# 13,13,14,14,15,15,16,16,17,17,18,18,19)
b.t2 = rev(c(seq(from=1,to=51, by=1)))
b.epi2 <- list()
b.epi2$itimes <- b.t2
b.epi2$dtimes <- rep(0, length(b.t2))
b.epi2$id <- b.id2
b.epi2$parent <- b.par2
b.epitree2 <- epi2tree(b.epi2)

b.examp2 <- phylogenetictree.trend(tree = b.epitree2)
b.x2 = b.examp2$num.tree
b.y2 = b.examp2$size.tree
b.reg2 <- lm(log(b.y2) ~ log(b.x2))
b.cozf2 = coef(b.reg2)

b.col.index.2 = colless(as.treeshape(b.epitree2))
b.sak.index.2 = sackin(as.treeshape(b.epitree2))

# 101 individuals

# one infection > two  or three infections after
b.id3 = c(seq(from=0,to=100, by=1)) # 39
b.par3 = c(-1,0,0,1,1,1,2,3,3,4,4,5,5,5,6,
         7,7,8,8,9,9,10,10,11,11,11,12,12,
         13,13,14,15,15,15,12,2,16,16,17,
         17,18,18,18,19,19,20,20,21,21,22,22,
         21,23,23,24,24,24,25,25,26,26,27,28,
         29,29,30,30,30,31,31,32,32,33,33,34,
         34,31,35,35,36,36,37,14,38,38,37,39,
         40,40,39,38,41,41,42,42,43,43,44,45,
         44,45) #,11,11,12,12,
# 13,13,14,14,15,15,16,16,17,17,18,18,19)
b.t3 = rev(c(seq(from=1,to=101, by=1)))
b.epi3 <- list()
b.epi3$itimes <- b.t3
b.epi3$dtimes <- rep(0, length(b.t3))
b.epi3$id <- b.id3
b.epi3$parent <- b.par3
b.epitree3 <- epi2tree(b.epi3)


b.examp3 <- phylogenetictree.trend(tree = b.epitree3)
b.x3 = b.examp3$num.tree
b.y3 = b.examp3$size.tree
b.reg3 <- lm(log(b.y3) ~ log(b.x3))
b.cozf3 = coef(b.reg3)

b.col.index.3 = colless(as.treeshape(b.epitree3))
b.sak.index.3 = sackin(as.treeshape(b.epitree3))


# 151 individuals

# one infection > two  or three infections after
b.id4 = c(seq(from=0,to=150, by=1)) # 39
b.par4 = c(-1,0,0,1,1,1,2,3,3,4,4,5,5,5,6,
         7,7,8,8,9,9,10,10,11,11,11,12,12,
         13,13,14,15,15,15,12,2,16,16,17,
         17,18,18,18,19,19,20,20,21,21,22,22,
         21,23,23,24,24,24,25,25,26,26,27,28,
         29,29,30,30,30,31,31,32,32,33,33,34,
         34,31,35,35,36,36,37,14,38,38,37,39,
         40,40,39,38,41,41,42,42,43,43,44,45,
         44,45,46,46,47,47,48,49,50,48,48,49,
         51,50,51,52,52,53,53,54,56,57,58,59,
         53,49,54,60,61,62,54,57,54,59,59,57,
         61,62,62,63,63,64,65,66,67,62,63,63,
         64,67,68,69) #,11,11,12,12,
# 13,13,14,14,15,15,16,16,17,17,18,18,19)
b.t4 = rev(c(seq(from=1,to=151, by=1)))
b.epi4 <- list()
b.epi4$itimes <- b.t4
b.epi4$dtimes <- rep(0, length(b.t4))
b.epi4$id <- b.id4
b.epi4$parent <- b.par4
b.epitree4 <- epi2tree(b.epi4)


b.examp4 <- phylogenetictree.trend(tree = b.epitree4)
b.x4 = b.examp4$num.tree
b.y4 = b.examp4$size.tree
b.reg4 <- lm(log(b.y4) ~ log(b.x4))
b.cozf4 = coef(b.reg4)

b.col.index.4 = colless(as.treeshape(b.epitree4))
b.sak.index.4 = sackin(as.treeshape(b.epitree4))


# 201 individuals

# one infection > two  or three infections after
b.id5 = c(seq(from=0,to=200, by=1)) # 39
b.par5 = c(-1,0,0,1,1,1,2,3,3,4,4,5,5,5,6,
         7,7,8,8,9,9,10,10,11,11,11,12,12,
         13,13,14,15,15,15,12,2,16,16,17,
         17,18,18,18,19,19,20,20,21,21,22,22,
         21,23,23,24,24,24,25,25,26,26,27,28,
         29,29,30,30,30,31,31,32,32,33,33,34,
         34,31,35,35,36,36,37,14,38,38,37,39,
         40,40,39,38,41,41,42,42,43,43,44,45,
         44,45,46,46,47,47,48,49,50,48,48,49,
         51,50,51,52,52,53,53,54,56,57,58,59,
         53,49,54,60,61,62,54,57,54,59,59,57,
         61,62,62,63,63,64,65,66,67,62,63,63,
         64,67,68,69,70,71,70,72,73,71,71,74,
         75,76,72,75,74,77,78,76,77,78,76,75,
         79,80,81,78,77,79,81,82,83,84,85,86,
         87,82,84,82,86,87,88,89,90,92,91,92,
         93,92,91,90,92,88) #,11,11,12,12,
# 13,13,14,14,15,15,16,16,17,17,18,18,19)
b.t5 = rev(c(seq(from=1,to=201, by=1)))
b.epi5 <- list()
b.epi5$itimes <- b.t5
b.epi5$dtimes <- rep(0, length(b.t5))
b.epi5$id <- b.id5
b.epi5$parent <- b.par5
b.epitree5 <- epi2tree(b.epi5)


b.examp5 <- phylogenetictree.trend(tree = b.epitree5)
b.x5 = b.examp5$num.tree
b.y5 = b.examp5$size.tree
b.reg5 <- lm(log(b.y5) ~ log(b.x5))
b.cozf5 = coef(b.reg5)

b.col.index.5 = colless(as.treeshape(b.epitree5))
b.sak.index.5 = sackin(as.treeshape(b.epitree5))


#########################################################
# 3. Super-spreader transmission > one to five or more ##
#########################################################

# 23 individuals

# one infection > two to five infections after
c.id1 = c(seq(from=0,to=22, by=1)) # 39
c.par1 = c(-1,0,0,1,1,1,2,1,3,4,1,4,5,5,6,
         7,4,8,8,4,9,9,4) #,11,11,12,12,
# 13,13,14,14,15,15,16,16,17,17,18,18,19)
c.t1 = rev(c(seq(from=1,to=23, by=1)))
c.epi1 <- list()
c.epi1$itimes <- c.t1
c.epi1$dtimes <- rep(0, length(c.t1))
c.epi1$id <- c.id1
c.epi1$parent <- c.par1
c.epitree1 <- epi2tree(c.epi1)

c.examp1 <- phylogenetictree.trend(tree = c.epitree1)
c.x1 = c.examp1$num.tree
c.y1 = c.examp1$size.tree
c.reg1 <- lm(log(c.y1) ~ log(c.x1))
c.cozf1 = coef(c.reg1)

c.col.index.1 = colless(as.treeshape(c.epitree1))
c.sak.index.1 = sackin(as.treeshape(c.epitree1))


# 51 individuals

# one infection > two  to five  infections after
c.id2 = c(seq(from=0,to=50, by=1)) # 39
c.par2 = c(-1,0,1,1,1,2,1,3,4,1,4,5,5,6,
         7,4,8,8,4,9,9,4,10,11,11,11,12,12,
         13,13,14,15,15,15,12,12,16,16,12,
         12,18,18,18,15,18,20,20,18,11,15,8) #,11,11,12,12,
# 13,13,14,14,15,15,16,16,17,17,18,18,19)
c.t2 = rev(c(seq(from=1,to=51, by=1)))
c.epi2 <- list()
c.epi2$itimes <- c.t2
c.epi2$dtimes <- rep(0, length(c.t2))
c.epi2$id <- c.id2
c.epi2$parent <- c.par2
c.epitree2 <- epi2tree(c.epi2)


c.examp2 <- phylogenetictree.trend(tree = c.epitree2)
c.x2 = c.examp2$num.tree
c.y2 = c.examp2$size.tree
c.reg2 <- lm(log(c.y2) ~ log(c.x2))
c.cozf2 = coef(c.reg2)

c.col.index.2 = colless(as.treeshape(c.epitree2))
c.sak.index.2 = sackin(as.treeshape(c.epitree2))



# 101 individuals

# one infection > two  to five  infections after
c.id3 = c(seq(from=0,to=100, by=1)) # 39
c.par3 = c(-1,0,1,1,1,2,1,3,4,1,4,5,5,6,
         7,4,8,8,4,9,9,4,10,11,11,11,12,12,
         13,13,14,15,14,15,12,12,16,16,12,
         12,18,18,18,15,18,20,20,18,11,15,8,
         21,23,23,24,24,24,25,25,26,26,27,28,
         21,38,39,21,21,32,32,33,33,34,12,34,
         34,21,35,32,35,33,32,15,38,38,14,39,
         43,43,39,38,41,41,42,42,43,41,41,43,
         43,41) #,11,11,12,12,
# 13,13,14,14,15,15,16,16,17,17,18,18,19)
c.t3 = rev(c(seq(from=1,to=101, by=1)))
c.epi3 <- list()
c.epi3$itimes <- c.t3
c.epi3$dtimes <- rep(0, length(c.t3))
c.epi3$id <- c.id3
c.epi3$parent <- c.par3
c.epitree3 <- epi2tree(c.epi3)


c.examp3 <- phylogenetictree.trend(tree = c.epitree3)
c.x3 = c.examp3$num.tree
c.y3 = c.examp3$size.tree
c.reg3 <- lm(log(c.y3) ~ log(c.x3))
c.cozf3 = coef(c.reg3)

c.col.index.3 = colless(as.treeshape(c.epitree3))
c.sak.index.3 = sackin(as.treeshape(c.epitree3))


# 151 individuals

# one infection > two  to five  infections after
c.id4 = c(seq(from=0,to=150, by=1)) # 39
c.par4 = c(-1,0,1,1,1,2,1,3,4,1,4,5,5,6,
         7,4,8,8,4,9,9,4,10,11,11,11,12,12,
         13,13,14,15,14,15,12,12,16,16,12,
         12,18,18,18,15,18,20,20,18,11,15,8,
         21,23,23,24,24,24,25,25,26,26,27,28,
         21,38,39,21,21,32,32,33,33,34,12,34,
         34,21,35,32,35,33,32,15,38,38,14,39,
         43,43,39,38,41,41,42,42,43,41,41,43,
         43,41,66,46,17,47,48,49,50,48,48,49,
         51,50,51,52,52,53,53,54,66,17,58,59,
         47,49,54,20,61,62,54,47,54,17,20,47,
         66,17,62,63,63,64,65,66,17,62,63,66,
         64,47,62,63) #,11,11,12,12,
# 13,13,14,14,15,15,16,16,17,17,18,18,19)
c.t4 = rev(c(seq(from=1,to=151, by=1)))
c.epi4 <- list()
c.epi4$itimes <- c.t4
c.epi4$dtimes <- rep(0, length(c.t4))
c.epi4$id <- c.id4
c.epi4$parent <- c.par4
c.epitree4 <- epi2tree(c.epi4)


c.examp4 <- phylogenetictree.trend(tree = c.epitree4)
c.x4 = c.examp4$num.tree
c.y4 = c.examp4$size.tree
c.reg4 <- lm(log(c.y4) ~ log(c.x4))
c.cozf4 = coef(c.reg4)

c.col.index.4 = colless(as.treeshape(c.epitree4))
c.sak.index.4 = sackin(as.treeshape(c.epitree4))


# 201 individuals

# one infection > two  to five  infections after
c.id5 = c(seq(from=0,to=200, by=1)) # 39
c.par5 = c(-1,0,1,1,1,2,1,3,4,1,4,5,5,6,
         7,4,8,8,4,9,9,4,10,11,11,11,12,12,
         13,13,14,15,14,15,12,12,16,16,12,
         12,18,18,18,15,18,20,20,18,11,15,8,
         21,23,23,24,24,24,25,25,26,26,27,28,
         21,38,39,21,21,32,32,33,33,34,12,34,
         34,21,35,32,35,33,32,15,38,38,14,39,
         43,43,39,38,41,41,42,42,43,41,41,43,
         43,41,66,46,17,47,48,49,50,48,48,49,
         51,50,51,52,52,53,53,54,66,17,58,59,
         47,49,54,20,61,62,54,47,54,17,20,47,
         66,17,62,63,63,64,65,66,17,62,63,66,
         64,47,62,63,70,71,70,72,73,71,71,74,
         75,76,72,75,74,77,78,76,77,78,76,75,
         81,80,81,78,70,79,81,82,81,84,81,86,
         87,82,84,82,86,87,88,89,70,92,81,92,
         92,92,81,70,92,88) #,11,11,12,12,
# 13,13,14,14,15,15,16,16,17,17,18,18,19)
c.t5 = rev(c(seq(from=1,to=201, by=1)))
c.epi5 <- list()
c.epi5$itimes <- c.t5
c.epi5$dtimes <- rep(0, length(c.t5))
c.epi5$id <- c.id5
c.epi5$parent <- c.par5
c.epitree5 <- epi2tree(c.epi5)


c.examp5 <- phylogenetictree.trend(tree = c.epitree5)
c.x5 = c.examp5$num.tree
c.y5 = c.examp5$size.tree
c.reg5 <- lm(log(c.y5) ~ log(c.x5))
c.cozf5 = coef(c.reg5)

c.col.index.5 = colless(as.treeshape(c.epitree5))
c.sak.index.5 = sackin(as.treeshape(c.epitree5))



# 1. Chain of transmission > one to one
inter.a <- c(a.cozf1[[1]]/a.sak.index.1,a.cozf2[[1]]/a.sak.index.2,a.cozf3[[1]]/a.sak.index.3,
             a.cozf4[[1]]/a.sak.index.4,a.cozf5[[1]]/a.sak.index.5)

inter.a.1 <- c(a.cozf1[[1]]/a.col.index.1,a.cozf2[[1]]/a.col.index.2,a.cozf3[[1]]/a.col.index.3,
             a.cozf4[[1]]/a.col.index.4,a.cozf5[[1]]/a.col.index.5)

slope.a <- c(a.cozf1[[2]],a.cozf2[[2]],a.cozf3[[2]],a.cozf4[[2]],a.cozf5[[2]])

# 2. Homogeneous transmission: one > two or three
inter.b <- c(b.cozf1[[1]]/b.sak.index.1,b.cozf2[[1]]/b.sak.index.2,b.cozf3[[1]]/b.sak.index.3,
             b.cozf4[[1]]/b.sak.index.4,b.cozf5[[1]]/b.sak.index.5)

inter.b.1 <- c(b.cozf1[[1]]/b.col.index.1,b.cozf2[[1]]/b.col.index.2,b.cozf3[[1]]/b.col.index.3,
             b.cozf4[[1]]/b.col.index.4,b.cozf5[[1]]/b.col.index.5)

slope.b <- -c(b.cozf1[[2]],b.cozf2[[2]],b.cozf3[[2]],b.cozf4[[2]],b.cozf5[[2]])

# 3. Super-spreader transmission > one to five or more
inter.c <- c(c.cozf1[[1]]/c.sak.index.1,c.cozf2[[1]]/c.sak.index.2,c.cozf3[[1]]/c.sak.index.3,
             c.cozf4[[1]]/c.sak.index.4,c.cozf5[[1]]/c.sak.index.5)

inter.c.1 <- c(c.cozf1[[1]]/c.col.index.1,c.cozf2[[1]]/c.col.index.2,c.cozf3[[1]]/c.col.index.3,
             c.cozf4[[1]]/c.col.index.4,c.cozf5[[1]]/c.col.index.5)

slope.c <- -c(c.cozf1[[2]],c.cozf2[[2]],c.cozf3[[2]],c.cozf4[[2]],c.cozf5[[2]])


# library(adephylo)
#
# c = distTips(c.epitree5, tips = "all", method = "Abouheif")
#
# f = as.matrix(c)


# # Interecepts and ratio between colless index and sakin' index
#
# inter.a <- c(a.cozf1[[1]],a.cozf2[[1]],a.cozf3[[1]],a.cozf4[[1]],a.cozf5[[1]], a.cozf6[[1]])
# inter.b <- c(b.cozf1[[1]],b.cozf2[[1]],b.cozf3[[1]],b.cozf4[[1]],b.cozf5[[1]])
# inter.c <- c(c.cozf1[[1]],c.cozf2[[1]],c.cozf3[[1]],c.cozf4[[1]],c.cozf5[[1]])
#
# ratio.a <- c(a.col.index.1/a.sak.index.1,a.col.index.2/a.sak.index.2,a.col.index.3/a.sak.index.3,
#              a.col.index.4/a.sak.index.4,a.col.index.5/a.sak.index.5, a.col.index.6/a.sak.index.6)
#
# ratio.b <- c(b.col.index.1/b.sak.index.1,b.col.index.2/b.sak.index.2,b.col.index.3/b.sak.index.3,
#              b.col.index.4/b.sak.index.4,b.col.index.5/b.sak.index.5)
#
# ratio.c <- c(c.col.index.1/c.sak.index.1,c.col.index.2/c.sak.index.2,c.col.index.3/c.sak.index.3,
#              c.col.index.4/c.sak.index.4,c.col.index.5/c.sak.index.5)
#
# xlim <- range(inter.a,inter.b,inter.c)
# ylim <- range(ratio.a,ratio.b,ratio.c)
#
# plot(inter.a,ratio.a, type="b", col="red", xlim=c(0, xlim[2]), ylim=c(0,ylim[2]), lwd=2)
# lines(inter.b,ratio.b, col='green', type='b', lwd=2)
# lines(inter.c, ratio.c, col='blue', type='b', lwd=2)


# Dig in deep

# Case 1: 1 infection > two

# 23 individuals

d.id1 = c(seq(from=0,to=22, by=1))
d.par1 = c(-1,0,0,1,1,2,2,3,3,4,4,
           5,5,6,6,7,7,8,8,9,9,10,10)

d.t1 = rev(c(seq(from=1,to=23, by=1)))
d.epi1 <- list()
d.epi1$itimes <- d.t1
d.epi1$dtimes <- rep(0, length(d.t1))
d.epi1$id <- d.id1
d.epi1$parent <- d.par1
d.epitree1 <- epi2tree(d.epi1)

d.examp1 <- phylogenetictree.trend(tree = d.epitree1)
d.x1 = d.examp1$num.tree
d.y1 = d.examp1$size.tree
d.reg1 <- lm(log(d.y1) ~ log(d.x1))
d.cozf1 = coef(d.reg1)

d.col.index.1 = colless(as.treeshape(d.epitree1))
d.sak.index.1 = sackin(as.treeshape(d.epitree1))


# 51 individuals

d.id2 = c(seq(from=0,to=50, by=1))
d.par2 = c(-1,0,0,1,1,2,2,3,3,4,4,
           5,5,6,6,7,7,8,8,9,9,10,10,
           11,11,12,12,13,13,14,14,15,15,
           16,16,17,17,18,18,19,19,20,20,
           21,21,22,22,23,23,24,24)

d.t2 = rev(c(seq(from=1,to=51, by=1)))
d.epi2 <- list()
d.epi2$itimes <- d.t2
d.epi2$dtimes <- rep(0, length(d.t2))
d.epi2$id <- d.id2
d.epi2$parent <- d.par2
d.epitree2 <- epi2tree(d.epi2)

d.examp2 <- phylogenetictree.trend(tree = d.epitree2)
d.x2 = d.examp2$num.tree
d.y2 = d.examp2$size.tree
d.reg2 <- lm(log(d.y2) ~ log(d.x2))
d.cozf2 = coef(d.reg2)

d.col.index.2 = colless(as.treeshape(d.epitree2))
d.sak.index.2 = sackin(as.treeshape(d.epitree2))

# 101 individuals

d.id3 = c(seq(from=0,to=100, by=1))
d.par3 = c(-1,0,0,1,1,2,2,3,3,4,4,
           5,5,6,6,7,7,8,8,9,9,10,10,
           11,11,12,12,13,13,14,14,15,15,
           16,16,17,17,18,18,19,19,20,20,
           21,21,22,22,23,23,24,24,25,25,
           26,26,27,27,28,28,29,29,30,30,
           31,31,32,32,33,33,34,34,35,35,
           36,36,37,37,38,38,39,39,40,40,
           41,41,42,42,43,43,44,44,45,45,
           46,46,47,47,48,48,49,49)
d.t3 = rev(c(seq(from=1,to=101, by=1)))
d.epi3 <- list()
d.epi3$itimes <- d.t3
d.epi3$dtimes <- rep(0, length(d.t3))
d.epi3$id <- d.id3
d.epi3$parent <- d.par3
d.epitree3 <- epi2tree(d.epi3)

d.examp3 <- phylogenetictree.trend(tree = d.epitree3)
d.x3 = d.examp3$num.tree
d.y3 = d.examp3$size.tree
d.reg3 <- lm(log(d.y3) ~ log(d.x3))
d.cozf3 = coef(d.reg3)

d.col.index.3 = colless(as.treeshape(d.epitree3))
d.sak.index.3 = sackin(as.treeshape(d.epitree3))


# 151 individuals

d.id4 = c(seq(from=0,to=150, by=1))
d.par4 = c(-1,0,0,1,1,2,2,3,3,4,4,
           5,5,6,6,7,7,8,8,9,9,10,10,
           11,11,12,12,13,13,14,14,15,15,
           16,16,17,17,18,18,19,19,20,20,
           21,21,22,22,23,23,24,24,25,25,
           26,26,27,27,28,28,29,29,30,30,
           31,31,32,32,33,33,34,34,35,35,
           36,36,37,37,38,38,39,39,40,40,
           41,41,42,42,43,43,44,44,45,45,
           46,46,47,47,48,48,49,49,50,50,
           51,51,52,52,53,53,54,54,55,55,
           56,56,57,57,58,58,59,59,60,60,
           61,61,62,62,63,63,64,64,65,65,
           66,66,67,67,68,68,69,69,70,70,
           71,71,72,72,73,73,74,74)
d.t4 = rev(c(seq(from=1,to=151, by=1)))
d.epi4 <- list()
d.epi4$itimes <- d.t4
d.epi4$dtimes <- rep(0, length(d.t4))
d.epi4$id <- d.id4
d.epi4$parent <- d.par4
d.epitree4 <- epi2tree(d.epi4)

d.examp4 <- phylogenetictree.trend(tree = d.epitree4)
d.x4 = d.examp4$num.tree
d.y4 = d.examp4$size.tree
d.reg4 <- lm(log(d.y4) ~ log(d.x4))
d.cozf4 = coef(d.reg4)

d.col.index.4 = colless(as.treeshape(d.epitree4))
d.sak.index.4 = sackin(as.treeshape(d.epitree4))


# 201 individuals

d.id5 = c(seq(from=0,to=200, by=1))
d.par5 = c(-1,0,0,1,1,2,2,3,3,4,4,
           5,5,6,6,7,7,8,8,9,9,10,10,
           11,11,12,12,13,13,14,14,15,15,
           16,16,17,17,18,18,19,19,20,20,
           21,21,22,22,23,23,24,24,25,25,
           26,26,27,27,28,28,29,29,30,30,
           31,31,32,32,33,33,34,34,35,35,
           36,36,37,37,38,38,39,39,40,40,
           41,41,42,42,43,43,44,44,45,45,
           46,46,47,47,48,48,49,49,50,50,
           51,51,52,52,53,53,54,54,55,55,
           56,56,57,57,58,58,59,59,60,60,
           61,61,62,62,63,63,64,64,65,65,
           66,66,67,67,68,68,69,69,70,70,
           71,71,72,72,73,73,74,74,75,75,
           76,76,77,77,78,78,79,79,80,80,
           81,81,82,82,83,83,84,84,85,85,
           86,86,87,87,89,89,90,90,91,91,
           92,92,93,93,94,94,95,95,96,96,
           97,97,98,98,99,99,100,100)
d.t5 = rev(c(seq(from=1,to=201, by=1)))
d.epi5 <- list()

d.epi5$itimes <- d.t5
d.epi5$dtimes <- rep(0, length(d.t5))
d.epi5$id <- d.id5
d.epi5$parent <- d.par5
d.epitree5 <- epi2tree(d.epi5)

d.examp5 <- phylogenetictree.trend(tree = d.epitree5)
d.x5 = d.examp5$num.tree
d.y5 = d.examp5$size.tree
d.reg5 <- lm(log(d.y5) ~ log(d.x5))
d.cozf5 = coef(d.reg5)

d.col.index.5 = colless(as.treeshape(d.epitree5))
d.sak.index.5 = sackin(as.treeshape(d.epitree5))


# Case 2: very super-spreader , 1 > 10

# 23 individuals

e.id1 = c(seq(from=0,to=22, by=1))
e.par1 = c(-1,0,0,1,1,1,1,1,1,1,1,1,1,
           2,2,2,2,2,2,2,2,2,2)

e.t1 = rev(c(seq(from=1,to=23, by=1)))
e.epi1 <- list()
e.epi1$itimes <- e.t1
e.epi1$dtimes <- rep(0, length(e.t1))
e.epi1$id <- e.id1
e.epi1$parent <- e.par1
e.epitree1 <- epi2tree(e.epi1)

e.examp1 <- phylogenetictree.trend(tree = e.epitree1)
e.x1 = d.examp1$num.tree
e.y1 = d.examp1$size.tree
e.reg1 <- lm(log(e.y1) ~ log(e.x1))
e.cozf1 = coef(d.reg1)

e.col.index.1 = colless(as.treeshape(e.epitree1))
e.sak.index.1 = sackin(as.treeshape(e.epitree1))



# 51 individuals

e.id2 = c(seq(from=0,to=50, by=1))
e.par2 = c(-1,0,0,1,1,1,1,1,1,1,1,1,1,
           2,2,2,2,2,2,2,2,2,2,3,3,3,3,
           3,3,3,3,3,3,4,4,4,4,4,4,4,4,
           4,4,5,5,5,5,5,5,5,5)

e.t2 = rev(c(seq(from=1,to=51, by=1)))
e.epi2 <- list()
e.epi2$itimes <- e.t2
e.epi2$dtimes <- rep(0, length(e.t2))
e.epi2$id <- e.id2
e.epi2$parent <- e.par2
e.epitree2 <- epi2tree(e.epi2)

e.examp2 <- phylogenetictree.trend(tree = e.epitree2)
e.x2 = d.examp2$num.tree
e.y2 = d.examp2$size.tree
e.reg2 <- lm(log(e.y2) ~ log(e.x2))
e.cozf2 = coef(d.reg2)

e.col.index.2 = colless(as.treeshape(e.epitree2))
e.sak.index.2 = sackin(as.treeshape(e.epitree2))


# 101 individuals

e.id3 = c(seq(from=0,to=100, by=1))
e.par3 = c(-1,0,0,1,1,1,1,1,1,1,1,1,1,
           2,2,2,2,2,2,2,2,2,2,3,3,3,3,
           3,3,3,3,3,3,4,4,4,4,4,4,4,4,
           4,4,5,5,5,5,5,5,5,5,6,6,6,6,
           6,6,6,6,6,6,7,7,7,7,7,7,7,7,
           7,7,8,8,8,8,8,8,8,8,8,8,9,9,
           9,9,9,9,9,9,9,9,10,10,10,10,
           10,10,10,10,10,10)

e.t3 = rev(c(seq(from=1,to=101, by=1)))

e.epi3 <- list()
e.epi3$itimes <- e.t3
e.epi3$dtimes <- rep(0, length(e.t3))
e.epi3$id <- e.id3
e.epi3$parent <- e.par3
e.epitree3 <- epi2tree(e.epi3)

e.examp3 <- phylogenetictree.trend(tree = e.epitree3)
e.x3 = d.examp3$num.tree
e.y3 = d.examp3$size.tree
e.reg3 <- lm(log(e.y3) ~ log(e.x3))
e.cozf3 = coef(d.reg3)

e.col.index.3 = colless(as.treeshape(e.epitree3))
e.sak.index.3 = sackin(as.treeshape(e.epitree3))

# 151 individuals

e.id4 = c(seq(from=0,to=150, by=1))
e.par4 = c(-1,0,0,1,1,1,1,1,1,1,1,1,1,
           2,2,2,2,2,2,2,2,2,2,3,3,3,3,
           3,3,3,3,3,3,4,4,4,4,4,4,4,4,
           4,4,5,5,5,5,5,5,5,5,6,6,6,6,
           6,6,6,6,6,6,7,7,7,7,7,7,7,7,
           7,7,8,8,8,8,8,8,8,8,8,8,9,9,
           9,9,9,9,9,9,9,9,10,10,10,10,
           10,10,10,10,10,10,11,11,11,11,
           11,11,11,11,11,11,12,12,12,12,
           12,12,12,12,12,12,13,13,13,13,
           13,13,13,13,13,13,14,14,14,14,
           14,14,14,14,14,14,15,15,15,15,
           15,15,15,15,15,15)

e.t4 = rev(c(seq(from=1,to=151, by=1)))

e.epi4 <- list()
e.epi4$itimes <- e.t4
e.epi4$dtimes <- rep(0, length(e.t4))
e.epi4$id <- e.id4
e.epi4$parent <- e.par4
e.epitree4 <- epi2tree(e.epi4)

e.examp4 <- phylogenetictree.trend(tree = e.epitree4)
e.x4 = d.examp4$num.tree
e.y4 = d.examp4$size.tree
e.reg4 <- lm(log(e.y4) ~ log(e.x4))
e.cozf4 = coef(d.reg4)

e.col.index.4 = colless(as.treeshape(e.epitree4))
e.sak.index.4 = sackin(as.treeshape(e.epitree4))


# 201 individuals

e.id5 = c(seq(from=0,to=200, by=1))
e.par5 = c(-1,0,0,1,1,1,1,1,1,1,1,1,1,
           2,2,2,2,2,2,2,2,2,2,3,3,3,3,
           3,3,3,3,3,3,4,4,4,4,4,4,4,4,
           4,4,5,5,5,5,5,5,5,5,6,6,6,6,
           6,6,6,6,6,6,7,7,7,7,7,7,7,7,
           7,7,8,8,8,8,8,8,8,8,8,8,9,9,
           9,9,9,9,9,9,9,9,10,10,10,10,
           10,10,10,10,10,10,11,11,11,11,
           11,11,11,11,11,11,12,12,12,12,
           12,12,12,12,12,12,13,13,13,13,
           13,13,13,13,13,13,14,14,14,14,
           14,14,14,14,14,14,15,15,15,15,
           15,15,15,15,15,15,16,16,16,16,
           16,16,16,16,16,16,17,17,17,17,
           17,17,17,17,17,17,18,18,18,18,
           18,18,18,18,18,18,19,19,19,19,
           19,19,19,19,19,19,20,20,20,20,
           20,20,20,20,20,20)

e.t5 = rev(c(seq(from=1,to=201, by=1)))

e.epi5 <- list()
e.epi5$itimes <- e.t5
e.epi5$dtimes <- rep(0, length(e.t5))
e.epi5$id <- e.id5
e.epi5$parent <- e.par5
e.epitree5 <- epi2tree(e.epi5)

e.examp5 <- phylogenetictree.trend(tree = e.epitree5)
e.x5 = d.examp5$num.tree
e.y5 = d.examp5$size.tree
e.reg5 <- lm(log(e.y5) ~ log(e.x5))
e.cozf5 = coef(d.reg5)

e.col.index.5 = colless(as.treeshape(e.epitree5))
e.sak.index.5 = sackin(as.treeshape(e.epitree5))



# Interecepts and ratio between colless index and sakin' index

inter.a <- c(a.cozf1[[1]],a.cozf2[[1]],a.cozf3[[1]],a.cozf4[[1]],a.cozf5[[1]])
inter.b <- c(b.cozf1[[1]],b.cozf2[[1]],b.cozf3[[1]],b.cozf4[[1]],b.cozf5[[1]])
inter.c <- c(c.cozf1[[1]],c.cozf2[[1]],c.cozf3[[1]],c.cozf4[[1]],c.cozf5[[1]])
inter.d <- c(d.cozf1[[1]],d.cozf2[[1]],d.cozf3[[1]],d.cozf4[[1]],d.cozf5[[1]])
inter.e <- c(e.cozf1[[1]],e.cozf2[[1]],e.cozf3[[1]],e.cozf4[[1]],e.cozf5[[1]])

ratio.a <- c(a.col.index.1/a.sak.index.1,a.col.index.2/a.sak.index.2,a.col.index.3/a.sak.index.3,
             a.col.index.4/a.sak.index.4,a.col.index.5/a.sak.index.5)

ratio.b <- c(b.col.index.1/b.sak.index.1,b.col.index.2/b.sak.index.2,b.col.index.3/b.sak.index.3,
             b.col.index.4/b.sak.index.4,b.col.index.5/b.sak.index.5)

ratio.c <- c(c.col.index.1/c.sak.index.1,c.col.index.2/c.sak.index.2,c.col.index.3/c.sak.index.3,
             c.col.index.4/c.sak.index.4,c.col.index.5/c.sak.index.5)

ratio.d <- c(d.col.index.1/d.sak.index.1,d.col.index.2/d.sak.index.2,d.col.index.3/d.sak.index.3,
             d.col.index.4/d.sak.index.4,d.col.index.5/d.sak.index.5)

ratio.e <- c(e.col.index.1/e.sak.index.1,e.col.index.2/e.sak.index.2,e.col.index.3/e.sak.index.3,
             e.col.index.4/e.sak.index.4,e.col.index.5/e.sak.index.5)

xlim <- range(inter.a,inter.b,inter.c, inter.d, inter.e)
ylim <- range(ratio.a,ratio.b,ratio.c, ratio.d, ratio.e)

plot(inter.a,ratio.a, type="b", col="red", xlim=c(0, xlim[2]), ylim=c(0,ylim[2]), lwd=2) # 1 > 1
lines(inter.b,ratio.b, col='green', type='b', lwd=2) # 1 > 2,3
lines(inter.c, ratio.c, col='blue', type='b', lwd=2) # 1 > 2,4,5
lines(inter.d, ratio.d, col='magenta2', type='b', lwd=2) # 1 > 2
lines(inter.e, ratio.e, col='yellow', type='b', lwd=2) # 1 > 10



# Test ladder like but add two more
a.id4.1 = c(seq(from=0,to=150, by=1))
a.par4.1 = c(seq(from=-1, to=148, by=1),1)
a.t4.1 = rev(c(seq(from=1,to=151, by=1)))
a.epi4.1 <- list()
a.epi4.1$itimes <- a.t4.1
a.epi4.1$dtimes <- rep(0, length(a.t4.1))
a.epi4.1$id <- a.id4.1
a.epi4.1$parent <- a.par4.1
a.epitree4.1 <- epi2tree(a.epi4.1)


a.examp4.1 <- phylogenetictree.trend(tree = a.epitree4.1)
a.x4.1 = a.examp4.1$num.tree
a.y4.1 = a.examp4.1$size.tree
a.reg4.1 <- lm(log(a.y4.1) ~ log(a.x4.1))
a.cozf4.1 = coef(a.reg4.1)

a.col.index.4.1 = colless(as.treeshape(a.epitree4.1))
a.sak.index.4.1 = sackin(as.treeshape(a.epitree4.1))

# Notes
ratio.a.4.1 <- a.col.index.4.1/a.sak.index.4.1 # similar to a's

# a.cozf4.1[[2]] exist, no NA as it is the case for a's
a.cozf4.1[[2]]


x <- c(23,51,101,151,201)
plot(x,ratio.a, type="b", col="red", lwd=2) # 1 > 1
lines(x, ratio.b, col='green', type='b', lwd=2)
lines(x, ratio.c, col='blue', type='b', lwd=2)
lines(x, ratio.d, col='magenta2', type='b', lwd=2)
lines(x,  ratio.e, col='yellow', type='b', lwd=2)

matplot(x, cbind(ratio.a, ratio.b,
                 ratio.c, ratio.d, ratio.e),type="b",
        col=c("red","green", "blue", "magenta2", "orange"),lwd=2, lty=c(1,1))
legend("topleft", legend = c("1 to 1",
                             "1 to 2 or 3",
                             "1 to 3 or 5",
                             "1 to 2",
                             "1 to 10"),
       col=c("red","green", "blue", "magenta2", "orange"), pch=1)



matplot(x, cbind(inter.a, inter.b,
                 inter.c, inter.d, inter.e),type="b",
        col=c("red","green", "blue", "orange", "yellow"),lwd=2, lty=c(1,1))

legend("topleft", legend = c("1 to 1",
                             "1 to 2 or 3",
                             "1 to 3 or 5",
                             "1 to 2",
                             "1 to 10"),
       col=c("red","green", "blue", "magenta2", "orange"), pch=1)



col.a.n <- c(a.col.index.1/23,a.col.index.2/51,a.col.index.3/101,
             a.col.index.4/151,a.col.index.5/201)

col.b.n <- c(b.col.index.1/23,b.col.index.2/51,b.col.index.3/101,
             b.col.index.4/151,b.col.index.5/201)

col.c.n <- c(c.col.index.1/23,c.col.index.2/51,c.col.index.3/101,
             c.col.index.4/151,c.col.index.5/201)

col.d.n <- c(d.col.index.1/23,d.col.index.2/51,d.col.index.3/101,
             d.col.index.4/101,d.col.index.5/201)

col.e.n <- c(e.col.index.1/23,e.col.index.2/51,e.col.index.3/101,
             e.col.index.4/151,e.col.index.5/201)


matplot(x, cbind(col.a.n, col.b.n,
                 col.c.n, col.d.n, col.e.n),type="b",
        col=c("red","green", "blue", "orange", "yellow"),lwd=2, lty=c(1,1))

legend("topleft", legend = c("1 to 1",
                             "1 to 2 or 3",
                             "1 to 3 or 5",
                             "1 to 2",
                             "1 to 10"),
       col=c("red","green", "blue", "magenta2", "orange"), pch=1)


sak.a.n <- c(a.sak.index.1/23,a.sak.index.2/51,a.sak.index.3/101,
             a.sak.index.4/151,a.sak.index.5/201)

sak.b.n <- c(b.sak.index.1/23,b.sak.index.2/51,b.sak.index.3/101,
             b.sak.index.4/151,b.sak.index.5/201)

sak.c.n <- c(c.sak.index.1/23,c.sak.index.2/51,c.sak.index.3/101,
             c.sak.index.4/151,c.sak.index.5/201)

sak.d.n <- c(d.sak.index.1/23,d.sak.index.2/51,d.sak.index.3/101,
             d.sak.index.4/151,d.sak.index.5/201)

sak.e.n <- c(e.sak.index.1/23,e.sak.index.2/51,e.sak.index.3/101,
             e.sak.index.4/151,e.sak.index.5/201)



matplot(x, cbind(sak.a.n, sak.b.n,
                 sak.c.n, sak.d.n, sak.e.n),type="b",
        col=c("red","green", "blue", "orange", "yellow"),lwd=2, lty=c(1,1))

legend("topleft", legend = c("1 to 1",
                             "1 to 2 or 3",
                             "1 to 3 or 5",
                             "1 to 2",
                             "1 to 10"),
       col=c("red","green", "blue", "magenta2", "orange"), pch=1)
