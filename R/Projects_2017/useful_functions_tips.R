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
