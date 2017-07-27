setwd("/home/david/Dropbox/Fitting_Simpact/FitEpiPhylo/")


# Compute target statistics

# growth rate: OK
# relationship rate: OK
# transmission rate: OK

# Slope of subtrees ~ size of subtrees

##################################
# 1. Growth rate: done by Trust ##
##################################

growth.rate <- pop.growth.calculator(datalist = chunk.datalist.test,
                                     timewindow = c(0, timewindow.max = end.time.wind))

###########################
# 2. Transmission rate ####
###########################

# Transmission rate in time interval
transmission.rate <- function(datalist = datalist,
                                  timewindow = c(0, 40)){

  Infec.pop.table <- datalist$ptable[InfectType==1]

  numb.infec.pop <- nrow(Infec.pop.table %>%
    subset(InfectTime <=timewindow[2] & InfectTime >= timewindow[1]))


  transm.rate <- numb.infec.pop / diff(timewindow)

  return(transm.rate)
}

# Transmission rate in time interval subdivided in different steps

# Transmission rate
transmission.rate.int <- function(datalist = datalist,
                              timewindow = c(0, 40), by=1){

  Infec.pop.table <- datalist$ptable[InfectType==1]

  upper.limit <- ceiling(diff(timewindow)/by)

  interv.time <- round(seq.int(timewindow[1], timewindow[2], length.out = upper.limit), digits = 0)

  infec.pop.int <- vector()
  trans.rate.int <- vector()


  for(i in 0:(upper.limit-2)){

    timewindow.int <- c(interv.time[1+i], interv.time[2+i])

    infec.pop.numb <- nrow(Infec.pop.table %>%
      subset(InfectTime <=timewindow.int[2] & InfectTime >= timewindow.int[1]))


    infec.pop.int <- c(infec.pop.int,  infec.pop.numb)


    trans.rate.int <- c(trans.rate.int, (infec.pop.numb / diff(timewindow.int)))

  }

  return(trans.rate.int)
}


# Transmission ratio for men and women

transm.gender.ratio <- function(datalist = datalist,
                              timewindow = c(0, 40)){

  pers.table.hiv <- datalist$ptable[InfectType==1]

  pers.table.hiv.window <- pers.table.hiv %>%
    subset(InfectTime <=timewindow[2] & InfectTime >= timewindow[1])

  numb.transm <- nrow(pers.table.hiv.window)

  numb.transm.men <- length(unique((pers.table.hiv[pers.table.hiv$Gender=="0"]$ID))) # Gender 0 men

  numb.transm.women <- length(unique((pers.table.hiv[pers.table.hiv$Gender=="0"]$ID)))# Gender 1 women

  transm.ratio <- vector("list", length(c(0,1)))

  for(i in 1:2){
    transm.ratio$men <- (numb.transm.men/numb.transm)  #/diff(timewindow)
    transm.ratio$women <- (numb.transm.women/numb.transm)  #/diff(timewindow)
  }
  return(transm.ratio)
}

###############################
# 3. Relationship rate ########
###############################

# Women, men and overall relationships rate in time interval
# Overall relationship rate
# Rste of women (men) who have been invlved in relationship udirng a time interval.s

rels.rate <- function(datalist = datalist,
                              timewindow = c(0, 40)){

  Rels.table <- datalist$rtable

  Rels.table.window <- Rels.table %>%
    subset(FormTime <=timewindow[2] & FormTime >= timewindow[1])

  numb.rels <- nrow(Rels.table.window)

  numb.rels.men <- length(unique(Rels.table.window$ID1)) # Gender 0 men & 1 women

  numb.rels.women <- length(unique(Rels.table.window$ID2))

  rels.rate <- vector("list", length(c(0,1,2)))

  for(i in 1:3){
    rels.rate$overall <- (numb.rels)/ diff(timewindow)
    rels.rate$men <- (numb.rels.men)/diff(timewindow)
    rels.rate$women <- (numb.rels.women)/diff(timewindow)
  }
  return(rels.rate)
}


# Relationship rate in time interval subdivided in different steps

# Relationship rate
rels.rate.int <- function(datalist = datalist,
                                  timewindow = c(0, 40), by=1){

  Rels.table <- datalist$rtable

  upper.limit <- ceiling(diff(timewindow)/by)

  interv.time <- round(seq.int(timewindow[1], timewindow[2], length.out = upper.limit), digits = 0)

  rels.int <- vector()
  rels.rate.int <- vector()


  for(i in 0:(upper.limit-2)){

    timewindow.int <- c(interv.time[1+i], interv.time[2+i])

    rels.numb <- nrow(Rels.table %>%
                             subset(FormTime <=timewindow.int[2] & FormTime >= timewindow.int[1]))


    rels.int <- c(rels.int,  rels.numb)


    rels.rate.int <- c(rels.rate.int, (rels.numb / diff(timewindow.int)))

  }

  return(rels.rate.int)
}


# Relationship ratio for men and women

rels.gender.ratio <- function(datalist = datalist,
                            timewindow = c(0, 40)){

  Rels.table <- datalist$rtable

  Rels.table.window <- Rels.table %>%
    subset(FormTime <=timewindow[2] & FormTime >= timewindow[1])

  numb.rels <- nrow(Rels.table.window)

  numb.rels.men <- length(unique(Rels.table.window$ID1)) # Gender 0 men & 1 women

  numb.rels.women <- length(unique(Rels.table.window$ID2))

  rels.ratio <- vector("list", length(c(0,1)))

  for(i in 1:2){
    rels.ratio$men <- (numb.rels.men/numb.rels)  #/diff(timewindow)
    rels.ratio$women <- (numb.rels.women/numb.rels)  #/diff(timewindow)
  }
  return(rels.ratio)
}



############################################
# 4. Slope of subtrees ~ size of subtrees ##
############################################

# fit.stat.slope
# library(apTreeshape)
fit.stat.slope <- function(tree=tree11){

  a <- as.phylo(tree)
  tips.labels <- a$tip.label
  numb.tips <- length(tips.labels)
  d <- rev(seq(from = 2, to = numb.tips, by = 1)) # sequence of sizes of subtrees (number of tips)

  tree1 <- as.treeshape(tree) # the tre must be a treeshape object
  s <- spectrum.treeshape(tree1) # sequence of number of subtrees

  # in d we read the number of tips for a given tree
  # in s we read the subtrees (trees) with corresponding (in s) tips

  # delete blank values
  nonzero.position = which(s != 0)
  s = s[nonzero.position]
  d = d[nonzero.position]

  reg <- lm(log(d) ~ log(s))
  cozf = coef(reg)

  # we expect the sizes of subtrees to decrease follwoing power-law hypothetically
  power.law.fit = function(x) exp(cozf[[1]] + cozf[[2]] * log(x))
  alpha = -cozf[[2]]
  R.square = summary(reg)$r.squared
  print(paste("Alpha =", round(alpha, 3)))
  print(paste("R square =", round(R.square, 3)))

  # # plot
  # plot(s ~ d, log = "xy", xlab = "Subtree size (log)", ylab = "Number of subtree (log)",
  #      col = 1, main = "Degree Distribution")
  # curve(power.law.fit, col = "red", add = T, n = length(d))

}

fit.stat.slope(tree=tree11)


## Tree imbalance

library(expoTree)
library(ape)
library(apTreeshape)

# one infection two infections after
id1 = c(seq(from=0,to=39, by=1))
par1 = c(-1,0,0,1,1,2,2,3,3,4,4,5,5,6,6,
       7,7,8,8,9,9,10,10,11,11,12,12,
       13,13,14,14,15,15,16,16,17,17,18,18,19)
t1 = rev(c(seq(from=1,to=40, by=1)))
epi1 <- list()
epi1$itimes <- t1
epi1$dtimes <- rep(0, length(t1))
epi1$id <- id1
epi1$parent <- par1
epitree1 <- epi2tree(epi1)
examp1 <- phylogenetictree.trend(tree = epitree1)
x1 = examp1$num.tree
y1 = examp1$size.tree
reg1 <- lm(log(y1) ~ log(x1))
cozf1 = coef(reg1)

# one infection one infection after

id2 = c(seq(from=0,to=39, by=1))
par2 = c(seq(from=-1, to=38, by=1))
t2 = rev(c(seq(from=1,to=40, by=1)))
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


# One infection many infections after more than 3
id3 = c(seq(from=0,to=39, by=1))
par3 = c(-1,0,0,0,1,1,1,1,1,4,4,1,5,5,6,
         7,7,8,8,9,9,10,10,5,0,5,5,
         3,3,7,14,7,15,16,3,17,17,3,18,5)
t3 = rev(c(seq(from=1,to=40, by=1)))
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



# Distances & tMRCA matrices

# function to compute differnces b/w elements

diff.ele <- function(v){
  w <- vector()
for(i in 1:(length(v)-1)){
  d <- v[i+1]-v[i]
  w <- c(w,d)
}
  return(w)
}

# Case 1: constant difference b/w elements

v1 = c(1,2,3,4,5,6,7,8,9,10)

d1 = diff.ele(v1)

# case 2: ascending differences

v2 = c(1,2,4,7,11,16,22,29,37,46)

d2 = diff.ele(v2)



# case 3: mix differences

v3 = c(1,2,4,5,7,11,13,17,18,19)

d3 = diff.ele(v3)


setwd("/home/david/google_drive/David_Project_drafts_2017/Simpact_Presentation/")

# 1 Function to get the MRCA time matrix

time.mrca.matrix <- function(tree = phylo.tree){

  Mytree <- phylo.tree


  # Use of mrca() and branching.time() functions from ape package

  # Symmatric metrix with tips which where entries are the internal nodes which represents the MRCA between two tips
  # from the phylogenetic tree

  Mytree.mrca.tips <- as.data.frame(mrca(Mytree, full = FALSE)) # MRCA for tips only >> full = FALSE


  # tips which represent individuals (sequences)
  inds <- Mytree$tip.label

  # branching time
  branch.time <- branching.times(Mytree)

  # Preparing the matrix of time of MRCA
  #######################################
  # Change the Mytree.mrca.tips matrix a bit to be easily handle for as igraph or network obejct
  # we have to make the diagonal 0
  n = as.matrix(Mytree.mrca.tips)
  mat.fun.tree <- function(n){
    diag(n) <- 0
    m <- n
    return(m)
  }

  m <- as.data.frame(mat.fun.tree(n))

  # Get internal nodes of the phyloegentic tree

  i.nodes <- NULL
  int.nodes <- function(m){
    for(i in 1:nrow(m)){
      for(j in 1:nrow(m)){
        i.nodes <- c(i.nodes,m[i,j])
      }
    }
    return(unique(i.nodes))
  }

  h = int.nodes(m)
  # remove the first element which is the 0 on diagonal
  g <- h[(2:length(h))] # these are really the internal nodes which represents the MRCA in the phylogenetic tree

  # MRCA times
  # each internal node has its branching time
  # thus, we have to replace the internal nodes in the matrix with their branching times respectively

  mrca.time.val <- as.data.frame(matrix(,nrow(m),nrow(m))) # initiate an n*n empty matrix (n tips)

  mrca.time.matrix <- function(g,m,branch.time){ # write a function which take as parms the internal nodes (g)
    # the matrix of internal nodes (MRCA) which is m and
    # the branching time which correspond with time where the two seq eveoled separelty
    k <- as.data.frame(branch.time)
    p <- k$branch.time

    for(l in 1:length(g)){
      for(i in 1:nrow(m)){
        for(j in 1:nrow(m)){
          if(m[i,j] == g[l]){
            mrca.time.val[i,j] <- p[l]# replace the number of the internal node with the time of branching
          }
        }
      }

    }
    return(mrca.time.val)
  }

  v = mrca.time.matrix(g,m,branch.time) # the diagonal is empty with entries NA

  # remove the NA in the diagonal and put 0 since no MRCA of a given sequence comparing at itself
  matrix.mrca.time.func <- function(v){
    diag(v) <- 0
    times.mrca.fine <- v
    return(times.mrca.fine)
  }

  names.inds <- names(m) # get the names of the tips
  mrca.times.final <- matrix.mrca.time.func(v)
  names(mrca.times.final) <- names.inds # put the names of the tips on the columns

  return(mrca.times.final)
}

# mrca matrix
library(ape)
phylo.tree <- read.nexus(file = "AnalyseAkaroreroRun.nex")
matrix.time <- time.mrca.matrix(tree = phylo.tree)


## 2. Function to get branch similarities

## Pairwise Distances from a Phylogenetic Tree

## cophenetic.phylo computes the pairwise distances between the pairs of tips
## from a phylogenetic tree using its branch lengths.
## dist.nodes does the same but between all nodes, internal and terminal, of the tree.
pair.similarity <- cophenetic.phylo(phylo.tree)

summary(phylo.tree)

num <- matrix.time[1:nrow(matrix.time),]

# Matrix elements sorted

matrix.elements <- function(M){
  vec0 <- vector()
  for(i in 1:nrow(M)){
    vec0 <- c(vec0,M[[i]])
  }

  veci <- sort(unique(vec0))
  vecf <- veci[-1]
  return(vecf)
}

matr.entries <- matrix.elements(matrix.time)

# Differences between matrix entries

diff.el.num <- sort(diff.ele(matr.entries))



# Fit with linear line the differences
X <- diff.el.num
hist(X, prob=TRUE)            # prob=TRUE for probabilities not counts
lines(density(X))             # add a density estimate with defaults
lines(density(X, adjust=2), lty="dotted")   # add another "smoother" density

dt <- diff.el.num
hist(X, prob=TRUE)
curve(dnorm(x, mean=mean(dt), sd=sd(dt)), add=TRUE)

# library(igraph)
net=graph.adjacency(as.matrix(matrix.time),mode="undirected",weighted=T,diag=FALSE)

stop.cut <- 1.9
net2 <- delete_edges(net, E(net)[weight>=stop.cut])
# get.all.shortest.paths(net2)
plot.igraph(net2)


matr.entries <- matrix.elements(matrix.time)
diff.el.num <- sort(diff.ele(matr.entries))
m <- mean(diff.el.num)

# graph.ls <- list()
# for(i in 1:length((E(net)$weight))){
#
#   if()
#
# }


mm <- m+m/(sd(diff.el.num)^2)

# > matr.entries
# [1] 0.001277305 0.574431157 0.599367627 0.736886814 1.299667274 1.311873478 1.773989454 1.884629179 1.970673869
# [10] 2.023170804 2.078535193 2.127039203 2.854873095 3.349539735 3.454737167 3.469222280 4.280282141 4.737336160
# [19] 5.266224504 6.760823423 7.641685171 8.681765011

stop.cut1 <- 0.736886814
net1 <- delete_edges(net, E(net)[weight>=stop.cut1])

stop.cut2 <- 1.773989454
net2 <- delete_edges(net, E(net)[weight>=stop.cut2])

stop.cut3 <- 2.023170804
net3 <- delete_edges(net, E(net)[weight>=stop.cut3])

net.uni <- (net1 %u% net2 %u% net3)

#
