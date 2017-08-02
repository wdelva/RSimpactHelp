# Phylogenetic trees comparison
##############################


# EX. 0

library(ape)

# same number of leaves is required for this
dist.topo(x, y = NULL, method = "PH85")
# This function computes the topological distance between two phylogenetic trees
# or among trees in a list (if y = NULL using different methods.

# The trees are always considered as unrooted.
# Note
#
# The geodesic distance of Billera et al. (2001) has been disabled:
# see the package distory on CRAN.
ta <- rtree(20)
tb <- rtree(20)
dist.topo(ta, ta) # 0
dist.topo(ta, tb) # unlikely to be 0


# EX. 1

# This function plots two trees face to face with the links if specified.
# It is possible to rotate the branches of each tree around the nodes by clicking.

# Usage

# cophyloplot(x, y, assoc=NULL, use.edge.length=FALSE,space=0,
#         length.line=1, gap=2, type="phylogram", rotate=FALSE,
#          col="red", show.tip.label=TRUE, font=3, ...)

#two random trees
tree1<-rtree(8) #random tree with 8 leaves
tree2<-rtree(4) #random tree with 4 leaves

#creation of the association matrix
association<-matrix(ncol=2, nrow=40)
association[,1]<-association[,2]<-tree2$tip.label

#plot
cophyloplot(tree1, tree2, assoc=association, length.line=4, space=28, gap=3)

#plot with rotations
## Not run:
cophyloplot(tree1, tree2, assoc=association, length.line=4, space=28, gap=3, rotate=TRUE)
## End(Not run)

# EX. 2

# ape

## S3 method for class 'phylo':
# all.equal(target, current, use.edge.length = TRUE,
#           use.tip.label = TRUE, index.return = FALSE,
#           tolerance = .Machine$double.eps ^ 0.5,
#           scale = NULL, ...)

### maybe the simplest example of two representations
### for the same rooted tree...:
t1 <- read.tree(text = "(a:1,b:1);")
t2 <- read.tree(text = "(b:1,a:1);")
all.equal(t1, t2)
### ... compare with this:
identical(t1, t2)
### one just slightly more complicated...:
t3 <- read.tree(text = "((a:1,b:1):1,c:2);")
t4 <- read.tree(text = "(c:2,(a:1,b:1):1);")
all.equal(t3, t4) # == all.equal.phylo(t3, t4)
### ... here we force the comparison as lists:
all.equal.list(t3, t4)

t5 <- read.tree(text = "(a:2,(c:1,b:1):1);")
### note that this does NOT return FALSE...:
all.equal(t3, t5)
### ... do this instead:
identical(all.equal(t3, t5), TRUE)

# EX. 3

library(distory)
#
# phylo.diff plots two trees side by side,
# highlighting edges unique to each tree in red.
#
# distinct.edges finds the edges present in the first argument
# not in the second.

phylo.diff(tree1,tree2)



# EX. 4

# same number of leaves is required for this

# dist.multiPhylo(x, method="geodesic", force.multi2di = FALSE, outgroup = NULL,
#                 convert.multifurcating = FALSE, use.random.resolution = FALSE, scale =
#                   NULL, verbose = FALSE)

# Computes the geodesic distance of a list of phylogenetic trees
# using a polynomial algorithm.

# Returns a distance matrix of class 'dist'
# representing the pairwise geodesic distances between all input trees.
# Keep in mind this distance matrix is not Euclidean. N/A values are
# provided in the case of an error in determining the distance.


t1 = rtree(8)
t2 = rtree(8)
dist.multiPhylo(c(t1,t2))

data(woodmouse)
otree <- root(nj(dist.dna(woodmouse)), "No305", resolve.root=TRUE)
breps = 250

trees <- boot.phylo(otree, woodmouse, B=breps, function(x)
  root(nj(dist.dna(x)), "No305", resolve.root=TRUE),trees=TRUE)

combined.trees <- c(list(otree), trees$trees)
tree.dists <- dist.multiPhylo(combined.trees)

mdres <- cmdscale(tree.dists, k=breps, add=TRUE)
plot(mdres$points[,1], mdres$points[,2], col = c("red", rep("black", breps)))
text(mdres$points[,1], mdres$points[,2], labels=1:(breps+1), cex=0.7, adj=c(0,2))


# Subtrees from a phylogenetic tree
####################################

library(ape)

# EX. 1

# subtrees()
# This function returns a list of all the subtrees of a phylogenetic tree.

# Value

# subtrees returns a list of trees of class "phylo" and returns invisibly for each subtree a list with the following components:
#
#   tip.label
# node.label
# Ntip
# Nnode

### Random tree with 12 leaves
phy<-rtree(12)
par(mfrow=c(4,3))
plot(phy, sub="Complete tree")


### Extract the subtrees
l<-subtrees(phy)

### plot all the subtrees
for (i in 1:11) plot(l[[i]], sub=paste("Node", l[[i]]$node.label[1]))
par(mfrow=c(1,1))
# like
# drop.tip()
# extract.clade() from ape, and cutreeshape() from apTreeshape packages

# spectrum.treeshape() from apTreeshape package
# This function returns a sequence containing the number of subtrees of size
# n, n-1, ..., 3, 2 where n is the size of the tree. The 'k'th element of the sequence
# is the number of subtrees of size n-k+1 in the tree, where n is the number of tips
# of the tree.


# EX. 2

library(apTreeshape)

# tipsubtree() from apTreeshape package
# Extract a subtree that contains pre-specified tip names or labels

## The universal tree of life provided in the data sets.
data(universal.treeshape)

## One might want to extract the tree containing the Animals, the Plants,
##      the Aquifex and the Microsporidia
tree1<-tipsubtree(universal.treeshape,tips=c("Animals", "Aquifex",
                                             "Microsporidia", "Plants"))
plot(universal.treeshape, tree1)

## Labels that do not appear in the tree are ignored
tree2<-tipsubtree(universal.treeshape,tips=c("Human", "Animals", "Aquifex",
                                             "Microsporidia", "Plants"))
plot(universal.treeshape, tree2)

tree3<-tipsubtree(universal.treeshape,tips=c(1,3,7), numeric=TRUE)
plot(universal.treeshape, tree3)

