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


# EX. 5

library(phangorn)

# treedist computes different tree distance methods and RF.dist the Robinson-Foulds
# or symmetric distance. The Robinson-Foulds distance only depends on the toplogy of the trees.
# If edge weights should be considered wRF.dist calculates
# the weighted RF distance (Robinson & Foulds 1981). and KF.dist calculates
# the branch score distance (Kuhner & Felsenstein 1994).
# path.dist computes the path difference metric as described in Steel and Penny 1993).
# sprdist computes the approximate SPR distance
# (Oliveira Martins et al. 2008, de Oliveira Martins 2016).

tree1 = rtree(21)
tree2 = rtree(21)

treedist(tree1, tree2, check.labels = TRUE)

sprdist(tree1, tree2)

SPR.dist(tree1, tree2)

RF.dist(tree1, tree2, normalize = FALSE, check.labels = TRUE,
        rooted = FALSE)

wRF.dist(tree1, tree2, normalize = FALSE, check.labels = TRUE,
         rooted = FALSE)

KF.dist(tree1, tree2 , check.labels = TRUE, rooted = FALSE)

path.dist(tree1, tree2, check.labels = TRUE, use.weight = FALSE)



# Subtrees from a phylogenetic tree, prune a tree
##################################################

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

# EX. 2

tree<-rtree(10)
write.tree(tree)

# Now, say we want to keep the species t2, t4, t6, t8, and t10
# in our pruned tree, we just put these tip names into a vector:
species<-c("t2","t4","t6","t8","t10")

pruned.tree<-drop.tip(tree,tree$tip.label[-match(species, tree$tip.label)])
# Or
pruned.tree<-drop.tip(tree, setdiff(tree$tip.label, species))

write.tree(pruned.tree)

# spectrum.treeshape() from apTreeshape package
# This function returns a sequence containing the number of subtrees of size
# n, n-1, ..., 3, 2 where n is the size of the tree. The 'k'th element of the sequence
# is the number of subtrees of size n-k+1 in the tree, where n is the number of tips
# of the tree.


# EX. 3

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


# Construct a phylogenetic tree
################################

# EX. 1

eg.phy <- list()

m<-matrix(nrow=18,ncol=2)
m[1,]<-c(11,12)
m[2,]<-c(12,13)
m[3,]<-c(13,14)
m[4,]<-c(14,15)
m[5,]<-c(15,1)
m[6,]<-c(15,2)
m[7,]<-c(14,16)
m[8,]<-c(16,3)
m[9,]<-c(16,4)
m[10,]<-c(13,17)
m[11,]<-c(17,18)
m[12,]<-c(18,5)
m[13,]<-c(18,6)
m[14,]<-c(17,7)
m[15,]<-c(12,8)
m[16,]<-c(11,19)
m[17,]<-c(19,9)
m[18,]<-c(19,10)

eg.phy$edge <- m

eg.phy$tip.label <- c("t9", "t4", "t5", "t10","t7", "t1", "t6","t8","t2", "t3")

eg.phy$edge.length <- c(0.76947, 0.33053, 0.88915, 0.41914, 0.09538,
                        0.92463, 0.80391, 0.62641, 0.85228, 0.11455,
                        0.39461, 0.21585, 0.33571, 0.39673, 0.82935,
                        0.21624, 0.91556, 0.59869)

eg.phy$Nnode <- length(unique(m[,1]))

class(eg.phy) <- "phylo"

reorder(eg.phy, order = "cladewise", index.only = FALSE)

# attr(eg.phy, "class") <- "phylo"
# attr(eg.phy, "order") <- "cladewise"

str(eg.phy)

tp <- vector()
for(i in 1:10){
  d <- paste(i)
  tp <- c(tp,d)
}


# EX. 2

library(apTreeshape)

# treeshape(nodes, names)
# Arguments
#
# nodes
# nodes is a n*2 matrix containing the node structure of the tree.
# names
# names is a vector which contains the names of the tips.


## Nodes will define the nodes of a five tips tree
nodes<-matrix(nrow=4,ncol=2)
nodes[1,]<-c(-5,-4)
nodes[2,]<-c(1,-1)
nodes[3,]<-c(-3,2)
nodes[4,]<-c(-2,3)

nodes1[5,]<-c(3,5)

## Now we can build the tree and plot it.
tree1<-treeshape(nodes)
plot(tree1)

## Computation of the sackin index for the tree :
sackin(tree1)

## Label will define the names of the tips
label=c("a", "b", "c", "d", "e")
tree2<-treeshape(nodes, label)
plot(tree1, tree2)


x1 <- c(0,1,1,1)
x2 <- c(1,0,2,1)
x3 <- c(1,2,0,2)
x4 <- c(1,1,2,0)

y = cbind(x1,x2,x3,x4)

x <- as.matrix(y)

d = hclust(dist(x))

plot(d)

g = as.phylo(d)

str(g)

plot(g)

# from trannet

d1 <- epi1$itimes
d2 <- epi1$dtimes

dt <- d2-d1


