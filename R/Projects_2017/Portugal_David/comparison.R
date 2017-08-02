# Phylogeneti trees comparison


# EX. 1

library(ape)

# This function plots two trees face to face with the links if specified.
# It is possible to rotate the branches of each tree around the nodes by clicking.

# Usage

# cophyloplot(x, y, assoc=NULL, use.edge.length=FALSE,space=0,
#         length.line=1, gap=2, type="phylogram", rotate=FALSE,
#          col="red", show.tip.label=TRUE, font=3, ...)

#two random trees
tree1<-rtree(8) #random tree with 40 leaves
tree2<-rtree(4) #random tree with 20 leaves

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
