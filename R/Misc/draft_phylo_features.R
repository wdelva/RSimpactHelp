

# Ape
library(ape)
tree <- rtree(4)

node.depth(tree)
node.depth.edgelength(tree)
node.height(tree)


# Trick to get the mean height of internal nodes
library(phytools)

# tree<-rtree(10)
plot(tree)
add.scale.bar(length=0.5)

H<-nodeHeights(tree) # similar to node.depth.edgelength(tree) 


# It's clear from a casual inspection of the matrix that each parent node height (in the right column) 
# is represented twice and only twice. Thus, if we exclude the root node (zero height), 
# we can just take the mean of H[,1].

hbar1<-mean(sort(H[,1])[3:nrow(H)]) # important
hbar1

t <- read.tree("~/RSimpactHelp/R/Projects_2017/Example3PaperSimpact2/calibratedTree_10.nwk")


library(phyloTop)

cherries(tree)

colless.phylo(tree, normalise = TRUE)

sackin.phylo(tree, normalise = FALSE)

getDepths(tree) # depth of tips and nodes

maxHeight(tree, normalise = FALSE)


library(nLTT)

tree1 <- rtree(7)
tree2 <- rtree(9)
nltt_diff(tree1, tree2, distance_method = "abs") # important

# nltt_diff_exact(tree1, tree2, distance_method = "abs", ignore_stem = FALSE)

nltt_lines(tree)


library(apTreeshape)

x <- as.treeshape(tree)

shape.statistic(x)

spectrum.treeshape(x)

summary.treeshape(x)


library(treetop)

treelabels(tree)           


labeldistance(treelabels(rtree(10)), treelabels(rtree(120)))

distunlab(rtree(10), rtree(120))

veclabeldistance(treelabels(rtree(10)), treelabels(rtree(120))) # important D2

dd1=multiDistUnlab(rmtree(10,12)) # D1

dd2=vecmultiDistUnlab(rmtree(10,12)) # D2     


            
