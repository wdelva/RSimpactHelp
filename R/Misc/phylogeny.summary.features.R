## Features from phylogenetic trees: applied for phylogenies with differeent # of tips

# 1. Difference of lineages through time between two phylogenies - based on branch length

library(nLTT)

tree1 <- rtree(7)
tree2 <- rtree(9)

feature1 <- nltt_diff(tree1, tree2, distance_method = "abs") # important


# 2. Tree shape difference between two phylogenies - based on labelling scheme

library(treetop)

feature2 <- veclabeldistance(treelabels(rtree(10)), treelabels(rtree(11))) # important D2


# 3. Mean height of internal nodes & Maximum height

# The height of a node is the number of edges on the longest path from the node to a leaf.
# A leaf node will have a height of 0.

library(phytools)

tree <- rtree(15)

H <- nodeHeights(tree) # similar to node.depth.edgelength(tree) 


# It's clear from a casual inspection of the matrix that each parent node height (in the right column) 
# is represented twice and only twice. Thus, if we exclude the root node (zero height), 
# we can just take the mean of H[,1].

feature3 <- mean(sort(H[,1])[3:nrow(H)]) # important

library(phyloTop)

feature4 <- maxHeight(tree, normalise = FALSE)


# Node depths

# The depth of a node is the number of edges from the node to the tree's root node.
# A root node will have a depth of 0.

depth.features <- getDepths(tree) # depth of tips and nodes


# 4. Phylogenetic tree balace  features

library(phyloTop)

feature5 <- colless.phylo(tree, normalise = TRUE)

feature6 <- sackin.phylo(tree, normalise = FALSE)


