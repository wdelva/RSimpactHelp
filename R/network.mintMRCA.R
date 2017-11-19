#' Compute transmission network from tMRCA of tips of a phylogenetic tree
#'

setwd("/home/david/Dropbox/ANALYSIS_NOVEMBER_2017/testTreeNetwork/")
# Sampling dates in calender time
dates.Transform.NamedVector  <- function(dates=dates){
  dates.val <- 1977+40-as.numeric(dates$V2) # dates datalist$itable$population.simtime[1] - dates$V2 + 1977
  names(dates.val) <- as.character(dates$V1) # names are the names of the tips
  return(dates.val)
}


# (i) tree construction
system(paste("./FastTree  <", paste("A.Epidemic32.Sequences.gene.pol.fasta", sep = ""), paste(">A.Epidemic32.Sequences.gene.pol.fasta.tree", sep = "")))


# (ii) internal node calibration
tree.fasttree.32 <- read.tree(paste("A.Epidemic32.Sequences.gene.pol.fasta.tree", sep = ""))

samp.dates.32 <- read.csv("~/Dropbox/ANALYSIS_NOVEMBER_2017/testTreeNetwork/samplingtimes_seed_number_32.csv")

time.samp.32 <- dates.Transform.NamedVector(dates=samp.dates.32)

tree.tips.32 <- as.numeric(tree.fasttree.32$tip.label)

Ord.tree.dates.32 <- vector() # order the dates according to tips order in the tree
for(i in 1:length(tree.tips.32)){
  for(j in 1:length(time.samp.32)){
    if(tree.tips.32[i] == samp.dates.32$V1[j]){
      Ord.tree.dates.32 <- c(Ord.tree.dates.32, time.samp.32[j])
    }
  }
}

# calibrate internal nodes
dater.tree.32 <- dater(tree.fasttree.32,
                       Ord.tree.dates.32,
                       s = 3012,
                       omega0 = 0.00475) # s is the length of sequence
d=node.age(dater.tree.32)

save(dater.tree.32, file = paste("dated.tree.A.Epidemic32.Sequences.gene.pol.Rdata", sep = ""))

tree.A.seed32 <- get(load("dated.tree.A.Epidemic32.Sequences.gene.pol.Rdata"))

simpact.trans.net.32 <- simpact.trans.net[[32]]


graph.net.22 <- as.data.frame(simpact.trans.net.32)

graph.build.22 <- graph.net.22[,3:4]

graph.build.22[,2] <- as.character(graph.build.22[,2]) # donors
graph.build.22[,1] <- as.character(graph.build.22[,1]) # recipients
graph.22 = as.matrix(graph.build.22)
graph.f.22 = graph.edgelist(graph.22[,1:2])



network.mintMRCA <- function(tree = tree){

  # Symmetric matrix with tips which where entries are the internal nodes which represents the MRCA between two tips
  # from the phylogenetic tree

  Mytree.mrca.tips <- as.data.frame(mrca(tree, full = FALSE)) # MRCA for tips only >> full = FALSE


  # Preparing the matrix of time to the MRCA
  ##########################################
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

  g.sorted <- sort(g)
  node.age.ordered <- tree$Ti
  ind.ids <- tree$tip.label


  for (i in 1:length(ind.ids)) {
    id <- ind.ids[i]
    for (j in 1:nrow(m)) {
      par.i <-
      for (k in 1:nrow(m)) {

        m.i <- m[i,]

      }
    }

  }
#
#   # Link individuals with same MRCA expressed by same internal node
#   mrca.time.val <- as.data.frame(matrix(,nrow(m),nrow(m))) # initiate an n*n empty matrix (n tips)
#
#   for(l in 1:length(g.sorted)){
#     for(i in 1:nrow(m)){
#       for(j in 1:nrow(m)){
#         if(m[i,j] == g.sorted[l]){
#           mrca.time.val[i,j] <- node.age.ordered[l]
#         }
#       }
#     }
#
#   }
#
#   # remove the NA in the diagonal and put 0 since no MRCA of a given sequence comparing at itself
#   matrix.mrca.time.func <- function(mrca.time.val){
#     diag(mrca.time.val) <- 0
#     times.mrca.fine <- mrca.time.val
#     return(times.mrca.fine)
#   }
#
#   names.inds <- names(m) # get the names of the tips
#   mrca.times.final <- matrix.mrca.time.func(mrca.time.val)
#   names(mrca.times.final) <- names.inds # put the names of the tips on the columns
#
#   mrca.times.done <- mrca.times.final


}
