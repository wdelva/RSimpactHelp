#' COnstruct and calibrate a phylogenetic tree from sequence data
#'
#' @param dir.tree Direcotry where we find the FastTree tool
#' @param dir.seq Direcotry wherewe find sequences and sampling times files
#' @param calendar.dates File containing named sampling times for each sequence
#' @param simseqfile File containing sequences
#' @param seeds.num  Seed number for reproducability
#' @param endsim End of simulation # place holder for sequence from simpact when assigning real calendar time

phylogenetic.tree.fasttree <- function(dir.seq = dirseqgen,
                                       dir.tree = dirfasttree,
                                       calendar.dates = "dates.csv",
                                       simseqfile = "C.Epidemic.sequence.fasta",
                                       endsim = 40){

  # Function to tranform dates in named vector to be handled by treedater

  dates.Transform.NamedVector  <- function(dates=dates){

    dates.val <- endsim - dates$V2 + 1977 # dates
    names(dates.val) <- as.character(dates$V1) # names are the names of the tips

    return(dates.val)
  }

  # 4.1. Construct phylogenetic tree

  # Many tools to build the trees: within R like ape, phangorn, and outside compiled tools like iq-tree and FastTree


  out.fast.tree.file <- paste0(dir.tree,"/", simseqfile ,".tree")

  system(paste0(dir.tree,"/",paste0("FastTree -gtr -nt < ", paste0(dir.seq,"/",simseqfile), paste0("> ", out.fast.tree.file))))



  samp.dates <- read.csv(paste0(calendar.dates))

  tree.const <- read.tree(out.fast.tree.file)

  time.samp <- dates.Transform.NamedVector(dates=samp.dates) # name the dates


  ### Match dates and phylogenetic tree leaves ###
  ################################################

  tree.tips <- as.character(tree.const$tip.label)

  Ord.tree.dates <- vector()

  for(i in 1:length(tree.tips)){
    for(j in 1:length(time.samp)){
      if(tree.tips[i] == samp.dates$V1[j]){
        Ord.tree.dates <- c(Ord.tree.dates, time.samp[j])
      }
    }
  }

  # 4.2. Calibrate the phylogenetic tree

  # Use of library(treedater) to calibrate internal nodes

  dater.tree <- dater(tree.const, Ord.tree.dates, s = 3000) # s is the length of sequence

  write.tree(dater.tree, file = "calibrated.tree.tree")

  print("End of time-stamped phylogenetic tree")

  return(dater.tree)

}
