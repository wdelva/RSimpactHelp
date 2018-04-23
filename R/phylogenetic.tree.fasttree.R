#' Construct and calibrate a phylogenetic tree from sequence data
#'
#' @param dir.tree Direcotry where we find the FastTree tool
#' @param dir.seq Direcotry where we find sequences and sampling times files, it may be in the folder where simulation of sequences by \code{\link{sequence.simulation.seqgen()}} taken place
#' @param fasttree.tool Name of the compiled file of FastTree
#' @param calendar.dates File containing named sampling calendar times for each sequence
#' @param simseqfile File containing sequences
#' @param count.start Calendar year when the simulation started
#' @param endsim Number of years when the simulation was done
#' @return A time-stamped phylogenetic  tree with annotation of internal nodes dates
#' @import ape
#' @import treedater
#' @export

phylogenetic.tree.fasttree <- function(dir.tree = dirfasttree,
                                       dir.seq = dirseqgen,
                                       fasttree.tool = "FastTree",
                                       calendar.dates = "dates.csv",
                                       simseqfile = "C.Epidemic.sequence.fasta",
                                       count.start = 1977,
                                       endsim = 40){

  # Function to tranform dates in named vector to be handled by treedater

  dates.Transform.NamedVector  <- function(dates=dates){

    dates.val <- endsim - dates$V2 + count.start # dates
    names(dates.val) <- as.character(dates$V1) # names are the names of the tips

    return(dates.val)
  }

  # 4.1. Construct phylogenetic tree

  # Many tools to build the trees: within R like ape, phangorn, and outside compiled tools like iq-tree and FastTree


  out.fast.tree.file <- paste0(dir.tree,"/", simseqfile ,".tree")

  print("Start construction of the phylogenetic tree with FastTree")

  system(paste0(dir.tree,"/",paste0(paste(fasttree.tool, "-gtr")," -nt < ", paste0(dir.seq,"/",simseqfile), paste0("> ", out.fast.tree.file))))


  print("End of construction of the phylogenetic tree  with FastTree")



  samp.dates <- read.csv(paste0(dirseqgen, "/", calendar.dates))

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


  print("Start of internal nodes calibration for the phylogenetic tree with treedater")

  # 4.2. Calibrate the phylogenetic tree

  # Use of library(treedater) to calibrate internal nodes

  dater.tree <- dater(tree.const, Ord.tree.dates, s = 3000) # s is the length of sequence

  write.tree(dater.tree, file = "calibrated.tree.tree")

  print("End of internal nodes calibration for the phylogenetic tree")


  print("Calibrated phylogenetic tree is saved in the working directory")

  return(dater.tree)

}
