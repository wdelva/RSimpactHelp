#' Construct phylogenetic tree with FastTree and calibrate same tree with treedater from sequence alignment
#'
#' @param dir.tree Direcotry where we find the FastTree tool
#' @param sub.dir.rename sub directory where we find sequences and sampling times files and where the simulation outputs will be stored
#' @param fasttree.tool Name of the compiled file of FastTree
#' @param calendar.dates File containing named sampling calendar times for each sequence
#' @param simseqfile File containing sequences
#' @param count.start Calendar year when the simulation started
#' @param endsim Number of years when the simulation was done
#' @param clust Logic parameter TRUE when running simulation on cluster and FALSE on personal computer
#' @return A time-stamped phylogenetic  tree with annotation of internal nodes dates
#' @import ape
#' @import treedater
#' @export


phylogenetic.tree.fasttree.par <- function(dir.tree = dirfasttree,
                                           sub.dir.rename = sub.dir.rename,
                                           fasttree.tool = "FastTree",
                                           calendar.dates = "samplingtimes.all.csv",
                                           simseqfile = "C.Epidemic_seed.seq.bis.sim.nwk.fasta",
                                           count.start = 1977,
                                           endsim = 40,
                                           clust = TRUE){

  # Function to tranform dates in named vector to be handled by treedater

  dates.Transform.NamedVector  <- function(dates=dates){

    dates.val <- endsim - dates$V2 + count.start # dates
    names(dates.val) <- as.character(dates$V1) # names are the names of the tips

    return(dates.val)
  }

  # 4.1. Construct phylogenetic tree

  # Many tools to build the trees: within R like ape, phangorn, and outside compiled tools like iq-tree and FastTree


  out.fast.tree.file <- paste0(sub.dir.rename,"/", simseqfile ,".nwk")

  print("Start construction of the phylogenetic tree with FastTree")

  if(clust==TRUE){

    system(paste0(paste0(paste(fasttree.tool, "-gtr")," -nt < ", paste0(sub.dir.rename,"/",simseqfile), paste0("> ", out.fast.tree.file))))

  }else{

    system(paste0(dir.tree,"/",paste0(paste(fasttree.tool, "-gtr")," -nt < ", paste0(sub.dir.rename,"/",simseqfile), paste0("> ", out.fast.tree.file))))

  }


  print("End of construction of the phylogenetic tree  with FastTree")


  samp.dates <- read.csv(paste0(sub.dir.rename, "/", calendar.dates))

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

  dater.tree <- treedater::dater(tree.const, Ord.tree.dates, s = 3000) # s is the length of sequence

  write.tree(dater.tree, file = paste0(sub.dir.rename, paste0("/calibrated.tree.",simseqfile,".tree")))

  print("End of internal nodes calibration for the phylogenetic tree")


  print("Calibrated phylogenetic tree is saved in the working directory")

  return(dater.tree)

}


