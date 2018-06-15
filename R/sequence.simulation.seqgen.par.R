#' Simulate sequence data of individuals in the transmissiion network using seq-gen tool
#'
#' @param dir.seq Direcotry where is the compiled seq-gen tool, sequences data may also be stored there
#' @param sub.dir.rename sub directory where the simulation outputs will be stored
#' @param simpact.trans.net Transmission network produced by \code{\link{transmission.network.builder()}}
#' @param seq.gen.tool Name of the file for compiled seq-gen
#' @param seeds.num  Seed number for reproducability
#' @param endpoint Only transmission events that took place before this point in simulation time
#' @param limitTransmEvents Minimum number of individuals in that transmission network (one transmission network per HIV seeding individual)
#' @param hiv.seq.file File containing seed HIV sequences in FASTA format, it might be in same repository as seq-gen
#' @return Files of sequences data, sampling dates, and transmission trees
#' @export


sequence.simulation.seqgen.par <- function(dir.seq = dirseqgen,
                                           sub.dir.rename = sub.dir.rename,
                                           simpact.trans.net = simpact.trans.net,
                                           seq.gen.tool = "seq-gen",
                                           seeds.num = seeds.num,
                                           endpoint = 40,
                                           limitTransmEvents = 7, # no less than 7
                                           hiv.seq.file = "hiv.seq.C.pol.j.fasta",
                                           clust = FALSE) # hiv.seq.file lodged in work.dir
{



  seedid <- seeds.num


  # source("R/transmission.network.builder.R")
  # source("R/epi2tree2.R")

  simpact.trans.net <- simpact.trans.net # transmission.network.builder(datalist = datalist, endpoint = endpoint)

  trans.net <- simpact.trans.net

  smallest.branches <- rep(NA, times = length(simpact.trans.net))
  for (list.element in 1:length(simpact.trans.net)){
    net.list <- simpact.trans.net[[list.element]]
    if(length(net.list$id) > 2){
      tree.tr <- epi2tree2(net.list)
      smallest.branch <- min(tree.tr$edge.length)
      smallest.branches[list.element] <- smallest.branch
    }
  }
  min.val <- min(smallest.branches, na.rm = TRUE)

  if(min.val <= 0){

    print("Some branch lengths are negative - someone transmitted infection after sampling - STOP")

  }else{

    TransmEventsCountVector <- vector()

    for(i in 1:length(trans.net)){
      trans.net.i.check <- as.data.frame(trans.net[[i]])
      if(nrow(trans.net.i.check)>=limitTransmEvents){
        TransmEventsCountVector <- c(TransmEventsCountVector, nrow(trans.net.i.check))
      }
    }


    ##### ######################
    #### Only One transmission network with at least limitTransmEvents each of them
    ############################
    if(length(TransmEventsCountVector)==1){



      num.trees <- vector() # ID of seeds # will be 0 because the transformation required to handle transmission tree from transmission network
      # constrained to rename IDs to -1, 0, 1, 2, ...
      num.i <- vector() # i_th seed in the list of seeds

      for(i in 1:length(trans.net)){

        tree.n <- trans.net[[i]] # transmission network for i^th seed

        if(nrow(as.data.frame(tree.n)) >= limitTransmEvents){

          # Construct transmission trees
          tree.i <- trans.network2tree(transnetwork = tree.n)
          num.trees <- c(num.trees,tree.n$id[1])
          num.i <- c(num.i,i)

          tree.j <- tree.i
          tree.j$tip.label <- paste(i,".", tree.j$tip.label, ".C", sep = "")

          # Save the transmission tree
          write.tree(tree.j, file = paste0(sub.dir.rename,"/tree.model1.seed_",i,".nwk"))

          # tr <- read.tree(file = paste("tree.model1.seed",i,".nwk", sep = "")


          # Keep sampling dates, add "A" on ID in order to handle all ID's as characters
          id.samplingtime <- as.data.frame(cbind(paste0(i,".",tree.n$id, ".C"), tree.n$dtimes)) # IDs and their sampling times in the transmission network

          write.csv(id.samplingtime,file=paste0(sub.dir.rename,"/samplingtimes_seed_",i,".csv"))

        }
      }

      IDs.transm <- num.i # vector of seeds with at least 2 transmissions

      p <- IDs.transm[1]

      resolved.combined.tree <- read.tree(file = paste0(sub.dir.rename,"/tree.model1.seed_",p,".nwk"))



      colitem <- read.csv(paste0(sub.dir.rename,"/samplingtimes_seed_",p,".csv"))


      write.csv(colitem, file=paste0(sub.dir.rename,"/samplingtimes.all.csv"))




      # Prepare to run the sequence simulation

      seq.rand <- 1 # first sequence

      n.tr <- 1 # number of transmission tree

      # # call the seed sequences - pool of viruses and rename the file
      file.copy(paste0(dir.seq,"/",hiv.seq.file),paste0(sub.dir.rename,"/seed.seq.bis.nwk"))

      # add the number of tree in the file and
      write(n.tr,file = paste0(sub.dir.rename,"/seed.seq.bis.nwk"), append = TRUE)  # n.tr
      # the tree, to prepare the file to simulate the evolution of the virus across the tree

      write.tree(resolved.combined.tree,file = paste0(sub.dir.rename,"/seed.seq.bis.nwk"), append = TRUE)

      file.rename(from = paste0(sub.dir.rename,"/seed.seq.bis.nwk"), to = paste0(sub.dir.rename,"/seed.seq.bis.sim.nwk"))

      out.seq.gen.file <- paste0(sub.dir.rename,"/C.Epidemic_seed.seq.bis.sim.nwk.fasta")

      in.seq.gen.file <- paste0(sub.dir.rename,"/seed.seq.bis.sim.nwk")

      system(paste0(paste0(dir.seq,"/"), paste0(paste0(seq.gen.tool)," -mGTR -f 0.3857, 0.1609, 0.2234, 0.2300  -a 0.9 -g 4 -i 0.5230  -r 2.9114, 12.5112, 1.2569, 0.8559, 12.9379, 1.0000 -s 0.00475  -n1 -k",seq.rand,"< ",in.seq.gen.file," -z",seedid," > ", out.seq.gen.file)))

    }


    ##### ######################
    #### Two or more transmissions networks with at least limitTransmEvents each of them
    ############################

    if(length(TransmEventsCountVector)>=2){ # at least two transmissions networks ith at least 6
      # transmissions each


      num.trees <- vector() # ID of seeds # will be 0 because the transformation required to handle transmission tree from transmission network
      # constrained to rename IDs to -1, 0, 1, 2, ...
      num.i <- vector() # i_th seed in the list of seeds

      for(i in 1:length(trans.net)){

        tree.n <- trans.net[[i]] # transmission network for i^th seed

        if(nrow(as.data.frame(tree.n)) >= limitTransmEvents){

          # Construct transmission trees
          tree.i <- trans.network2tree(transnetwork = tree.n)
          num.trees <- c(num.trees,tree.n$id[1])
          num.i <- c(num.i,i)

          tree.j <- tree.i
          tree.j$tip.label <- paste(i,".", tree.j$tip.label, ".C", sep = "")

          # Save the transmission tree
          write.tree(tree.j, file = paste0(sub.dir.rename,"/tree.model1.seed_",i,".nwk"))

          # tr <- read.tree(file = paste("tree.model1.seed",i,".nwk", sep = "")


          # Keep sampling dates, add "A" on ID in order to handle all ID's as characters
          id.samplingtime <- as.data.frame(cbind(paste0(i,".",tree.n$id, ".C"), tree.n$dtimes)) # IDs and their sampling times in the transmission network

          write.csv(id.samplingtime,file=paste0(sub.dir.rename,"/samplingtimes_seed_",i,".csv"))

        }
      }

      IDs.transm <- num.i # vector of seeds with at least 2 transmissions



      ### Binding all transmission trees together ###
      ###############################################

      # Make a list of transmission trees

      trees <- list() # list of all transmission trees

      for(j in 1:length(IDs.transm)){
        p <- IDs.transm[j]
        tr <- read.tree(file = paste0(sub.dir.rename,"/tree.model1.seed_",p,".nwk"))
        trees[[j]] <- tr
      }

      class(trees)<-"multiPhylo"

      # print(trees,details=TRUE)

      # Function to bind all transmission trees in the list above

      bind.trees<-function(trees){
        if(length(trees)==2) return(bind.tree(trees[[1]],trees[[2]], where = "root", position = 0))
        else {
          trees<-c(bind.tree(trees[[1]],trees[[2]]),
                   if(length(trees)>2) trees[3:length(trees)])
          trees<-bind.trees(trees) ## this is the recursive part
          return(trees)
        }
      }

      combined.tree<-bind.trees(trees) # This is a polytomy tree

      # is.binary.tree(combined.tree)

      # resolve polytomies of the combined transmission tree

      resolved.combined.tree <- multi2di(combined.tree) # resolve polytomies

      # is.binary.tree(resolved.combined.tree)


      ### Bind sampling dates ###
      ###########################

      for (i in IDs.transm){
        if(i==IDs.transm[1]){
          colitem <- read.csv(paste0(sub.dir.rename,"/samplingtimes_seed_",i,".csv"))
        }
        else{

          ritem <- read.csv(paste0(sub.dir.rename,"/samplingtimes_seed_",i,".csv"))
          colitem <- rbind(colitem, ritem)
        }
      }

      write.csv(colitem, file=paste0(sub.dir.rename,"/samplingtimes.all.csv"))


      # Prepare to run the sequence simulation

      seq.rand <- 1 # first sequence

      n.tr <- 1 # number of transmission tree

      # # call the seed sequences - pool of viruses and rename the file
      file.copy(paste0(dir.seq,"/",hiv.seq.file),paste0(sub.dir.rename,"/seed.seq.bis.nwk"))

      # add the number of tree in the file and
      write(n.tr,file = paste0(sub.dir.rename,"/seed.seq.bis.nwk"), append = TRUE)  # n.tr
      # the tree, to prepare the file to simulate the evolution of the virus across the tree

      write.tree(resolved.combined.tree,file = paste0(sub.dir.rename,"/seed.seq.bis.nwk"), append = TRUE)

      file.rename(from = paste0(sub.dir.rename,"/seed.seq.bis.nwk"), to = paste0(sub.dir.rename,"/seed.seq.bis.sim.nwk"))

      out.seq.gen.file <- paste0(sub.dir.rename,"/C.Epidemic_seed.seq.bis.sim.nwk.fasta") # will be input for tree construction
      in.seq.gen.file <- paste0(sub.dir.rename,"/seed.seq.bis.sim.nwk")

      system(paste0(paste0(dir.seq,"/"), paste0(paste0(seq.gen.tool)," -mGTR -f 0.3857, 0.1609, 0.2234, 0.2300  -a 0.9 -g 4 -i 0.5230  -r 2.9114, 12.5112, 1.2569, 0.8559, 12.9379, 1.0000 -s 0.00475  -n1 -k",seq.rand,"< ",in.seq.gen.file," -z",seeds.num," > ", out.seq.gen.file)))

    }

  }

  print("Sequence simulation finished!")

}

# sequence.simulation.seqgen.par to merge

