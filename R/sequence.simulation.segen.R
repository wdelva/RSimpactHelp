#' Simulate sequences data of individual in the transmissiion network using seq-gen
#'
#' @param dir Direcotry where the simulations of sequence will be performed, thre might be compiled seq-gen tool
#' @param datalist The datalist that is produced by \code{\link{readthedata()}}
#' @param seeds.num  Seed number for reproducability
#' @param endopint Only transmission events that took place before this point in simulation time
#' @param seed.file File containing seed HIV sequences


sequence.simulation.seqgen <- function(dir = dir, datalist = datalist, seeds.num = 123, endpoint = 40, seed.file = "hiv.seq.A.pol.j.fasta"){



  simpact.trans.net <- transmNetworkBuilder.diff2(datalist = datalist, endpoint = endpoint)

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
    print("Some branch lengths are negative - someone transmitted infection after sampling")
  }else{
    current.dir <- getwd()

    dir <- setwd(dir)

    seed <- seeds.num

    trans.net <- simpact.trans.net

    num.i <- vector() # i_th seed in the list of seeds
    for(i in 1:length(trans.net)){

      tree.n <- trans.net[[i]] # transmission network for i^th seed

      if(nrow(as.data.frame(tree.n)) >= 3){
        tree.i <- trans.network2tree(transnetwork = tree.n)
        #num.trees <- c(num.trees,tree.n$id[1])
        num.i <- c(num.i,i)

        # Save the transmission tree
        write.tree(tree.i, file = paste("tree.model1.seed",i,".nwk", sep = ""))

        tr <- read.tree(file = paste("tree.model1.seed",i,".nwk", sep = ""))

        # # count number of trees generated (normally one)
        # numb.tr <- function(tree=tree){
        #     if(length(tree) == 4){
        #         return(1)
        #     }else{
        #         return(length(tree))
        #     }
        # }
        #
        # # Simulate sequences, this require compiled toold seq-gen
        # # random sequence which will be simulated is choosen in the pool
        # seq.rand <- sample(1:30,1) # chose one sequence in first 30's in the seed sequence pool named "seed.seq.fasta"
        seq.rand <- 1
        #
        # tr <- read.tree(file = paste("tree.model1.seed",i,".nwk", sep = ""))
        # n.tr <- numb.tr(tree = tr)

        n.tr <- 1

        seed.id <- tree.n$id[1]

        # call the seed sequences and rename the file; paste("hiv.seq.A.pol.j.fasta", sep = "")
        file.copy(seed.file,paste("hiv.seq.A.pol.j.fasta",i,".nwk", sep = ""))
        # add the number of tree in the file and
        write(n.tr,file = paste("hiv.seq.A.pol.j.fasta",i,".nwk", sep = ""), append = TRUE) # 1 is the number of tree across which we simulate the sequences
        # the tree, to prepare the file to simulate the evolution of the virus across the tree
        write.tree(tr,file = paste("hiv.seq.A.pol.j.fasta",i,".nwk", sep = ""), append = TRUE)
        file.rename(from = paste("hiv.seq.A.pol.j.fasta",i,".nwk", sep = ""), to = paste("hiv.seq.A.pol.j.fasta",i,"seed",i,".nwk", sep = ""))

        system(paste("./seq-gen -mGTR -f 0.3857, 0.1609, 0.2234, 0.2300 -a 0.9410 -i 0.80 -g 4 -r 2.2228, 10.7771, 1.0675, 1.2966, 12.6824, 1.0000 -s 0.00475 -n1 -k 1 <hiv.seq.A.pol.j.fasta",i,"seed",i,".nwk  -z",seed," >A.Epidemic",i,".Sequences.gene.pol.fasta", sep = ""))


        # a: shape parameter of Gamma > Gamma Rate Heterogeneity
        # g: category of Gamma > Discrete Gamma Rate Heterogeneity
        # r: rate matrix
        # z: seed

        # Keep sampling dates
        id.samplingtime <- as.data.frame(cbind(tree.n$id, tree.n$dtimes)) # IDs and their samling times in the transmission network

        write.csv(id.samplingtime,file=paste("samplingtimes_seed_number_",i,".csv", sep = ""))

      }
    }

    setwd(current.dir)
  }


}
