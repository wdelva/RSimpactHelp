#' Simulate evolutionary sequences accross transmission networks
#'
#' By considering past demographic events of viral population
#' phylogenetic trees of virus allow to simulate sequences at infection and sampling times
#' Within the working directory must be ms and seq-gen compiled tools and a file of alignment of seeds sequences
#' Inputs is a transmission networks with viral load, rate of change of viral load, and time from infection to sampling
#'
#'
#' Considering the within host dynamics, simulations of viruses
#'  hisotry are done per each individual by Wrish-Fisher model by considering,
#'  individual j is infected by one virus from a pool of individual i,
#'  at sampling times, a consensus sequence for each individual is obtained from
#'  a pool of viruses of that individual

############################################################
# Remember to use rates and frquencies of the chosen gene! #
############################################################

# Transmission network in a raw form with suppllement data includes viral load, CD4, etc.
seed=123

source("/home/david/RSimpactHelp/R/Projects_2017/transNetworkPastDemographic.R")

simpact.output.raw <- transNetworkPastDemographic(datalist=master.datalist, endpoint=40)

# Note: will consider transmission network with at least 3 individuals

for(v in 1:length(simpact.output.raw)){

  seed2.dat <- as.data.frame(simpact.output.raw[[v]])

  if(nrow(seed2.dat) >= 3){

    dat.process <- seed2.dat


    ## Within host dynamics per individual: ms > produce coalescent trees
    ## When sampled, each individual, many copies are taken and these copies have different demographics events
    ## When an individual i transmit an infection to a=individual j, we assume one virus seed the infection in new host
    for(i in 1:length(seed2.dat$rec.IDs)){

      don <- seed2.dat$don.IDs[i] #seed2.dat$parent[i] # further we do not consider -1 (universal infector)
      rec <- seed2.dat$rec.IDs[i] #seed2.dat$id[i]

      # Viral load at infection times and sampling times

      vl.inf <- 10*round(seed2.dat$infec.vload[i]) # viral load of the donors at time point of infecting someone else
      vl.samp <- 10*round(seed2.dat$samp.vload[i]) # viral load (of the recipients) at sampling time
      t.inf <- seed2.dat$t.infec[i] # infection time for recipients, time interval from when the current donor j
      # got the infection from a donor i (or seeding event) until when current donor j transmits to recipient k
      t.samp <- seed2.dat$t.samp[i] #sampling time for recipients, time when recipient k has been sampled (diadnosed or removed by death)

      past.demographic.inf <- c(round(seed2.dat$rate.2back.infec.vload[i], digits = 3), round(seed2.dat$rate.1back.infec.vload[i], digits = 3)) #, 0.01) # a vector of change rates of viral loads from infection time to time to transmit the infection

      past.demographic.samp <- c(round(seed2.dat$rate.2back.samp.vload[i], digits = 3), round(seed2.dat$rate.1back.samp.vload[i], digits = 3)) #, 0.01) # or sampling time

      # Note: just two values in the vectors

      seeds <- c(1+seed, 2+seed, 3+seed)

      seeds1 <- seeds[1]
      seeds2 <- seeds[2]
      seeds3 <- seeds[3]

      # Viral load models and coalescent: chosen "model 3"


      # Models 3: Consider viral load without islands models with past demographic events which
      #           are higlihthed by the viral load variation with different infection stages - growth rate set to zero
      # system('./msa 30 1 -T -t 4.0 -eN 0.2 .02')
      # specifies that the population size was constant at size N 0 from the present back to
      # time 0.2*4N0, and farther back in time the population size was 0.02*4N0  .
      # The population size change was instantaneous and occurred at time 0.2*4N0
      # generations before the present.

      # Note: growth rate set to zero is the side-effect of -eN and -en commands in ms
      #       but is nour case, this coincide of simpact modelling approach of the viral load

      topaste.inf <- paste("./msa", vl.inf," 1 -T -t",t.inf)
      topaste.samp <- paste("./msa", vl.samp," 1 -T -t",t.samp)
      f.inf <- paste(paste(paste(topaste.inf,paste("-eN",
                                                   paste(past.demographic.inf[1], past.demographic.inf[2], collapse = ",")), collapse = ""),
                           paste("-seeds",paste(seeds1,seeds2,seeds3, collapse = ","), collapse = ","), collapse = ""),
                     ">tree_at_from_",don,"_to_",rec,"_net_",v,".nwk", sep = "")  # transmission

      f.samp <- paste(paste(paste(topaste.inf,paste("-eN", paste(past.demographic.inf[1],
                                                                 past.demographic.inf[2], collapse = ",")), collapse = ""),
                            paste("-seeds",paste(seeds1,seeds2,seeds3, collapse = ","), collapse = ","), collapse = ""),
                      ">tree_at_for_",rec,"_samp_net_",v,".nwk", sep = "")  # sampling

      system(f.inf[1])
      system(f.samp[1])

    }



    ## Generate sequences for each individual at sampling time ##
    #############################################################

    ## The seed or (and) its fist recipients are responsible of the subsequent infections


    # (i) Transmissions by the seed [1] in XXX$id

    seed.id <- seed2.dat$rec.IDs[1] # Seed individual [1] in XXX$id
    # rec.seed <- seed2.dat$id[2] # First transmission by the seed

    rec.ids <- seed2.dat$rec.IDs
    par.ids <- seed2.dat$don.IDs

    for(i in 2:length(seed2.dat$rec.IDs)){

      rec.id <- rec.ids[i]
      par.id <- par.ids[i]
      rec.seed <- rec.id # recipient from seed donor

      if(par.id == seed.id){ # this means this individual is infected by the seed element

        # donor coalescent tree at transmission times by the seed
        seed.tree.trans <- read.tree(paste("tree_at_from_",seed.id,"_to_",rec.seed,"_net_",v,".nwk", sep = ""))

        # compute number of trees in the file
        numb.tr <- function(tree=tree){
          if(length(tree) == 4){
            return(1)
          }else{
            return(length(tree))
          }
        }
        numb.trees.seed <- numb.tr(tree=seed.tree.trans) # count number of trees generated (normally one)

        file.copy(paste("HIV.Pol.gene.fasta", sep = ""),paste("HIV.Pol.gene.bis.nwk", sep = "")) # call the seed sequences - pool of viruses

        write(numb.trees.seed,file = "HIV.Pol.gene.bis.nwk", append = TRUE) # add the number of tree in the file and
        write.tree(seed.tree.trans,file = "HIV.Pol.gene.bis.nwk", append = TRUE) # the tree, to prepare the file to simulate the evolution of the virus across the tree
        file.rename(from = "HIV.Pol.gene.bis.nwk", to = paste("HIV.Pol.gene.bis",seed.id,".nwk", sep = ""))

        seq.rand <- sample(1:10,1) # random number correponds to random sequence chosed in the pool of viruses

        # Sequence of the seed at the first tramsmission event
        system(paste("./seq-gen -mGTR -f 0.3887782, 0.1646746, 0.2277556, 0.2187915 -a 0.9 -g 4 -r 0.00125, 0.00125, 0.00125, 0.00125, 0.00125, 0.00125 -n1 -k",seq.rand,"< HIV.Pol.gene.bis",seed.id,".nwk -z", seed,"> sequence_from_",seed.id,"_to_",rec.seed,"_net_",v,".fasta",sep = ""))

        # a: shape parameter of Gamma > Gamma Rate Heterogeneity
        # g: category of Gamma > Discrete Gamma Rate Heterogeneity
        # r: rate matrix
        # z: seed
      }
    }




    # (ii) Sampling of the seed [1] in XXX$id

    # donor coalescent tree at sampling time (by the seed)
    seed.tree.samp <- read.tree(paste("tree_at_for_",seed.id,"_samp_net_",v,".nwk", sep = "")) # tree at sampling time of the seed
    numb.tree.samp <- numb.tr(tree=seed.tree.samp) # count number of trees generated (normally one)

    file.copy(paste("HIV.Pol.gene.fasta", sep = ""),paste("HIV.Pol.gene.bis.nwk", sep = "")) # call the seed sequences - pool of viruses

    write(numb.tree.samp,file = "HIV.Pol.gene.bis.nwk", append = TRUE) # add the number of tree in the file and
    write.tree(seed.tree.samp,file = "HIV.Pol.gene.bis.nwk", append = TRUE) # the tree, to prepare the file to simulate the evolution of the virus across the tree
    file.rename(from = "HIV.Pol.gene.bis.nwk", to = paste("HIV.Pol.gene.bis",seed.id,".nwk", sep = ""))

    seq.rand <- sample(1:10,1) # random number correponds to random sequence chosed in the pool of viruses

    # Sequence of the seed at the sampling (diagnosis/removal) event

    system(paste("./seq-gen -mGTR -f 0.3887782, 0.1646746, 0.2277556, 0.2187915 -a 0.9 -g 4 -r 0.00125, 0.00125, 0.00125, 0.00125, 0.00125, 0.00125 -n1 -k",seq.rand,"< HIV.Pol.gene.bis",seed.id,".nwk -z",seed,"> sequence_at_samp_",seed.id,"_net_",v,".fasta",sep = ""))


    # (iii) Sampling of recipients from the seed donor in XXX$id[2]

    rec.ids <- seed2.dat$rec.IDs
    par.ids <- seed2.dat$don.IDs

    for(i in 2:length(seed2.dat$rec.IDs)){

      rec.id <- rec.ids[i]
      par.id <- par.ids[i]
      rec.seed <- rec.id

      if(par.id == seed.id){

        # donor coalescent tree at sampling time for the first recipient
        seed.tree.samp <- read.tree(paste("tree_at_for_",rec.seed,"_samp_net_",v,".nwk", sep = "")) # tree at sampling time of the seed
        numb.tree.samp <- numb.tr(tree=seed.tree.samp) # count number of trees generated (normally one)

        file.copy(paste("sequence_from_",seed.id,"_to_",rec.seed,"_net_",v,".fasta", sep = ""),paste("HIV.Pol.gene.bis.nwk", sep = "")) # call the seed sequences - pool of viruses
        write(numb.tree.samp,file = "HIV.Pol.gene.bis.nwk", append = TRUE) # add the number of tree in the file and
        write.tree(seed.tree.samp,file = "HIV.Pol.gene.bis.nwk", append = TRUE) # the tree, to prepare the file to simulate the evolution of the virus across the tree
        file.rename(from = "HIV.Pol.gene.bis.nwk", to = paste("HIV.Pol.gene.bis",seed.id,".nwk", sep = ""))

        seq.rand <- sample(1:10,1) # random number correponds to random sequence chosed in the pool of viruses

        # Sequence of the first recipient at the sampling (diagnosis/removal) event
        system(paste("./seq-gen -mGTR -f 0.3887782, 0.1646746, 0.2277556, 0.2187915 -a 0.9 -g 4 -r 0.00125, 0.00125, 0.00125, 0.00125, 0.00125, 0.00125 -n1 -k",seq.rand,"< HIV.Pol.gene.bis",seed.id,".nwk -z",seed,"> sequence_at_samp_",rec.seed,"_net_",v,".fasta",sep = ""))
      }
    }


    # Others

    # We know that individual j got infection from i which was transmitted by one virus from i
    # having sequences of i at transmission time, choose only one from the pool which will
    # be evolving in j



    # (iv) Sequences at infetion time >> sequence of the donors

    # Find parent of the current donor here "h": due to above mentionned argument
    parent.find <- function(h){
      for(l in 1:length(seed2.dat$rec.IDs)){
        if(seed2.dat$rec.IDs[l] == h){
          par <- seed2.dat$don.IDs[l]
        }
      }
      return(par)
    }

    for(i in 3:length(seed2.dat$rec.IDs)){ # start at 3 bcz the seed and 1st transmission already simulated
      # random number for random sequence to infect new individual
      seq.rand <- sample(1:10,1)
      # recipients - k and donors - h
      k <- seed2.dat$rec.IDs[i] # 219 recipients
      h <- seed2.dat$don.IDs[i] # 1408 donors

      prev.par <- parent.find(h)
      if(prev.par!=-1){
        # tree of the recipient
        tr.ms <- read.tree(paste("tree_at_from_",h,"_to_",k,"_net_",v,".nwk", sep = "")) # tree of donor at infecting time
        numb.trees <- numb.tr(tr.ms)
        file.copy(paste("sequence_from_",prev.par,"_to_",h,"_net_",v,".fasta", sep = ""),paste("Sequence_",h,".bis.nwk", sep = "")) # h is going to give one of what (s)he got previously from prev.par
        write(numb.trees,file = paste("Sequence_",h,".bis.nwk", sep = ""), append = TRUE)
        write.tree(tr.ms,file = paste("Sequence_",h,".bis.nwk", sep = ""), append = TRUE)
        file.rename(from = paste("Sequence_",h,".bis.nwk", sep = ""), to = paste("HIV.Pol.gene.bis",h,".nwk", sep = ""))
        # file.remove()
        system(paste("./seq-gen -mGTR -f 0.3887782, 0.1646746, 0.2277556, 0.2187915 -a 0.9 -g 4 -r 0.00125, 0.00125, 0.00125, 0.00125, 0.00125, 0.00125 -n1 -k",seq.rand,"< HIV.Pol.gene.bis",h,".nwk > sequence_from_",h,"_to_",k,"_net_",v,".fasta",sep = ""))
      }
    }

    # (v) Sequences at sampling time >> sequence of recipients

    # Find parent of the current donor here "h": due to above mentionned argument
    parent.find <- function(h){
      for(l in 1:length(seed2.dat$rec.IDs)){
        if(seed2.dat$rec.IDs[l] == h){
          par <- seed2.dat$don.IDs[l]
        }
      }
      return(par)
    }


    for(i in 3:length(seed2.dat$rec.IDs)){ # start at 3 bcz the seed and 1st transmission already simulated
      # random number for random sequence to infect new individual
      seq.rand <- sample(1:10,1)
      # recipients - k and donors - h
      k <- seed2.dat$rec.IDs[i] # [3]=219, recipients
      h <- seed2.dat$don.IDs[i] # [3]=1408, donors

      prev.par <- parent.find(h)

      if(prev.par!=-1){
        # tree of the recipient
        tr.ms <- read.tree(paste("tree_at_for_",k,"_samp_net_",v,".nwk", sep = "")) # tree of the recipient at sampling time
        numb.trees <- numb.tr(tr.ms)
        file.copy(paste("sequence_from_",prev.par,"_to_",h,"_net_",v,".fasta", sep = ""),paste("Sequence_",h,".bis.nwk", sep = "")) # pool of seq of the donor to the current donor
        write(numb.trees,file = paste("Sequence_",h,".bis.nwk", sep = ""), append = TRUE)
        write.tree(tr.ms,file = paste("Sequence_",h,".bis.nwk", sep = ""), append = TRUE)
        file.rename(from = paste("Sequence_",h,".bis.nwk", sep = ""), to = paste("HIV.Pol.gene.bis",h,".nwk", sep = ""))
        # file.remove()
        system(paste("./seq-gen -mGTR -f 0.3887782, 0.1646746, 0.2277556, 0.2187915 -a 0.9 -g 4 -r 0.00125, 0.00125, 0.00125, 0.00125, 0.00125, 0.00125 -n1 -k",seq.rand,"< HIV.Pol.gene.bis",h,".nwk > sequence_at_samp_",k,"_net_",v,".fasta",sep = ""))
      }
    }


    ## Construct then consensus sequences ##
    ########################################
    #
    # library(DECIPHER)
    # library(phyclust)

    label.names <- vector()
    for(i in 1:length(seed2.dat$rec.IDs)){
      k <- seed2.dat$rec.IDs[i]
      dna.samp <- read.phylip(file = paste("sequence_at_samp_",k,"_net_",v,".fasta", sep = ""), sep = ",")
      matrix.dna <- dna.samp$org.code
      dna.vec <- vector()
      for(i in 1:nrow(matrix.dna)){
        dna.vec <- c(dna.vec,matrix.dna[i,])
      }
      consensus.seq <- ConsensusSequence(DNAStringSet(dna.vec), threshold=0.8)
      write.dna(consensus.seq,file = paste("consensus.seq.raw_net_",v,".fas", sep = ""), format = "fasta", nbcol=-1,
                append = TRUE) # the name of seq all are "1"

      label.names <- c(label.names,k) # Individuals IDs
    }

    # read the consensus sequence all with name "1" and rename them accoridng to individuals IDs
    d <- read.FASTA(file = paste("consensus.seq.raw_net_",v,".fas", sep = ""))
    names(d) <- label.names

    write.dna(d,file = paste("consensus.seq.final_net_",v,"_model2.fas", sep = ""), format = "fasta", nbcol=-1)


  } # end loop of chosen transmission networks (with at least 3 people)

} # end loop of transmission networks

