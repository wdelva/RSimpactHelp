# ## For sub-optimal sequence coverage
#

setwd("/home/david/RSimpactHelp/R/Projects_2017/Report_1/")

## Load required packages

pacman::p_load(devtools, Rcpp, ape, expoTree, data.table, phylosim, RSimpactCyan,
               RSimpactHelper, readr, phangorn, Biostrings, dplyr, adephylo,
               phyclust, DECIPHER,treedater,geiger,picante)


## Run Simpact for specific parameter combination

age.distr <- agedistr.creator(shape = 5, scale = 65)
cfg.list <- input.params.creator(population.eyecap.fraction = 0.2, #0.21,#1,
                                 population.simtime = 15, #20, #40,  #25 for validation. 20 for calibration
                                 population.nummen = 500, # 3800, #2500,
                                 population.numwomen = 700, #4200, #2500,
                                 hivseed.time = 5, # 10,
                                 hivseed.type = "amount",
                                 hivseed.amount = 20, #30,
                                 hivseed.age.min = 20,
                                 hivseed.age.max = 50,
                                 hivtransmission.param.a = -1,
                                 hivtransmission.param.b = -90,
                                 hivtransmission.param.c = 0.5,
                                 hivtransmission.param.f1 = log(2), #log(inputvector[2]) , #log(2),
                                 hivtransmission.param.f2 = log(log(1.4) / log(2)) / 5, #log(log(sqrt(inputvector[2])) / log(inputvector[2])) / 5, #log(log(1.4) / log(2)) / 5,
                                 formation.hazard.agegapry.gap_factor_man_age = -0.01, #-0.01472653928518528523251061,
                                 formation.hazard.agegapry.gap_factor_woman_age = -0.01, #-0.0726539285185285232510561,
                                 formation.hazard.agegapry.meanage = -0.025,
                                 formation.hazard.agegapry.gap_factor_man_const = 0,
                                 formation.hazard.agegapry.gap_factor_woman_const = 0,
                                 formation.hazard.agegapry.gap_factor_man_exp = -1, #-6,#-1.5,
                                 formation.hazard.agegapry.gap_factor_woman_exp = -1, #-6,#-1.5,
                                 formation.hazard.agegapry.gap_agescale_man = 0.25, #inputvector[3], # 0.25,
                                 formation.hazard.agegapry.gap_agescale_woman = 0.25, #inputvector[3], # 0.25,#-0.30000007,#-0.03,
                                 debut.debutage = 15,
                                 conception.alpha_base = -2.5#inputvector[14]#-2.5#,
                                 #person.art.accept.threshold.dist.fixed.value = 0
)


cfg.list["formation.hazard.agegapry.baseline"] <- 2
cfg.list["mortality.aids.survtime.C"] <- 65
cfg.list["mortality.aids.survtime.k"] <- -0.2
cfg.list["monitoring.fraction.log_viralload"] <- 0.3
cfg.list["dropout.interval.dist.uniform.min"] <- 100
cfg.list["dropout.interval.dist.uniform.max"] <- 200

cfg.list["person.survtime.logoffset.dist.type"] <- "normal"
cfg.list["person.survtime.logoffset.dist.normal.mu"] <- 0
cfg.list["person.survtime.logoffset.dist.normal.sigma"] <- 0.1

cfg.list["person.agegap.man.dist.type"] <- "normal" #fixed
#cfg.list["person.agegap.man.dist.fixed.value"] <- -6
cfg.list["person.agegap.woman.dist.type"] <- "normal" #"fixed"
#cfg.list["person.agegap.woman.dist.fixed.value"] <- -6

cfg.list["mortality.aids.survtime.C"] <- 65
cfg.list["mortality.aids.survtime.k"] <- -0.2
cfg.list["monitoring.cd4.threshold"] <- 0






# Let's introduce ART, and evaluate whether the HIV prevalence drops less  rapidly
art.intro <- list()
art.intro["time"] <- 25
art.intro["person.art.accept.threshold.dist.fixed.value"] <- 0.5 # inputvector[4] ######### 0.5
art.intro["diagnosis.baseline"] <- 100 # 0#100
art.intro["monitoring.cd4.threshold"] <- 100 # 1200

### add something about diagnosis
art.intro["diagnosis.agefactor"] <- 0
art.intro["diagnosis.genderfactor"] <- 0
art.intro["diagnosis.diagpartnersfactor"] <- 0
art.intro["diagnosis.isdiagnosedfactor"] <- 0
### end of add-on about diagnosis



#art.intro["monitoring.interval.piecewise.cd4s"] <- "0,1300"


# Gradual increase in CD4 threshold. in 2007:200. in 2010:350. in 2013:500

art.intro2 <- list()
art.intro2["time"] <- 25 + 5 # inputvector[5] ######### 30
art.intro2["monitoring.cd4.threshold"] <- 200

art.intro3 <- list()
art.intro3["time"] <- 33 # inputvector[4] + inputvector[5] + inputvector[6] ########### 33
art.intro3["monitoring.cd4.threshold"] <- 350

art.intro4 <- list()
art.intro4["time"] <- 3 # inputvector[4] + inputvector[5] + inputvector[6] + inputvector[7] ########### 36
art.intro4["monitoring.cd4.threshold"] <- 500

art.intro5 <- list()
art.intro5["time"] <- 38
art.intro5["monitoring.cd4.threshold"] <- 5000 # This is equivalent to immediate access
art.intro5["person.art.accept.threshold.dist.fixed.value"] <- 0.5 # inputvector[8] ########### 0.75

# tasp.indicator <- inputvector[9] # 1 if the scenario is TasP, 0 if the scenario is current status

interventionlist <- list(art.intro, art.intro2, art.intro3, art.intro4, art.intro5)

intervention <- interventionlist # scenario(interventionlist, tasp.indicator)




inputvector <- c(123,1.1, 0.25, 0, 3, 0.23, 0.23, 45, 45, -0.5, 2.8, -0.2, -0.2, -2.5, -0.52, -0.05)


cfg.list["hivtransmission.param.f1"] = log(inputvector[2])
cfg.list["hivtransmission.param.f2"] = log(log(sqrt(inputvector[2])) / log(inputvector[2])) / 5
cfg.list["formation.hazard.agegapry.gap_agescale_man"] = inputvector[3]
cfg.list["formation.hazard.agegapry.gap_agescale_woman"] = inputvector[3]
cfg.list["person.agegap.man.dist.normal.mu"] <- inputvector[4]
cfg.list["person.agegap.woman.dist.normal.mu"] <- inputvector[4]
cfg.list["person.agegap.man.dist.normal.sigma"] <- inputvector[5]
cfg.list["person.agegap.woman.dist.normal.sigma"] <- inputvector[5]
cfg.list["person.eagerness.man.dist.gamma.a"] <- inputvector[6]
cfg.list["person.eagerness.woman.dist.gamma.a"] <- inputvector[7]
cfg.list["person.eagerness.man.dist.gamma.b"] <- inputvector[8]
cfg.list["person.eagerness.woman.dist.gamma.b"] <- inputvector[9]




cfg <- cfg.list

cfg["population.maxevents"] <- as.numeric(cfg.list["population.simtime"][1]) * as.numeric(cfg.list["population.nummen"][1]) * 3
cfg["monitoring.fraction.log_viralload"] <- 0.3
cfg["person.vsp.toacute.x"] <- 5 # See Bellan PLoS Medicine

seedid <- inputvector[1]
#cfg["person.agegap.man.dist.fixed.value"] <- -2 # inputvector[2]
#cfg["person.agegap.woman.dist.fixed.value"] <- -2 # inputvector[2]
cfg["formation.hazard.agegapry.gap_factor_man_exp"] <- inputvector[10] ######### -0.5
cfg["formation.hazard.agegapry.gap_factor_woman_exp"] <- inputvector[10] ######### -0.5
cfg["formation.hazard.agegapry.baseline"] <- inputvector[11]

cfg["formation.hazard.agegapry.numrel_man"] <- inputvector[12]
cfg["formation.hazard.agegapry.numrel_woman"] <- inputvector[13]
cfg["conception.alpha_base"] <- inputvector[14] #is conception.alpha.base (higher up)
cfg["dissolution.alpha_0"] <- inputvector[15]
cfg["dissolution.alpha_4"] <- inputvector[16]


# Run Simpact
results <- simpact.run(configParams = cfg,
                       destDir = "temp",
                       agedist = age.distr,
                       seed = seedid,
                       intervention = intervention)

datalist <- readthedata(results)


# Prevlence per age group and time at a given time-point
prevalence.df <- prevalence.calculator(datalist = datalist,
                                       agegroup = c(15, 30),
                                       timepoint = 10)

prevalence.df.plot <-prevalence.plotter(datalist = datalist, agegroup = c(15, 30))


# Incidence per age group in a given time-window
incidence.df <- incidence.calculator(datalist = datalist,
                                     agegroup = c(15, 30), timewindow = c(5, 15))


# Resource required RSimpactHelp function in my branch
source("/home/david/RSimpactHelp/R/transmNetworkBuilder.diff2.R")

simpact.trans.net <- transmNetworkBuilder.diff2(datalist = datalist, endpoint = 40)


## Prepare parameters for sequence simulations scenarios

hiv.seq.A.pol <- read_file("~/RSimpactHelp/R/Projects_2017/Report_1/HIV_1_A_single_pol.fas")
hiv.seq.C.pol <- read_file("~/RSimpactHelp/R/Projects_2017/Report_1/HIV_1_C_single_pol.fas")
hiv.seq.D.pol <- read_file("~/RSimpactHelp/R/Projects_2017/Report_1/HIV_1_D_single_pol.fas")
hiv.seq.G.pol <- read_file("~/RSimpactHelp/R/Projects_2017/Report_1/HIV_1_G_single_pol.fas")
## Remove the break line in the string of DNA
clean.hiv.seq.A.pol <-  gsub("\n", "", hiv.seq.A.pol)
clean.hiv.seq.C.pol <-  gsub("\n", "", hiv.seq.C.pol)
clean.hiv.seq.D.pol <-  gsub("\n", "", hiv.seq.D.pol)
clean.hiv.seq.G.pol <-  gsub("\n", "", hiv.seq.G.pol)


freq.clean.hiv.seq.A.pol <- letterFrequency(clean.hiv.seq.A.pol, letters = c("A", "C", "G", "T"))/nchar(clean.hiv.seq.A.pol)
freq.clean.hiv.seq.C.pol <- letterFrequency(clean.hiv.seq.C.pol, letters = c("A", "C", "G", "T"))/nchar(clean.hiv.seq.C.pol)
freq.clean.hiv.seq.D.pol <- letterFrequency(clean.hiv.seq.D.pol, letters = c("A", "C", "G", "T"))/nchar(clean.hiv.seq.D.pol)
freq.clean.hiv.seq.G.pol <- letterFrequency(clean.hiv.seq.G.pol, letters = c("A", "C", "G", "T"))/nchar(clean.hiv.seq.G.pol)


###################################################################################################
# Scenario 1
#
# one subtype of the virus (HIV-1-C) for all seeds
# complete sampling for a transmission network of one seed
# same sampling time interval (e.g.: five or three years) for a transmission network of one seed



## (i) Sequence simulation

# Use external tool seq-gen it is fast more than phylosim embeded in RSimpactHelp

# Note: transmission network with less than 3 individuals will not be considered

seed=123
trans.net <- simpact.trans.net # all transmission networks
num.trees <- vector() # ID of seeds # will be 0 because the transformation required to handle transmission tree from transmission network
# constrained to rename IDs to -1, 0, 1, 2, ...
num.i <- vector() # i_th seed in the list of seeds
for(i in 1:length(trans.net)){

  tree.n <- trans.net[[i]] # transmission network for i^th seed

  if(nrow(as.data.frame(tree.n)) >= 3){
    tree.i <- trans.network2tree(transnetwork = tree.n)
    num.trees <- c(num.trees,tree.n$id[1])
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
    seq.rand <- 3
    #
    # tr <- read.tree(file = paste("tree.model1.seed",i,".nwk", sep = ""))
    # n.tr <- numb.tr(tree = tr)

    n.tr <- 1

    seed.id <- tree.n$id[1]
    # # call the seed sequences - pool of viruses and rename the file
    file.copy(paste("HIV_1_C_single_pol.fas", sep = ""),paste("seed.seq.bis",i,".nwk", sep = ""))
    # add the number of tree in the file and
    write(n.tr,file = paste("seed.seq.bis",i,".nwk",sep = ""), append = TRUE)  # n.tr
    # the tree, to prepare the file to simulate the evolution of the virus across the tree
    write.tree(tr,file = paste("seed.seq.bis",i,".nwk", sep = ""), append = TRUE)
    file.rename(from = paste("seed.seq.bis",i,".nwk", sep = ""), to = paste("seed.seq.bis",i,"Id",seed.id,".nwk", sep = ""))

    system(paste("./seq-gen -mGTR -f 0.3887782, 0.1646746, 0.2277556, 0.2187915 -a 0.9 -g 4 -r 0.00125, 0.00125, 0.00125, 0.00125, 0.00125, 0.00125  -n1 -k",seq.rand,"< seed.seq.bis",i,"Id",seed.id,".nwk -z",seed," > sequence_all_seed_number_",i,"_model1.fasta",sep = ""))

    # a: shape parameter of Gamma > Gamma Rate Heterogeneity
    # g: category of Gamma > Discrete Gamma Rate Heterogeneity
    # r: rate matrix
    # z: seed

    # Keep sampling dates
    id.samplingtime <- as.data.frame(cbind(tree.n$id, tree.n$dtimes)) # IDs and their samling times in the transmission network

    write.csv(id.samplingtime,file=paste("samplingtimes_seed_number_",i,".csv", sep = ""))

  }
}

# Chosen transmission networks with at least 3 individuals

IDs.transm <- num.i # vector of of seeds chosen in the list of seeds



## (ii) Construct time stamped phylogenetic trees


# Function to tranform dates in named vector to be handled by treedater

dates.Transform.NamedVector  <- function(dates=dates){
  dates.val <- 1990+dates$V2 # dates
  names(dates.val) <- as.character(dates$V1) # names are the names of the tips
  return(dates.val)
}


## (ii.a) Construct phylogenetic trees

# Many tools to build the trees: within R like ape, phangorn, and outside compiled tools like iq-tree and FastTree

for (j in 1:length(IDs.transm)){

  id.trans <- IDs.transm[j]
  # External IQ-TREE
  # system(" ./iqtree-omp -s HIV.Env.gene.fasta -m GTR+R4 -nt AUTO -alrt 1000 -bb 1000")
  # Consensus tree written to HIV.Env.gene.fasta.contree
  # Reading input trees file HIV.Env.gene.fasta.contree
  # Log-likelihood of consensus tree: -10565.685
  #
  # Analysis results written to:
  #   IQ-TREE report:                HIV.Env.gene.fasta.iqtree
  # Maximum-likelihood tree:       HIV.Env.gene.fasta.treefile
  # Likelihood distances:          HIV.Env.gene.fasta.mldist
  #
  # Ultrafast bootstrap approximation results written to:
  #   Split support values:          HIV.Env.gene.fasta.splits.nex
  # Consensus tree:                HIV.Env.gene.fasta.contree
  # Screen log file:               HIV.Env.gene.fasta.log

  system(paste("./iqtree-omp -s", paste("sequence_all_seed_number_",id.trans,"_model1.fasta", sep = ""), "-m GTR+R4 -nt AUTO -alrt 1000 -bb 1000"))


}



## (ii.b) Internal node optimisation


# Named the vector of sampling times by IDs (crucil to match tips IDs and sampling date)
for (j in 1:length(IDs.transm)){

  id.trans <- IDs.transm[j]

  tree.const <- read.tree(paste("sequence_all_seed_number_",id.trans,"_model1.fasta.treefile", sep = ""))

  samp.dates <- read.csv(paste("samplingtimes_seed_number_",id.trans,".csv", sep = ""))

  time.samp <- dates.Transform.NamedVector(dates=samp.dates)

  tree.tips <- as.numeric(tree.const$tip.label)

  Ord.tree.dates <- vector()
  for(i in 1:length(tree.tips)){
    for(j in 1:length(time.samp)){
      if(tree.tips[i] == samp.dates$V1[j]){
        Ord.tree.dates <- c(Ord.tree.dates, time.samp[j])
      }
    }
  }

  # Use of library(treedater) to calibrate internal nodes
  dater.tree <- dater(tree.const, Ord.tree.dates, s = 1000) # s is the length of sequence

  # Save the tree
  write.tree(dater.tree, file = paste("calibratedTree_",id.trans,".nwk", sep = ""))
  ###########################################################################################
  ###########################################################################################
  #
  #   # With changed label names, make sure the vector of sampling dates is renamed as wells
  #
  #   # Rename the tree by adding sampling date
  # g <- tree.const
  # a <- g$tip.label # label names
  # d <- as.character(round(Ord.tree.dates, digits = 2)) # sampling dates
  # g$tip.label <- paste(a, "_", d, sep = "") # label_names+sampling_dates
  #
  # renamed.Ord.tree.dates <- Ord.tree.dates
  # names(renamed.Ord.tree.dates) <- paste(a, "_", d, sep = "") # label_names+sampling_dates
  #
  # # Use of library(treedater) to calibrate internal nodes
  # dater.tree.g <- dater(g, renamed.Ord.tree.dates, s = 1000) # s is the length of sequence

  #   # Save the tree
  #   write.tree(dater.tree.g, file = paste("calibratedTree_",id.trans,"_pref.nwk", sep = ""))

}





# Scenario 2
#
# one subtype of the virus (HIV-1-C) for all seeds
# complete sampling for a transmission network of all seeds
# same sampling time interval (e.g.: five or three years) for a transmission network of all seeds


## (i) Sequence simulation

# Use external tool seq-gen it is fast more than phylosim embeded in RSimpactHelp

# Note: transmission network with less than 3 individuals will not be considered

seed=123
trans.net <- simpact.trans.net # all transmission networks
num.trees <- vector() # ID of seeds # will be 0 because the transformation required to handle transmission tree from transmission network
# constrained to rename IDs to -1, 0, 1, 2, ...
num.i <- vector() # i_th seed in the list of seeds
for(i in 1:length(trans.net)){

  tree.n <- trans.net[[i]] # transmission network for i^th seed

  if(nrow(as.data.frame(tree.n)) >= 3){
    tree.i <- trans.network2tree(transnetwork = tree.n)
    num.trees <- c(num.trees,tree.n$id[1])
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
    seq.rand <- 3
    #
    # tr <- read.tree(file = paste("tree.model1.seed",i,".nwk", sep = ""))
    # n.tr <- numb.tr(tree = tr)

    n.tr <- 1

    seed.id <- tree.n$id[1]
    # # call the seed sequences - pool of viruses and rename the file
    file.copy(paste("HIV_1_C_single_pol.fas", sep = ""),paste("seed.seq.bis",i,".nwk", sep = ""))
    # add the number of tree in the file and
    write(n.tr,file = paste("seed.seq.bis",i,".nwk",sep = ""), append = TRUE)  # n.tr
    # the tree, to prepare the file to simulate the evolution of the virus across the tree
    write.tree(tr,file = paste("seed.seq.bis",i,".nwk", sep = ""), append = TRUE)
    file.rename(from = paste("seed.seq.bis",i,".nwk", sep = ""), to = paste("seed.seq.bis",i,"Id",seed.id,".nwk", sep = ""))

    system(paste("./seq-gen -mGTR -f 0.3887782, 0.1646746, 0.2277556, 0.2187915 -a 0.9 -g 4 -r 0.00125, 0.00125, 0.00125, 0.00125, 0.00125, 0.00125  -n1 -k",seq.rand,"< seed.seq.bis",i,"Id",seed.id,".nwk -z",seed," > sequence_all_seed_number_",i,"_model1.fasta",sep = ""))

    # a: shape parameter of Gamma > Gamma Rate Heterogeneity
    # g: category of Gamma > Discrete Gamma Rate Heterogeneity
    # r: rate matrix
    # z: seed

    # Keep sampling dates
    id.samplingtime <- as.data.frame(cbind(tree.n$id, tree.n$dtimes)) # IDs and their samling times in the transmission network

    write.csv(id.samplingtime,file=paste("samplingtimes_seed_number_",i,".csv", sep = ""))

  }
}

# Chosen transmission networks with at least 3 individuals

IDs.transm <- num.i # vector of of seeds chosen in the list of seeds



## (ii) Construct time stamped phylogenetic trees


# Function to tranform dates in named vector to be handled by treedater

dates.Transform.NamedVector  <- function(dates=dates){
  dates.val <- 1990+dates$V2 # dates
  names(dates.val) <- as.character(dates$V1) # names are the names of the tips
  return(dates.val)
}


## (ii.a) Construct phylogenetic trees

# Many tools to build the trees: within R like ape, phangorn, and outside compiled tools like iq-tree and FastTree

for (j in 1:length(IDs.transm)){

  id.trans <- IDs.transm[j]
  # External IQ-TREE
  # system(" ./iqtree-omp -s HIV.Env.gene.fasta -m GTR+R4 -nt AUTO -alrt 1000 -bb 1000")
  # Consensus tree written to HIV.Env.gene.fasta.contree
  # Reading input trees file HIV.Env.gene.fasta.contree
  # Log-likelihood of consensus tree: -10565.685
  #
  # Analysis results written to:
  #   IQ-TREE report:                HIV.Env.gene.fasta.iqtree
  # Maximum-likelihood tree:       HIV.Env.gene.fasta.treefile
  # Likelihood distances:          HIV.Env.gene.fasta.mldist
  #
  # Ultrafast bootstrap approximation results written to:
  #   Split support values:          HIV.Env.gene.fasta.splits.nex
  # Consensus tree:                HIV.Env.gene.fasta.contree
  # Screen log file:               HIV.Env.gene.fasta.log

  system(paste("./iqtree-omp -s", paste("sequence_all_seed_number_",id.trans,"_model1.fasta", sep = ""), "-m GTR+R4 -nt AUTO -alrt 1000 -bb 1000"))


}



## (ii.b) Internal node optimisation


# Named the vector of sampling times by IDs (crucil to match tips IDs and sampling date)
for (j in 1:length(IDs.transm)){

  id.trans <- IDs.transm[j]

  tree.const <- read.tree(paste("sequence_all_seed_number_",id.trans,"_model1.fasta.treefile", sep = ""))

  samp.dates <- read.csv(paste("samplingtimes_seed_number_",id.trans,".csv", sep = ""))

  time.samp <- dates.Transform.NamedVector(dates=samp.dates)

  tree.tips <- as.numeric(tree.const$tip.label)

  Ord.tree.dates <- vector()
  for(i in 1:length(tree.tips)){
    for(j in 1:length(time.samp)){
      if(tree.tips[i] == samp.dates$V1[j]){
        Ord.tree.dates <- c(Ord.tree.dates, time.samp[j])
      }
    }
  }

  # Use of library(treedater) to calibrate internal nodes
  dater.tree <- dater(tree.const, Ord.tree.dates, s = 1000) # s is the length of sequence

  # Save the tree
  write.tree(dater.tree, file = paste("calibratedTree_",id.trans,".nwk", sep = ""))
  ###########################################################################################
  ###########################################################################################
  #
  #   # With changed label names, make sure the vector of sampling dates is renamed as wells
  #
  #   # Rename the tree by adding sampling date
  # g <- tree.const
  # a <- g$tip.label # label names
  # d <- as.character(round(Ord.tree.dates, digits = 2)) # sampling dates
  # g$tip.label <- paste(a, "_", d, sep = "") # label_names+sampling_dates
  #
  # renamed.Ord.tree.dates <- Ord.tree.dates
  # names(renamed.Ord.tree.dates) <- paste(a, "_", d, sep = "") # label_names+sampling_dates
  #
  # # Use of library(treedater) to calibrate internal nodes
  # dater.tree.g <- dater(g, renamed.Ord.tree.dates, s = 1000) # s is the length of sequence

  #   # Save the tree
  #   write.tree(dater.tree.g, file = paste("calibratedTree_",id.trans,"_pref.nwk", sep = ""))

}



# Scenario 3
#
# different subtypes of the virus (HIV-1-A-C-D-G) for all seeds
# complete sampling for a transmission network of one seed
# same sampling time interval (e.g.: five or three years) for a transmission network of one seed


## (i) Sequence simulation

# Use external tool seq-gen it is fast more than phylosim embeded in RSimpactHelp

# Note: transmission network with less than 3 individuals will not be considered

seed=123
trans.net <- simpact.trans.net # all transmission networks
num.trees <- vector() # ID of seeds # will be 0 because the transformation required to handle transmission tree from transmission network
# constrained to rename IDs to -1, 0, 1, 2, ...
num.i <- vector() # i_th seed in the list of seeds
for(i in 1:length(trans.net)){

  tree.n <- trans.net[[i]] # transmission network for i^th seed

  if(nrow(as.data.frame(tree.n)) >= 3){
    tree.i <- trans.network2tree(transnetwork = tree.n)
    num.trees <- c(num.trees,tree.n$id[1])
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
    seq.rand <- 3
    #
    # tr <- read.tree(file = paste("tree.model1.seed",i,".nwk", sep = ""))
    # n.tr <- numb.tr(tree = tr)

    n.tr <- 1

    seed.id <- tree.n$id[1]
    # # call the seed sequences - pool of viruses and rename the file
    file.copy(paste("HIV_1_A_C_D_G_pol.fas", sep = ""),paste("seed.seq.bis",i,".nwk", sep = ""))
    # add the number of tree in the file and
    write(n.tr,file = paste("seed.seq.bis",i,".nwk",sep = ""), append = TRUE)  # n.tr
    # the tree, to prepare the file to simulate the evolution of the virus across the tree
    write.tree(tr,file = paste("seed.seq.bis",i,".nwk", sep = ""), append = TRUE)
    file.rename(from = paste("seed.seq.bis",i,".nwk", sep = ""), to = paste("seed.seq.bis",i,"Id",seed.id,".nwk", sep = ""))

    system(paste("./seq-gen -mGTR -f 0.3887782, 0.1646746, 0.2277556, 0.2187915 -a 0.9 -g 4 -r 0.00125, 0.00125, 0.00125, 0.00125, 0.00125, 0.00125  -n1 -k",seq.rand,"< seed.seq.bis",i,"Id",seed.id,".nwk -z",seed," > sequence_all_seed_number_",i,"_model1.fasta",sep = ""))

    # a: shape parameter of Gamma > Gamma Rate Heterogeneity
    # g: category of Gamma > Discrete Gamma Rate Heterogeneity
    # r: rate matrix
    # z: seed

    # Keep sampling dates
    id.samplingtime <- as.data.frame(cbind(tree.n$id, tree.n$dtimes)) # IDs and their samling times in the transmission network

    write.csv(id.samplingtime,file=paste("samplingtimes_seed_number_",i,".csv", sep = ""))

  }
}

# Chosen transmission networks with at least 3 individuals

IDs.transm <- num.i # vector of of seeds chosen in the list of seeds



## (ii) Construct time stamped phylogenetic trees


# Function to tranform dates in named vector to be handled by treedater

dates.Transform.NamedVector  <- function(dates=dates){
  dates.val <- 1990+dates$V2 # dates
  names(dates.val) <- as.character(dates$V1) # names are the names of the tips
  return(dates.val)
}


## (ii.a) Construct phylogenetic trees

# Many tools to build the trees: within R like ape, phangorn, and outside compiled tools like iq-tree and FastTree

for (j in 1:length(IDs.transm)){

  id.trans <- IDs.transm[j]
  # External IQ-TREE
  # system(" ./iqtree-omp -s HIV.Env.gene.fasta -m GTR+R4 -nt AUTO -alrt 1000 -bb 1000")
  # Consensus tree written to HIV.Env.gene.fasta.contree
  # Reading input trees file HIV.Env.gene.fasta.contree
  # Log-likelihood of consensus tree: -10565.685
  #
  # Analysis results written to:
  #   IQ-TREE report:                HIV.Env.gene.fasta.iqtree
  # Maximum-likelihood tree:       HIV.Env.gene.fasta.treefile
  # Likelihood distances:          HIV.Env.gene.fasta.mldist
  #
  # Ultrafast bootstrap approximation results written to:
  #   Split support values:          HIV.Env.gene.fasta.splits.nex
  # Consensus tree:                HIV.Env.gene.fasta.contree
  # Screen log file:               HIV.Env.gene.fasta.log

  system(paste("./iqtree-omp -s", paste("sequence_all_seed_number_",id.trans,"_model1.fasta", sep = ""), "-m GTR+R4 -nt AUTO -alrt 1000 -bb 1000"))


}



## (ii.b) Internal node optimisation


# Named the vector of sampling times by IDs (crucil to match tips IDs and sampling date)
for (j in 1:length(IDs.transm)){

  id.trans <- IDs.transm[j]

  tree.const <- read.tree(paste("sequence_all_seed_number_",id.trans,"_model1.fasta.treefile", sep = ""))

  samp.dates <- read.csv(paste("samplingtimes_seed_number_",id.trans,".csv", sep = ""))

  time.samp <- dates.Transform.NamedVector(dates=samp.dates)

  tree.tips <- as.numeric(tree.const$tip.label)

  Ord.tree.dates <- vector()
  for(i in 1:length(tree.tips)){
    for(j in 1:length(time.samp)){
      if(tree.tips[i] == samp.dates$V1[j]){
        Ord.tree.dates <- c(Ord.tree.dates, time.samp[j])
      }
    }
  }

  # Use of library(treedater) to calibrate internal nodes
  dater.tree <- dater(tree.const, Ord.tree.dates, s = 1000) # s is the length of sequence

  # Save the tree
  write.tree(dater.tree, file = paste("calibratedTree_",id.trans,".nwk", sep = ""))
  ###########################################################################################
  ###########################################################################################
  #
  #   # With changed label names, make sure the vector of sampling dates is renamed as wells
  #
  #   # Rename the tree by adding sampling date
  # g <- tree.const
  # a <- g$tip.label # label names
  # d <- as.character(round(Ord.tree.dates, digits = 2)) # sampling dates
  # g$tip.label <- paste(a, "_", d, sep = "") # label_names+sampling_dates
  #
  # renamed.Ord.tree.dates <- Ord.tree.dates
  # names(renamed.Ord.tree.dates) <- paste(a, "_", d, sep = "") # label_names+sampling_dates
  #
  # # Use of library(treedater) to calibrate internal nodes
  # dater.tree.g <- dater(g, renamed.Ord.tree.dates, s = 1000) # s is the length of sequence

  #   # Save the tree
  #   write.tree(dater.tree.g, file = paste("calibratedTree_",id.trans,"_pref.nwk", sep = ""))

}




# Scenario 4
#
# different subtypes of the virus (HIV-1-A-C-D-G) for all seeds
# complete sampling for a transmission network of all seeds
# same sampling time interval (e.g.: five or three years) for a transmission network of all seeds


## (i) Sequence simulation

# Use external tool seq-gen it is fast more than phylosim embeded in RSimpactHelp

# Note: transmission network with less than 3 individuals will not be considered

seed=123
trans.net <- simpact.trans.net # all transmission networks
num.trees <- vector() # ID of seeds # will be 0 because the transformation required to handle transmission tree from transmission network
# constrained to rename IDs to -1, 0, 1, 2, ...
num.i <- vector() # i_th seed in the list of seeds
for(i in 1:length(trans.net)){

  tree.n <- trans.net[[i]] # transmission network for i^th seed

  if(nrow(as.data.frame(tree.n)) >= 3){
    tree.i <- trans.network2tree(transnetwork = tree.n)
    num.trees <- c(num.trees,tree.n$id[1])
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
    seq.rand <- 3
    #
    # tr <- read.tree(file = paste("tree.model1.seed",i,".nwk", sep = ""))
    # n.tr <- numb.tr(tree = tr)

    n.tr <- 1

    seed.id <- tree.n$id[1]
    # # call the seed sequences - pool of viruses and rename the file
    file.copy(paste("HIV_1_A_C_D_G_pol.fas", sep = ""),paste("seed.seq.bis",i,".nwk", sep = ""))
    # add the number of tree in the file and
    write(n.tr,file = paste("seed.seq.bis",i,".nwk",sep = ""), append = TRUE)  # n.tr
    # the tree, to prepare the file to simulate the evolution of the virus across the tree
    write.tree(tr,file = paste("seed.seq.bis",i,".nwk", sep = ""), append = TRUE)
    file.rename(from = paste("seed.seq.bis",i,".nwk", sep = ""), to = paste("seed.seq.bis",i,"Id",seed.id,".nwk", sep = ""))

    system(paste("./seq-gen -mGTR -f 0.3887782, 0.1646746, 0.2277556, 0.2187915 -a 0.9 -g 4 -r 0.00125, 0.00125, 0.00125, 0.00125, 0.00125, 0.00125  -n1 -k",seq.rand,"< seed.seq.bis",i,"Id",seed.id,".nwk -z",seed," > sequence_all_seed_number_",i,"_model1.fasta",sep = ""))

    # a: shape parameter of Gamma > Gamma Rate Heterogeneity
    # g: category of Gamma > Discrete Gamma Rate Heterogeneity
    # r: rate matrix
    # z: seed

    # Keep sampling dates
    id.samplingtime <- as.data.frame(cbind(tree.n$id, tree.n$dtimes)) # IDs and their samling times in the transmission network

    write.csv(id.samplingtime,file=paste("samplingtimes_seed_number_",i,".csv", sep = ""))

  }
}

# Chosen transmission networks with at least 3 individuals

IDs.transm <- num.i # vector of of seeds chosen in the list of seeds



## (ii) Construct time stamped phylogenetic trees


# Function to tranform dates in named vector to be handled by treedater

dates.Transform.NamedVector  <- function(dates=dates){
  dates.val <- 1990+dates$V2 # dates
  names(dates.val) <- as.character(dates$V1) # names are the names of the tips
  return(dates.val)
}


## (ii.a) Construct phylogenetic trees

# Many tools to build the trees: within R like ape, phangorn, and outside compiled tools like iq-tree and FastTree

for (j in 1:length(IDs.transm)){

  id.trans <- IDs.transm[j]
  # External IQ-TREE
  # system(" ./iqtree-omp -s HIV.Env.gene.fasta -m GTR+R4 -nt AUTO -alrt 1000 -bb 1000")
  # Consensus tree written to HIV.Env.gene.fasta.contree
  # Reading input trees file HIV.Env.gene.fasta.contree
  # Log-likelihood of consensus tree: -10565.685
  #
  # Analysis results written to:
  #   IQ-TREE report:                HIV.Env.gene.fasta.iqtree
  # Maximum-likelihood tree:       HIV.Env.gene.fasta.treefile
  # Likelihood distances:          HIV.Env.gene.fasta.mldist
  #
  # Ultrafast bootstrap approximation results written to:
  #   Split support values:          HIV.Env.gene.fasta.splits.nex
  # Consensus tree:                HIV.Env.gene.fasta.contree
  # Screen log file:               HIV.Env.gene.fasta.log

  system(paste("./iqtree-omp -s", paste("sequence_all_seed_number_",id.trans,"_model1.fasta", sep = ""), "-m GTR+R4 -nt AUTO -alrt 1000 -bb 1000"))


}



## (ii.b) Internal node optimisation


# Named the vector of sampling times by IDs (crucil to match tips IDs and sampling date)
for (j in 1:length(IDs.transm)){

  id.trans <- IDs.transm[j]

  tree.const <- read.tree(paste("sequence_all_seed_number_",id.trans,"_model1.fasta.treefile", sep = ""))

  samp.dates <- read.csv(paste("samplingtimes_seed_number_",id.trans,".csv", sep = ""))

  time.samp <- dates.Transform.NamedVector(dates=samp.dates)

  tree.tips <- as.numeric(tree.const$tip.label)

  Ord.tree.dates <- vector()
  for(i in 1:length(tree.tips)){
    for(j in 1:length(time.samp)){
      if(tree.tips[i] == samp.dates$V1[j]){
        Ord.tree.dates <- c(Ord.tree.dates, time.samp[j])
      }
    }
  }

  # Use of library(treedater) to calibrate internal nodes
  dater.tree <- dater(tree.const, Ord.tree.dates, s = 1000) # s is the length of sequence

  # Save the tree
  write.tree(dater.tree, file = paste("calibratedTree_",id.trans,".nwk", sep = ""))
  ###########################################################################################
  ###########################################################################################
  #
  #   # With changed label names, make sure the vector of sampling dates is renamed as wells
  #
  #   # Rename the tree by adding sampling date
  # g <- tree.const
  # a <- g$tip.label # label names
  # d <- as.character(round(Ord.tree.dates, digits = 2)) # sampling dates
  # g$tip.label <- paste(a, "_", d, sep = "") # label_names+sampling_dates
  #
  # renamed.Ord.tree.dates <- Ord.tree.dates
  # names(renamed.Ord.tree.dates) <- paste(a, "_", d, sep = "") # label_names+sampling_dates
  #
  # # Use of library(treedater) to calibrate internal nodes
  # dater.tree.g <- dater(g, renamed.Ord.tree.dates, s = 1000) # s is the length of sequence

  #   # Save the tree
  #   write.tree(dater.tree.g, file = paste("calibratedTree_",id.trans,"_pref.nwk", sep = ""))

}
