
setwd("/home/david/RSimpactHelp/R/Projects_2017/Example3PaperSimpact2/")

# ## For sub-optimal sequence coverage


## Load required packages

pacman::p_load(devtools, Rcpp, ape, expoTree, data.table, phylosim, RSimpactCyan,
               RSimpactHelper, readr, phangorn, Biostrings, dplyr, adephylo,igraph,
               phyclust, DECIPHER,treedater,geiger,picante)

#######################
# Step 1: Run Simpact #
#######################

## Run Simpact for specific parameter combination

age.distr <- agedistr.creator(shape = 5, scale = 65)
cfg.list <- input.params.creator(population.eyecap.fraction = 0.2, #0.21,#1,
                                 population.msm = "no",
                                 population.simtime = 40, #20, #40,  #25 for validation. 20 for calibration
                                 population.nummen = 600, # 3800, #2500,
                                 population.numwomen = 600, #4200, #2500,
                                 hivseed.time = 10, # 10,
                                 hivseed.type = "amount",
                                 hivseed.amount = 20, #30,
                                 hivseed.age.min = 20,
                                 hivseed.age.max = 50,
                                 hivtransmission.param.a = -1, # -1,
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
cfg.list["monitoring.fraction.log_viralload"] <- 0 #0.3
cfg.list["dropout.interval.dist.uniform.min"] <- 1000
cfg.list["dropout.interval.dist.uniform.max"] <- 2000

cfg.list["person.survtime.logoffset.dist.type"] <- "normal"
cfg.list["person.survtime.logoffset.dist.normal.mu"] <- 0
cfg.list["person.survtime.logoffset.dist.normal.sigma"] <- 0.1

cfg.list["person.agegap.man.dist.type"] <- "normal" #fixed
#cfg.list["person.agegap.man.dist.fixed.value"] <- -6
cfg.list["person.agegap.woman.dist.type"] <- "normal" #"fixed"
#cfg.list["person.agegap.woman.dist.fixed.value"] <- -6

cfg.list["mortality.aids.survtime.C"] <- 65
cfg.list["mortality.aids.survtime.k"] <- -0.2
cfg.list["monitoring.cd4.threshold"] <- 0 # 0 means nobody qualifies for ART
cfg.list["person.art.accept.threshold.dist.fixed.value"] <- 0.4
cfg.list["diagnosis.baseline"] <- -2






# Let's introduce ART, and evaluate whether the HIV prevalence drops less  rapidly
art.intro <- list()
art.intro["time"] <- 20
art.intro["diagnosis.baseline"] <- -2 # 0#100
art.intro["monitoring.cd4.threshold"] <- 100 # 1200

### add something about diagnosis
art.intro["diagnosis.agefactor"] <- 0
art.intro["diagnosis.genderfactor"] <- 0
art.intro["diagnosis.diagpartnersfactor"] <- 0
art.intro["diagnosis.isdiagnosedfactor"] <- 0
### end of add-on about diagnosis



#art.intro["monitoring.interval.piecewise.cd4s"] <- "0,1300"


# Gradual increase in CD4 threshold. in 2007:200. in 2010:350. in 2013:500
art.intro1 <- list()
art.intro1["time"] <- 22
art.intro1["diagnosis.baseline"] <- -2 # 0#100
art.intro1["monitoring.cd4.threshold"] <- 150 # 1200


art.intro2 <- list()
art.intro2["time"] <- 25 # inputvector[5] ######### 30
art.intro2["monitoring.cd4.threshold"] <- 200

art.intro3 <- list()
art.intro3["time"] <- 30 # inputvector[4] + inputvector[5] + inputvector[6] ########### 33
art.intro3["monitoring.cd4.threshold"] <- 350

art.intro4 <- list()
art.intro4["time"] <- 33 # inputvector[4] + inputvector[5] + inputvector[6] + inputvector[7] ########### 36
art.intro4["monitoring.cd4.threshold"] <- 500

art.intro5 <- list()
art.intro5["time"] <- 36
art.intro5["monitoring.cd4.threshold"] <- 700 # This is equivalent to immediate access

# tasp.indicator <- inputvector[9] # 1 if the scenario is TasP, 0 if the scenario is current status

interventionlist <- list(art.intro, art.intro1, art.intro2, art.intro3, art.intro4, art.intro5)

intervention <- interventionlist # scenario(interventionlist, tasp.indicator)




inputvector <- c(123, 1.05, 0.25, 0, 3, 0.23, 0.23, 45, 45, -0.7, 2.8,
                 -0.3, -0.3,
                 -2.7, # conception
                 -0.52, -0.05)


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
# cfg["monitoring.fraction.log_viralload"] <- 0.3
cfg["person.vsp.toacute.x"] <- 5 # See Bellan PLoS Medicine

seedid <- inputvector[1] + j
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


#
# # # # Run Simpact
results <- simpact.run(configParams = cfg,
                       destDir = "temp",
                       agedist = age.distr,
                       seed = seedid,
                       intervention = intervention)

datalist <- readthedata(results)
# table(datalist$etable$eventname)
#
# ## Save the output
# save(datalist, file = "MasterModelSubOptimalSeqCovearge.datalist.RData")

# Read saved output data set
# datalist <- get(load("MasterModelSubOptimalSeqCovearge.datalist.RData"))

# table(datalist$etable$eventname) # check events

# Prevlence per age group and time at a given time-point
# prevalence.df <- prevalence.calculator(datalist = datalist,
#                                        agegroup = c(15, 30),
#                                        timepoint = 10)

# prevalence.df.plot <-prevalence.plotter(datalist = datalist, agegroup = c(15, 50))


# Incidence per age group in a given time-window
# incidence.df <- incidence.calculator(datalist = datalist,
#                                      agegroup = c(15, 30), timewindow = c(10, 20))



###########################################
# Step 2: Construct transmission networks #
###########################################


# Resource required RSimpactHelp function in my branch
source("/user/data/gent/vsc400/vsc40070/phylo/transmNetworkBuilder.diff2.R")
source("/home/david/RSimpactHelp/R/transmNetworkBuilder.diff3.R")
source("/home/david/RSimpactHelp/R/trans.network2tree.R")
source("/home/david/RSimpactHelp/R/epi2tree2.R")


simpact.trans.net <- transmNetworkBuilder.diff2(datalist = datalist, endpoint = 40)

smallest.branches <- rep(NA, times = length(simpact.trans.net))
for (list.element in 1:length(simpact.trans.net)){
  net.list <- simpact.trans.net[[list.element]]
  if(length(net.list$id) > 2){
    tree.tr <- epi2tree2(net.list)
    smallest.branch <- min(tree.tr$edge.length)
    smallest.branches[list.element] <- smallest.branch
  }
}
#min(smallest.branches, na.rm = TRUE) #
## seeds and transmission network sizes:
# 2>>10,  3>>212,  8>>20, 10>>10, 12>>3, 16>>13, 19>>1741
#which(smallest.branches!="NA")

#epi.dt <- simpact.trans.net[[3]] #
#tree.tr <- trans.network2tree(epi.dt)
#summary(tree.tr) # some branch lengths are negative with Gabriel function
#plot(tree.tr, root.edge = TRUE, show.tip.label=FALSE)
#axisPhylo(backward = FALSE)


###############################
# Step 3: Sequence simulation #
###############################


# Use external tool seq-gen it is fast more than phylosim embeded in RSimpactHelp

# Note: transmission network with less than 3 individuals will not be considered

seed=123
trans.net <- simpact.trans.net # all transmission networks
num.trees <- vector() # ID of seeds # will be 0 because the transformation required to handle transmission tree from transmission network
# constrained to rename IDs to -1, 0, 1, 2, ...
num.i <- vector() # i_th seed in the list of seeds
for(i in 1:length(trans.net)){

  tree.n <- trans.net[[i]] # transmission network for i^th seed

  if(nrow(as.data.frame(tree.n)) >= 50){
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
    seq.rand <- 1
    #
    # tr <- read.tree(file = paste("tree.model1.seed",i,".nwk", sep = ""))
    # n.tr <- numb.tr(tree = tr)

    n.tr <- 1

    seed.id <- tree.n$id[1]
    # # call the seed sequences - pool of viruses and rename the file
    file.copy(paste("hiv.seq.B.pol.j.fasta", sep = ""),paste("seed.seq.bis",i,".nwk", sep = ""))
    # add the number of tree in the file and
    write(n.tr,file = paste("seed.seq.bis",i,".nwk",sep = ""), append = TRUE)  # n.tr
    # the tree, to prepare the file to simulate the evolution of the virus across the tree
    write.tree(tr,file = paste("seed.seq.bis",i,".nwk", sep = ""), append = TRUE)
    file.rename(from = paste("seed.seq.bis",i,".nwk", sep = ""), to = paste("seed.seq.bis",i,"Id",seed.id,".nwk", sep = ""))

    system(paste("seq-gen -mGTR -f 0.3857, 0.1609, 0.2234, 0.2300  -a 0.9 -g 4 -i 0.5230 -r 2.9114, 12.5112, 1.2569, 0.8559, 12.9379, 1.0000 -s 0.00475  -n1 -k",seq.rand,"< seed.seq.bis",i,"Id",seed.id,".nwk -z",seed," > B.EpidemicSequences_seed_",i,".fasta",sep = ""))

    # a: shape parameter of Gamma > Gamma Rate Heterogeneity
    # g: category of Gamma > Discrete Gamma Rate Heterogeneity
    # r: rate matrix
    # s: scale which is the substitution rate of pol gene
    # z: seed

    # # Keep sampling dates
    # id.samplingtime <- as.data.frame(cbind(tree.n$id, tree.n$dtimes)) # IDs and their samling times in the transmission network

    # Keep sampling dates, add "A" on ID-tip- labels to handle them as charcters
    id.samplingtime <- as.data.frame(cbind(paste(i,".",tree.n$id, ".A", sep = ""), tree.n$dtimes)) # IDs and their samling times in the transmission network


    write.csv(id.samplingtime,file=paste("samplingtimes_seed_",i,".csv", sep = ""))

  }
}

# Chosen transmission networks with at least 3 individuals

IDs.transm <- num.i # vector of of seeds chosen in the list of seeds


#####################################################
# Step 4: Construct time stamped phylogenetic trees #
#####################################################


# Function to tranform dates in named vector to be handled by treedater

dates.Transform.NamedVector  <- function(dates=dates){
  dates.val <- datalist$itable$population.simtime[1] - dates$V2 + 1977 # dates
  names(dates.val) <- as.character(dates$V1) # names are the names of the tips
  return(dates.val)
}

# 4.1. Many tools to build the trees: within R like ape, phangorn, and outside compiled tools like iq-tree and FastTree

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

  # system(paste("./iqtree-omp -s", paste("B.EpidemicSequences_seed_",id.trans,".fasta", sep = ""), " -nt AUTO -alrt 1000 -bb 1000"))

  # Compiling FastTree
  # gcc -DUSE_DOUBLE -O3 -finline-functions -funroll-loops -Wall -o FastTree FastTree.c -lm

  system(paste("FastTree -gtr -nt <", paste("B.EpidemicSequences_seed_",id.trans,".fasta", sep = ""), paste(">B.EpidemicSequences_seed_",id.trans,".fasta.tree", sep = "")))


}


# 4.2. Internal node optimisation requires sampled dates


for (j in 1:length(IDs.transm)){

  id.trans <- IDs.transm[j]

  tree.const <- read.tree(paste("B.EpidemicSequences_seed_",id.trans,".fasta.tree", sep = ""))

  samp.dates <- read.csv(paste("samplingtimes_seed_",id.trans,".csv", sep = ""))

  time.samp <- dates.Transform.NamedVector(dates=samp.dates)

  tree.tips <- as.character(tree.const$tip.label)

  Ord.tree.dates <- vector()
  for(i in 1:length(tree.tips)){
    for(j in 1:length(time.samp)){
      if(tree.tips[i] == samp.dates$V1[j]){
        Ord.tree.dates <- c(Ord.tree.dates, time.samp[j])
      }
    }
  }

  # Use of library(treedater) to calibrate internal nodes
  dater.tree <- dater(tree.const, Ord.tree.dates, s = 3012) # s is the length of sequence

  save(dater.tree, file = paste("dated.tree.object_seed_",id.trans,".Rdata", sep = ""))

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

# IDs.transm <- c(2,3,8,10,12,16,19)

#test.t.dated.2 <- get(load("dated.tree.object_seed_2.Rdata"))
#N.2 <- node.age(test.t.dated.2)
#int.node.age.2 <- N.2$Ti # internal nodes ages
#latest.samp.2 <- N.2$timeToMRCA+N.2$timeOfMRCA # latest sampling date
#
#plot(test.t.dated.2, label.offset=0.1)  #the label.offset argument moves the species names a bit to the right (try different values!)
#nodelabels() #add node numbers
#tiplabels() #add tip numbers



#test.t.dated.3 <- get(load("dated.tree.object_seed_3.Rdata"))
#N.3 <- node.age(test.t.dated.3)
#int.node.age.3 <- N.3$Ti # internal nodes ages
#latest.samp.3 <- N.3$timeToMRCA+N.3$timeOfMRCA # latest sampling date

#test.t.dated.8 <- get(load("dated.tree.object_seed_8.Rdata"))
#N.8 <- node.age(test.t.dated.8)
#int.node.age.8 <- N.8$Ti # internal nodes ages
#latest.samp.8 <- N.8$timeToMRCA+N.8$timeOfMRCA # latest sampling date

#test.t.dated.10 <- get(load("dated.tree.object_seed_10.Rdata"))
#N.10 <- node.age(test.t.dated.10)
#int.node.age.10 <- N.10$Ti # internal nodes ages
#latest.samp.10 <- N.10$timeToMRCA+N.10$timeOfMRCA # latest sampling date

#test.t.dated.12 <- get(load("dated.tree.object_seed_12.Rdata"))
#N.12 <- node.age(test.t.dated.12)
#int.node.age.12 <- N.12$Ti # internal nodes ages
#latest.samp.12 <- N.12$timeToMRCA+N.12$timeOfMRCA # latest sampling date

#test.t.dated.16 <- get(load("dated.tree.object_seed_16.Rdata"))
#N.16 <- node.age(test.t.dated.16)
#int.node.age.16 <- N.16$Ti # internal nodes ages
#latest.samp.16 <- N.16$timeToMRCA+N.16$timeOfMRCA # latest sampling date

#test.t.dated.19 <- get(load("dated.tree.object_seed_19.Rdata"))
#N.19 <- node.age(test.t.dated.19)
#int.node.age.19 <- N.19$Ti # internal nodes ages
#latest.samp.19 <- N.19$timeToMRCA+N.19$timeOfMRCA # latest sampling date


#### We choose seed 5 for the example

id.trans <- 5

tree.const <- read.tree(paste("B.EpidemicSequences_seed_",id.trans,".fasta.tree", sep = ""))

samp.dates <- read.csv(paste("samplingtimes_seed_",id.trans,".csv", sep = ""))

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
dater.tree <- dater(tree.const, Ord.tree.dates, s = 3012) # s is the length of sequence

# Save the tree
write.tree(dater.tree, file = paste("calibratedTree_",id.trans,".nwk", sep = ""))

# Node age with picante package

N <- node.age(dater.tree)

# potential transmission times
int.node.age <- N$Ti # internal nodes ages
##########################################

# begin 1.12.2017
dated.tree.5 <- get(load("dated.tree.object_seed_5.Rdata"))
pb <- parboot.treedater(dated.tree.5)
# end 1.12.2017

latest.samp <- N$timeToMRCA+N$timeOfMRCA # latest sampling date


## Transmission network of seed 5
trans.net <- simpact.trans.net
#int.node.age <- int.node.age.5
tra.net.5 <- trans.net[[5]]

tra.net.5$dtimes <- abs(tra.net.5$dtimes-40)+1977 #(endpoint=40)
tra.net.5$itimes <- abs(tra.net.5$itimes-40)+1977 #(endpoint=40) 10>1990, -> +1980

min.val = 1977
max.val = round(max(tra.net.5$itimes))


step.int=1
d <- (max.val-min.val)/step.int

dat.f.trans <- as.data.frame(tra.net.5)
dt.node.age.dt <- int.node.age

numb.tra <- vector()
i.vec <- vector()
int.node.vec <- vector()
for (i in 1:d) {
  inf <- 1976+i
  sup <- 1977+i
  dat.f.trans.i <- dat.f.trans[which(dat.f.trans$itimes <= sup & dat.f.trans$itimes  >= inf),]
  numb.i <- nrow(dat.f.trans.i)
  numb.tra <- c(numb.tra, numb.i)
  i.vec <- c(i.vec, sup)
  int.node.age.i <- int.node.age[int.node.age <= sup & dt.node.age.dt >= inf]
  int.node.vec <- c(int.node.vec,length(int.node.age.i))
}



graph.build <- as.data.frame(trans.net[[5]])

graph.build[,4] <- as.character(graph.build$parent) # donors
graph.build[,3] <- as.character(graph.build$id) # recipients
gag = as.matrix(graph.build)
ga.graph = graph.edgelist(gag[,4:3])

V(ga.graph)$color <- "red"

transNet.yrs.Old <- delete.vertices(ga.graph, "-1")



### Plot figures
#################

# 1. Transmission network from simpact                           # 1 #
#plot.igraph(transNet.yrs.Old, edge.arrow.size=0.1, vertex.size=5,
#            edge.color="black",
#            asp = 1,
#            xlim = c(-1, 2),
#            ylim = c(-0.5,0.5),
#            vertex.frame.color="black",
#            vertex.label.color="black",
#            #vertex.label = NULL,
#            layout = layout_with_kk,
#            edge.width = 1,
#            vertex.label.cex=0.1,
#            vertex.label.dist=0.0, edge.curved=0.0, asp=0, margin=0
#            #main = "True transmission network"
#)

# 2. Phylogenetic tree

#plot(dater.tree, show.tip.label=FALSE,
#     edge.width=1,
#     edge.color="blue") # Try a few different settings!
#axisPhylo(backward = FALSE)


# 3. Transmission events and internal nodes
#x <- i.vec
#plot(x, int.node.vec, type="b", col="red", lwd=2,
#     xlab = "Calendar time",
#     ylab = "Count") # 1 > 1
#lines(x, numb.tra, col='green3', type='b', lwd=2)
#legend("topleft", legend = c("Internal nodes", "Transmission events"),
#       col=c("red","green3"), pch=1)


# Object for plotting by ggplot
SimpactPaperPhyloExample <- list()
SimpactPaperPhyloExample$transNet.yrs.Ord <- transNet.yrs.Old
SimpactPaperPhyloExample$dater.tree <- dater.tree
SimpactPaperPhyloExample$i.vec <- i.vec
SimpactPaperPhyloExample$int.node.vec <- int.node.vec
SimpactPaperPhyloExample$numb.tra <- numb.tra
save(SimpactPaperPhyloExample, file = "SimpactPaperPhyloExample.RData")

#fig.obj <- get(load("SimpactPaperPhyloExample.RData"))
