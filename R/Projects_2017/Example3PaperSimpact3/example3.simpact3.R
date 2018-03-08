
setwd("/home/david/RSimpactHelp/R/Projects_2017/Example3PaperSimpact3/")

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
                                 population.nummen = 1500, # 3800, #2500,
                                 population.numwomen = 2000, #4200, #2500,
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
cfg.list["monitoring.fraction.log_viralload"] <- 0 # 0.3
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
cfg.list["monitoring.cd4.threshold"] <- 10000
cfg.list["person.art.accept.threshold.dist.fixed.value"] <- 1
cfg.list["diagnosis.baseline"] <- -2






# Let's introduce ART, and evaluate whether the HIV prevalence drops less  rapidly
art.intro <- list()
art.intro["time"] <- 0.0001 # 25
art.intro["person.art.accept.threshold.dist.fixed.value"] <- 1 # 0.5 # inputvector[4] ######### 0.5
art.intro["diagnosis.baseline"] <- -2 # 0#100
art.intro["monitoring.cd4.threshold"] <- 10000 # 1200

### add something about diagnosis
art.intro["diagnosis.agefactor"] <- 0
art.intro["diagnosis.genderfactor"] <- 0
art.intro["diagnosis.diagpartnersfactor"] <- 0
art.intro["diagnosis.isdiagnosedfactor"] <- 0
### end of add-on about diagnosis



#art.intro["monitoring.interval.piecewise.cd4s"] <- "0,1300"


# Gradual increase in CD4 threshold. in 2007:200. in 2010:350. in 2013:500

# art.intro2 <- list()
# art.intro2["time"] <- 25 + 5 # inputvector[5] ######### 30
# art.intro2["monitoring.cd4.threshold"] <- 200
#
# art.intro3 <- list()
# art.intro3["time"] <- 33 # inputvector[4] + inputvector[5] + inputvector[6] ########### 33
# art.intro3["monitoring.cd4.threshold"] <- 350
#
# art.intro4 <- list()
# art.intro4["time"] <- 3 # inputvector[4] + inputvector[5] + inputvector[6] + inputvector[7] ########### 36
# art.intro4["monitoring.cd4.threshold"] <- 500
#
# art.intro5 <- list()
# art.intro5["time"] <- 38
# art.intro5["monitoring.cd4.threshold"] <- 5000 # This is equivalent to immediate access
# art.intro5["person.art.accept.threshold.dist.fixed.value"] <- 0.5 # inputvector[8] ########### 0.75

# tasp.indicator <- inputvector[9] # 1 if the scenario is TasP, 0 if the scenario is current status

interventionlist <- list(art.intro) #, art.intro2, art.intro3, art.intro4, art.intro5)

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
# cfg["monitoring.fraction.log_viralload"] <- 0.3
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

# # #
# # # # # Run Simpact
results <- simpact.run(configParams = cfg,
                       destDir = "temp",
                       agedist = age.distr,
                       seed = seedid,
                       intervention = intervention)
#
datalist <- readthedata(results)
# #
# # table(datalist$etable$eventname)
# #
# # ## Save the output
# save(datalist, file = "MasterModelSubOptimalSeqCovearge.datalist.RData")

# Read saved output data set
# datalist <- get(load("MasterModelSubOptimalSeqCovearge.datalist.RData"))

# table(datalist$etable$eventname) # check events

# # Prevlence per age group and time at a given time-point
# prevalence.df <- prevalence.calculator(datalist = datalist,
#                                        agegroup = c(15, 30),
#                                        timepoint = 10)
#
# prevalence.df.plot <-prevalence.plotter(datalist = datalist, agegroup = c(15, 50))
#
#
# # Incidence per age group in a given time-window
# incidence.df <- incidence.calculator(datalist = datalist,
#                                      agegroup = c(15, 30), timewindow = c(10, 20))



###########################################
# Step 2: Construct transmission networks #
###########################################


# Resource required RSimpactHelp function in my branch
source("/home/david/RSimpactHelp/R/transmNetworkBuilder.diff3.R")
source("/home/david/RSimpactHelp/R/trans.network2tree.R")
source("/home/david/RSimpactHelp/R/epi2tree2.R")


simpact.trans.net <- transmNetworkBuilder.diff3(datalist = datalist, endpoint = 40)


# save(simpact.trans.net, file = "simpact.trans.net.RData")

# Read saved output data set
# simpact.trans.net <- get(load("simpact.trans.net.RData"))

# Check trees with negative branch lengths and at least with 3 leaves

# smallest.branches <- rep(NA, times = length(simpact.trans.net))
# for (list.element in 1:length(simpact.trans.net)){
#   net.list <- simpact.trans.net[[list.element]]
#   if(length(net.list$id) > 2){
#     tree.tr <- epi2tree2(net.list)
#     smallest.branch <- min(tree.tr$edge.length)
#     smallest.branches[list.element] <- smallest.branch
#   }
# }
# min(smallest.branches, na.rm = TRUE) #
# ## seeds and transmission network sizes:
# # 2>>10,  3>>212,  8>>20, 10>>10, 12>>3, 16>>13, 19>>1741
# # which(smallest.branches!="NA")



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

  if(nrow(as.data.frame(tree.n)) >= 3){

    # Construct transmission trees
    tree.i <- trans.network2tree(transnetwork = tree.n)
    num.trees <- c(num.trees,tree.n$id[1])
    num.i <- c(num.i,i)

    tree.j <- tree.i
    tree.j$tip.label <- paste(i,".", tree.j$tip.label, ".A", sep = "")

    # Save the transmission tree
    write.tree(tree.j, file = paste("tree.model1.seed",i,".nwk", sep = ""))

    # tr <- read.tree(file = paste("tree.model1.seed",i,".nwk", sep = "")


    # Keep sampling dates
    id.samplingtime <- as.data.frame(cbind(paste(i,".",tree.n$id, ".A", sep = ""), tree.n$dtimes)) # IDs and their samling times in the transmission network

    write.csv(id.samplingtime,file=paste("samplingtimes_seed_",i,".csv", sep = ""))

  }
}

IDs.transm <- num.i # vector of seeds with at least 2 transmissions



### Binding all transmission trees together ###
###############################################

# Make a list of transmission trees

trees <- list() # list of all transmission trees

for(j in 1:length(IDs.transm)){
  p <- IDs.transm[j]
  tr <- read.tree(file = paste("tree.model1.seed",p,".nwk", sep = ""))
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
    colitem <- read.csv(paste("samplingtimes_seed_",i,".csv", sep = ""))
  }
  else{

    ritem <- read.csv(paste("samplingtimes_seed_",i,".csv", sep = ""))
    colitem <- rbind(colitem, ritem)
  }
}


# Prepare to run the sequence simulation

seq.rand <- 1

n.tr <- 1

seed.id <- 123

# # call the seed sequences - pool of viruses and rename the file
file.copy(paste("hiv.seq.B.pol.j.fasta", sep = ""),paste("seed.seq.bis.nwk", sep = ""))
# add the number of tree in the file and
write(n.tr,file = paste("seed.seq.bis.nwk",sep = ""), append = TRUE)  # n.tr
# the tree, to prepare the file to simulate the evolution of the virus across the tree

write.tree(resolved.combined.tree,file = paste("seed.seq.bis.nwk", sep = ""), append = TRUE)

file.rename(from = paste("seed.seq.bis.nwk", sep = ""), to = paste("seed.seq.bis.sim.nwk", sep = ""))

system(paste("./seq-gen -mGTR -f 0.3857, 0.1609, 0.2234, 0.2300  -a 0.9 -g 4 -i 0.5230  -r 2.9114, 12.5112, 1.2569, 0.8559, 12.9379, 1.0000 -s 0.00475  -n1 -k",seq.rand,"< seed.seq.bis.sim.nwk -z",seed," > B.Epidemic_seed.seq.bis.sim.nwk.fasta",sep = ""))

# a: shape parameter of Gamma > Gamma Rate Heterogeneity
# g: category of Gamma > Discrete Gamma Rate Heterogeneity
# r: rate matrix
# s: scale which is the substitution rate of pol gene
# z: seed


##################################################################################################
# Step 4: Construct time stamped phylogenetic tree of the epidemic with simulated sequences data #
##################################################################################################


# Function to tranform dates in named vector to be handled by treedater

dates.Transform.NamedVector  <- function(dates=dates){

  dates.val <- datalist$itable$population.simtime[1] - dates$V2 + 1977 # dates
  names(dates.val) <- as.character(dates$V1) # names are the names of the tips

  return(dates.val)
}

# 4.1. Construct phylogenetic tree

# Many tools to build the trees: within R like ape, phangorn, and outside compiled tools like iq-tree and FastTree

system(paste("./FastTree  -nt <", paste("B.Epidemic_seed.seq.bis.sim.nwk.fasta"), paste(">B.Epidemic_seed.seq.bis.sim.nwk.fasta.tree", sep = "")))


samp.dates <- colitem # call the data table of dates

tree.const <- read.tree(paste("B.Epidemic_seed.seq.bis.sim.nwk.fasta.tree", sep = ""))

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

dater.tree <- dater(tree.const, Ord.tree.dates, s = 3012) # s is the length of sequence

# save(dater.tree, file = paste("dated.tree.object.Rdata", sep = ""))

# test.t.dated <- dater.tree

dater.tree <- get(load("dated.tree.object.Rdata"))

N <- node.age(dater.tree)

int.node.age <- N$Ti # internal nodes ages

latest.samp <- N$timeToMRCA+N$timeOfMRCA # latest sampling date


### Lineage Through Time ###
############################

# Estimating confidence intervals for rates and dates using a parametric bootstrap

pb <- parboot.treedater(dater.tree,
                        nreps = 100,
                        level = .95 ) # Lineage Through Time

# plot.parboot.ltt( pb ) # This function computes the lineages through time given bootstrap replicate trees.
# The pseudo-maximum likelihood estimate is plotted alongside CIs based on bootstrap trees.

# Function to compute lineages through time with confidence intervals

plot.parboot.ltt.dat <- function (pbtd, t0 = NA, res = 100, ...)
{
  t1 <- max(pbtd$td$sts, na.rm = T)
  if (is.na(t0))
    t0 <- min(sapply(pbtd$trees, function(tr) tr$timeOf))
  times <- seq(t0, t1, l = res)
  ltt <- cbind(times = times, t(sapply(times, function(t) {
    c(pml = sum(pbtd$td$sts > t) - sum(pbtd$td$Ti > t), setNames(quantile(sapply(pbtd$trees,
                                                                                 function(tre) sum(tre$sts > t) - sum(tre$Ti > t)),
                                                                          probs = c(0.025, 0.5, 0.975)), c("lb", "median",
                                                                                                           "ub")))
  })))
  pl.df <- as.data.frame(ltt)
  return(pl.df)
  # p <- ggplot(pl.df) + geom_ribbon(aes(x = times, ymin = lb,
  #                                      ymax = ub), fill = "blue", col = "blue", alpha = 0.1)
  # p <- p + geom_path(aes(x = times, y = pml))
  # (p <- p + ylab("Lineages through time") + xlab("Time"))
}

LTT <- plot.parboot.ltt.dat(pb)



###################################################
# Step 5: Visualise and Use of simulation outputs #
###################################################



## Count transmission events per each transmission network ##
#############################################################


trans.net <- simpact.trans.net

transmissions.vec.i <- list()
time.calendar.vec.i <- list()

for (i in 1:length(IDs.transm)){

  p <- IDs.transm[i]

  trans.net.i <- simpact.trans.net[[p]]

  Infec.time.i <- trans.net.i$InfecTime + 1987

  min.val = 1987
  max.val = round(max(trans.net.i$itimes)) + 1987

  step.int=1
  d <- (max.val-min.val)/step.int


  dat.f.trans.i <- as.data.frame(trans.net.i)
  dat.f.trans.i$itimes <- abs(dat.f.trans.i$itimes-40)+1977

  numb.tra <- vector()
  i.vec <- vector()
  for (j in 1:d){
    inf <- 1986+j
    sup <- 1987+j
    dat.f.trans.count <- dat.f.trans.i[which(dat.f.trans.i$itimes <= sup & dat.f.trans.i$itimes  >= inf),]
    numb.i <- nrow(dat.f.trans.count)
    numb.tra <- c(numb.tra, numb.i)
    i.vec <- c(i.vec, sup)
  }

  transmissions.vec.i[[i]] <- numb.tra
  time.calendar.vec.i[[i]] <- i.vec

}

## Summing transmission events in all transmission networks ##
##############################################################

for (i in 1:length(IDs.transm)){
  if(i==1){
    trans.sum <- transmissions.vec.i[[i]]
  }
  else{

    read.trans.sum <- transmissions.vec.i[[i]]
    trans.sum <- trans.sum + read.trans.sum
  }
}



## Summing internal nodes in same age interval ##
#################################################

dt.node.age.dt <- int.node.age

min.val = 1987
max.val = 2017

step.int=1
d <- (max.val-min.val)/step.int

int.node.vec <- vector()
for (i in 1:d) {
  inf <- 1986+i
  sup <- 1987+i
  int.node.age.i <- int.node.age[int.node.age <= sup & dt.node.age.dt >= inf]
  int.node.vec <- c(int.node.vec,length(int.node.age.i))
}


# Transmissison events and internal nodes within same time intervals

numb.tra <- trans.sum
int.node.vec <- int.node.vec



#### PLOTTING #####


## Transmission network plotting: Union of all transmission networks ##
#######################################################################

network.list <- list()

for (i in 1:length(IDs.transm)){

  p <- IDs.transm[i]

  trans.net.i <- as.data.frame(simpact.trans.net[[p]])

  trans.net.i <- trans.net.i[-1,]

  trans.net.i$id <- paste("A.",p,".",trans.net.i$id, sep = "")
  trans.net.i$parent <- paste("A.",p,".",trans.net.i$parent, sep = "")

  graph.build <- trans.net.i

  graph.build[,4] <- as.character(graph.build$parent) # donors
  graph.build[,3] <- as.character(graph.build$id) # recipients
  gag = as.matrix(graph.build)
  ga.graph = graph.edgelist(gag[,4:3])

  V(ga.graph)$color <- "red"

  # transNet.yrs.Old <- delete.vertices(ga.graph, "-1")

  network.list[[i]] <- ga.graph

}


for (i in 1:length(IDs.transm)){
  if(i==1){
    trans.network.union <- network.list[[i]]
  }
  else{

    read.trans.network.union <- network.list[[i]]
    trans.network.union <- trans.network.union + read.trans.network.union
  }
}


# 1. Transmission network from simpact                           # 1 #
plot.igraph(trans.network.union, edge.arrow.size=0.1, vertex.size=5,
            edge.color="black",
            asp = 1,
            xlim = c(-1, 2),
            ylim = c(-0.5,0.5),
            vertex.frame.color="black",
            vertex.label.color="black",
            #vertex.label = NULL,
            layout = layout_with_kk,
            edge.width = 1,
            vertex.label.cex=0.1,
            vertex.label.dist=0.0, edge.curved=0.0, asp=0, margin=0
            #main = "True transmission network"
)

# 2. Phylogenetic tree

plot(dater.tree, show.tip.label=FALSE,
     edge.width=1,
     edge.color="blue") # Try a few different settings!
axisPhylo(backward = FALSE)


# 3. Transmission events and internal nodes
x <- i.vec
plot(x, int.node.vec, type="b", col="red", lwd=2,
     xlab = "Calendar time",
     ylab = "Count") # 1 > 1
lines(x, numb.tra, col='green3', type='b', lwd=2)
legend("topleft", legend = c("Internal nodes", "Transmission events"),
       col=c("red","green3"), pch=1)

# 4. Plot of lineages through time with confidence intervals

plot.parboot.ltt( pb )

# Or

pl.df <- LTT

p <- ggplot(pl.df) + geom_ribbon(aes(x = times, ymin = lb,
                                     ymax = ub), fill = "blue", col = "blue", alpha = 0.1)
p <- p + geom_path(aes(x = times, y = pml))
(p <- p + ylab("Lineages through time") + xlab("Calendar time"))


# Object for plotting by ggplot
SimpactPaperPhyloExample <- list()
SimpactPaperPhyloExample$transm.network <- trans.network.union
SimpactPaperPhyloExample$dater.tree <- dater.tree
SimpactPaperPhyloExample$years.vec <- i.vec
SimpactPaperPhyloExample$int.node.vec <- int.node.vec
SimpactPaperPhyloExample$numb.trans.vec <- numb.tra
save(SimpactPaperPhyloExample, file = "SimpactPaperPhyloExample.RData")
