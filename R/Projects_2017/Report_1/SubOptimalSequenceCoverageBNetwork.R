# ## For sub-optimal sequence coverage
#

setwd("/home/david/RSimpactHelp/R/Projects_2017/Report_1/")


## Load required packages

pacman::p_load(devtools, Rcpp, ape, expoTree, data.table, phylosim, RSimpactCyan,
               RSimpactHelper, readr, phangorn, Biostrings, dplyr, adephylo,
               phyclust, DECIPHER,treedater,geiger,picante)


datalist <- get(load("MasterModelSubOptimalSeqCovearge.datalistB.RData"))

# simpact.trans.net <- transmNetworkBuilder.diff2(datalist = datalist, endpoint = 40)




# Sampling dates in calender time
dates.Transform.NamedVector  <- function(dates=dates){
  dates.val <- 1977+40-as.numeric(dates$V2) # dates datalist$itable$population.simtime[1] - dates$V2 + 1977
  names(dates.val) <- as.character(dates$V1) # names are the names of the tips
  return(dates.val)
}




setwd("/home/david/Dropbox/Analysis_for_Manuscripts_December2017_Sequence_Coverage")

sim.data <-  get(load("MasterModelSubOptimalSeqCovearge.datalistB.RData"))

trans.net.data <- get(load("simpact.trans.net.B.RData"))

trans.net.data.12 <- get(load("filtered.trans.net.seed12.Rdata"))
trans.net.data.22 <- get(load("filtered.trans.net.seed22.Rdata"))
trans.net.data.24 <- get(load("filtered.trans.net.seed24.Rdata"))

transmit_rate <- (nrow(trans.net.data.12) + nrow(trans.net.data.22) + nrow(trans.net.data.24))/7


###########################################################
# Scenario 1 Analysis: Network - one seed - one subtype A #
###########################################################

setwd("/home/david/Dropbox/ANALYSIS_NOVEMBER_2017/Analysis1_B/Complete_Transmission_Seed_22/Network_1/")


net <- get(load("simpact.trans.net.B.RData"))
net.22 <- net[[22]]

transmit_rate <- (length(net.22$id))/30

dated.tree.A.seed.22.seq.20.606 <- get(load("dated.tree.A.seed.22.seq.20.606.Rdata"))
write.tree(dated.tree.A.seed.22.seq.20.606, file = "dated.tree.A.seed.22.seq.20.606.nwk")
t.20 <- read.tree("dated.tree.A.seed.22.seq.20.606.nwk")
remove_rate.20 <-  (length(t.20$tip.label))/30

dated.tree.A.seed.22.seq.30.910 <- get(load("dated.tree.A.seed.22.seq.30.910.Rdata"))
write.tree(dated.tree.A.seed.22.seq.30.910, file = "dated.tree.A.seed.22.seq.30.910.nwk")
t.30 <- read.tree("dated.tree.A.seed.22.seq.30.910.nwk")
remove_rate.30 <-  (length(t.30$tip.label))/30

dated.tree.A.seed.22.seq.40.1213 <- get(load("dated.tree.A.seed.22.seq.40.1213.Rdata"))
write.tree(dated.tree.A.seed.22.seq.40.1213, file = "dated.tree.A.seed.22.seq.40.1213.nwk")
t.40 <- read.tree("dated.tree.A.seed.22.seq.40.1213.nwk")
remove_rate.40 <-  (length(t.40$tip.label))/30

dated.tree.A.seed.22.seq.50.1516 <- get(load("dated.tree.A.seed.22.seq.50.1516.Rdata"))
write.tree(dated.tree.A.seed.22.seq.50.1516, file = "dated.tree.A.seed.22.seq.50.1516.nwk")
t.50 <- read.tree("dated.tree.A.seed.22.seq.50.1516.nwk")
remove_rate.50 <-  (length(t.50$tip.label))/30

dated.tree.A.seed.22.seq.60.1819 <- get(load("dated.tree.A.seed.22.seq.60.1819.Rdata"))
write.tree(dated.tree.A.seed.22.seq.60.1819, file = "dated.tree.A.seed.22.seq.60.1819.nwk")
t.60 <- read.tree("dated.tree.A.seed.22.seq.60.1819.nwk")
remove_rate.60 <-  (length(t.60$tip.label))/30

dated.tree.A.seed.22.seq.70.2122 <- get(load("dated.tree.A.seed.22.seq.70.2122.Rdata"))
write.tree(dated.tree.A.seed.22.seq.70.2122, file = "dated.tree.A.seed.22.seq.70.2122.nwk")
t.70 <- read.tree("dated.tree.A.seed.22.seq.70.2122.nwk")
remove_rate.70 <-  (length(t.70$tip.label))/30

dated.tree.A.seed.22.seq.80.2426 <- get(load("dated.tree.A.seed.22.seq.80.2426.Rdata"))
write.tree(dated.tree.A.seed.22.seq.80.2426, file = "dated.tree.A.seed.22.seq.80.2426.nwk")
t.80 <- read.tree("dated.tree.A.seed.22.seq.80.2426.nwk")
remove_rate.80 <-  (length(t.80$tip.label))/30

dated.tree.A.seed.22.seq.90.2426 <- get(load("dated.tree.A.seed.22.seq.90.2729.Rdata"))
write.tree(dated.tree.A.seed.22.seq.90.2426, file = "dated.tree.A.seed.22.seq.90.2729.nwk")
t.90 <- read.tree("dated.tree.A.seed.22.seq.90.2729.nwk")
remove_rate.90 <-  (length(t.90$tip.label))/30

dated.tree.A.seed.22.seq.100.3032 <- get(load("dated.tree.A.seed.22.seq.100.3032.Rdata"))
write.tree(dated.tree.A.seed.22.seq.100.3032, file = "dated.tree.A.seed.22.seq.100.3032.nwk")
t.100 <- read.tree("dated.tree.A.seed.22.seq.100.3032.nwk")
remove_rate.100 <-  (length(t.100$tip.label))/30

##########################################################
# Scenario 2 Analysis: Network - 3 seeds - one subtype A #
##########################################################


setwd("/home/david/Dropbox/ANALYSIS_NOVEMBER_2017/Analysis2_B/7yrs_past_all_seeds/Network_2/")

dated.tree.A.seed.12.22.24seq.p20.381 <- get(load("dated.tree.A.seed.12.22.24seq.p20.381.Rdata"))
write.tree(dated.tree.A.seed.12.22.24seq.p20.381, file = "dated.tree.A.seed.12.22.24seq.p20.381.nwk")
t.20 <- read.tree("dated.tree.A.seed.12.22.24seq.p20.381.nwk")
remove_rate.20 <-  (length(t.20$tip.label))/7

dated.tree.A.seed.12.22.24seq.p30.571 <- get(load("dated.tree.A.seed.12.22.24seq.p30.571.Rdata"))
write.tree(dated.tree.A.seed.12.22.24seq.p30.571, file = "dated.tree.A.seed.12.22.24seq.p30.571.nwk")
t.30 <- read.tree("dated.tree.A.seed.12.22.24seq.p30.571.nwk")
remove_rate.30 <-  (length(t.30$tip.label))/7

dated.tree.A.seed.12.22.24seq.p40.762 <- get(load("dated.tree.A.seed.12.22.24seq.p40.762.Rdata"))
write.tree(dated.tree.A.seed.12.22.24seq.p40.762, file = "dated.tree.A.seed.12.22.24seq.p40.762.nwk" )
t.40 <- read.tree("dated.tree.A.seed.12.22.24seq.p40.762.nwk")
remove_rate.40 <-  (length(t.40$tip.label))/7

dated.tree.A.seed.12.22.24seq.p50.952 <- get(load("dated.tree.A.seed.12.22.24seq.p50.952.Rdata"))
write.tree(dated.tree.A.seed.12.22.24seq.p50.952, file = "dated.tree.A.seed.12.22.24seq.p50.952.nwk")
t.50 <- read.tree("dated.tree.A.seed.12.22.24seq.p50.952.nwk")
remove_rate.50 <-  (length(t.50$tip.label))/7

dated.tree.A.seed.12.22.24seq.p60.1142 <- get(load("dated.tree.A.seed.12.22.24seq.p60.1142.Rdata"))
write.tree(dated.tree.A.seed.12.22.24seq.p60.1142, file = "dated.tree.A.seed.12.22.24seq.p60.1142.nwk")
t.60 <- read.tree("dated.tree.A.seed.12.22.24seq.p60.1142.nwk")
remove_rate.60 <-  (length(t.60$tip.label))/7

dated.tree.A.seed.12.22.24seq.p70.1333 <- get(load("dated.tree.A.seed.12.22.24seq.p70.1333.Rdata"))
write.tree(dated.tree.A.seed.12.22.24seq.p70.1333, file = "dated.tree.A.seed.12.22.24seq.p70.1333.nwk")
t.70 <- read.tree("dated.tree.A.seed.12.22.24seq.p70.1333.nwk")
remove_rate.70 <-  (length(t.70$tip.label))/7

dated.tree.A.seed.12.22.24seq.p80.1523 <- get(load("dated.tree.A.seed.12.22.24seq.p80.1523.Rdata"))
write.tree(dated.tree.A.seed.12.22.24seq.p80.1523, file = "dated.tree.A.seed.12.22.24seq.p80.1523.nwk")
t.80 <- read.tree("dated.tree.A.seed.12.22.24seq.p80.1523.nwk")
remove_rate.80 <-  (length(t.80$tip.label))/7

dated.tree.A.seed.12.22.24Seq.p90.1714 <- get(load("dated.tree.A.seed.12.22.24Seq.p90.1714.Rdata"))
write.tree(dated.tree.A.seed.12.22.24Seq.p90.1714, file = "dated.tree.A.seed.12.22.24Seq.p90.1714.nwk")
t.90 <- read.tree("dated.tree.A.seed.12.22.24Seq.p90.1714.nwk")
remove_rate.90 <-  (length(t.90$tip.label))/7

dated.tree.A.seed.12.22.24seq.p100.1904 <- get(load("dated.tree.A.seed.12.22.24seq.p100.1904.Rdata"))
write.tree(dated.tree.A.seed.12.22.24Seq.p100.1904, file = "dated.tree.A.seed.12.22.24Seq.p100.1904.nwk")
t.100 <- read.tree("dated.tree.A.seed.12.22.24Seq.p100.1904.nwk")
remove_rate.100 <-  (length(t.100$tip.label))/7


#############################################################
# Scenario 3 Analysis: Network - 3 seeds - 2 subtypes A & G #
#############################################################

setwd("/home/david/Dropbox/ANALYSIS_NOVEMBER_2017/Analysis3_B/7_yrs_all_seeds_diff_subtype/Network_3/")

dated.tree.G.A.G.seed.12.22.24seq.p20.381 <- get(load("dated.tree.G.A.G.seed.12.22.24seq.p20.381.Rdata"))
write.tree(dated.tree.G.A.G.seed.12.22.24seq.p20.381, file = "dated.tree.G.A.G.seed.12.22.24seq.p20.381.nwk")

dated.tree.G.A.G.seed.12.22.24seq.p30.571 <- get(load("dated.tree.G.A.G.seed.12.22.24seq.p30.571.Rdata"))
write.tree(dated.tree.G.A.G.seed.12.22.24seq.p30.571, file = "dated.tree.G.A.G.seed.12.22.24seq.p30.571.nwk")

dated.tree.G.A.G.seed.12.22.24seq.p40.762 <- get(load("dated.tree.G.A.G.seed.12.22.24seq.p40.762.Rdata"))
write.tree(dated.tree.G.A.G.seed.12.22.24seq.p40.762, file = "dated.tree.G.A.G.seed.12.22.24seq.p40.762.nwk")

dated.tree.G.A.G.seed.12.22.24seq.p50.952 <- get(load("dated.tree.G.A.G.seed.12.22.24seq.p50.952.Rdata"))
write.tree(dated.tree.G.A.G.seed.12.22.24seq.p50.952, file = "dated.tree.G.A.G.seed.12.22.24seq.p50.952.nwk")

dated.tree.G.A.G.seed.12.22.24seq.p60.1142 <- get(load("dated.tree.G.A.G.seed.12.22.24seq.p60.1142.Rdata"))
write.tree(dated.tree.G.A.G.seed.12.22.24seq.p60.1142, file = "dated.tree.G.A.G.seed.12.22.24seq.p60.1142.nwk")

dated.tree.G.A.G.seed.12.22.24seq.p70.1333 <- get(load("dated.tree.G.A.G.seed.12.22.24seq.p70.1333.Rdata"))
write.tree(dated.tree.G.A.G.seed.12.22.24seq.p70.1333, file = "dated.tree.G.A.G.seed.12.22.24seq.p70.1333.nwk")

dated.tree.G.A.G.seed.12.22.24seq.p80.1523 <- get(load("dated.tree.G.A.G.seed.12.22.24seq.p80.1523.Rdata"))
write.tree(dated.tree.G.A.G.seed.12.22.24seq.p80.1523, file = "dated.tree.G.A.G.seed.12.22.24seq.p80.1523.nwk")

dated.tree.G.A.G.seed.12.22.24seq.p90.1714 <- get(load("dated.tree.G.A.G.seed.12.22.24seq.p90.1714.Rdata"))
write.tree(dated.tree.G.A.G.seed.12.22.24seq.p90.1714, file = "dated.tree.G.A.G.seed.12.22.24seq.p90.1714.nwk")

dated.tree.G.A.G.seed.12.22.24seq.p100.1904 <- get(load("dated.tree.G.A.G.seed.12.22.24seq.p100.1904.Rdata"))
write.tree(dated.tree.G.A.G.seed.12.22.24seq.p100.1904, file = "dated.tree.G.A.G.seed.12.22.24seq.p100.1904.nwk")

