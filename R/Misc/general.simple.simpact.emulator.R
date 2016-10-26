#Read the libraries
rm(list=ls())
pacman::p_load(RSimpactCyan, RSimpactHelper, data.table, dplyr, magrittr, exactci,
               nlme, ggplot2,survival, KMsurv, tidyr, expoTree, sna, intergraph,
               igraph,lhs, GGally, emulator, multivator, tidyr)

## mem_used() Check what the memory is. library(pryr)
## gc(TRUE)
## sort(sapply(mget(ls()),object.size)) size of objects that are in the environment

dirname <- getwd()
## Simple Simpact emulator

#file.name.csv <-"~RSimpactHelp/data/inANDout.df2016-09-27.csv" #6param.varied
#file.name.csv <-"~/RSimpactHelp/data/inANDout.df12016-10-01.csv" #2param.varied
#file.name.csv <-"~/RSimpactHelp/data/succRows.df.sim16.2016-09-28.csv" #2param.varied
file.name.csv <-paste(dirname, "/data/RowUpdate-2826Points10Par_Partial2016-10-19.csv", sep = "") #10param.varied

#Read the output file from running simpact many times.
inANDout.df <- read.csv(file = file.name.csv, header = TRUE, sep = ",")

# The x variables (model parameters) that were varied:
#x.variables = c("eagerness.a", "eagerness.b")
x.variables <- c("person.eagerness.man.dist.gamma.a", "person.eagerness.man.dist.gamma.b", "conception.alpha_base",
                 "formation.hazard.agegapry.numrel_man", "formation.hazard.agegapry.eagerness_diff",
                 "formation.hazard.agegapry.gap_factor_man_exp", "person.agegap.man.dist.normal.mu",
                 "person.agegap.woman.dist.normal.mu","person.agegap.man.dist.normal.sigma",
                 "person.agegap.woman.dist.normal.sigma")
x.variables.boundaries <- list(person.eagerness.man.dist.gamma.a.min =0.1, person.eagerness.man.dist.gamma.a.max =2,
                    person.eagerness.man.dist.gamma.b.min = 5, person.eagerness.man.dist.gamma.b.max = 60,
                    conception.alpha_base.min = -3.6, conception.alpha_base.max = -1.2,
                    formation.hazard.agegapry.numrel_man.min = -1.5, formation.hazard.agegapry.numrel_man.max = -0.1,
                    formation.hazard.agegapry.eagerness_diff.min = -0.1, formation.hazard.agegapry.eagerness_diff.max = 0,
                    formation.hazard.agegapry.gap_factor_woman_exp.min = -1.5, formation.hazard.agegapry.gap_factor_woman_exp.max =-0.4,
                    person.agegap.man.dist.normal.mu.min = 0, person.agegap.man.dist.normal.mu.max = 4,
                    person.agegap.woman.dist.normal.mu.min =0, person.agegap.woman.dist.normal.mu.max = 4,
                    person.agegap.man.dist.normal.sigma.min = 0.5, person.agegap.man.dist.normal.sigma.max =2,
                    person.agegap.woman.dist.normal.sigma.min =0.5, person.agegap.woman.dist.normal.sigma.max =2)

# The "z" variables (model summary statistics), for which we will define targets:
z.variables <- c("growth.rate", "median.AD", "Q1.AD", "Q3.AD", "prev.men.15.25", "prev.men.25.50",
                 "ART.cov.15.50")

#Set the targets for the summary statistics.
targets <- c(0.014, 3, 2, 5, 0.08, 0.25, 0.3)

#, "incid.wom.15.30", "frac.degreeGT1.wom.15.30"

## Decide if we want to keep only rows without NA or
simpact.inANDout.df <- dplyr::filter(inANDout.df,
                                     complete.cases(inANDout.df[,z.variables]),
                                     prev.15.25 > 0.05)

#Select the first 250 testing the optiomal_paras()
#You can also select a fraction of simulated dataset set round(dim(simpact.inANDout.df)[1]*0.10, digits=0)
simpact.inANDout.df <- head(simpact.inANDout.df, 250)

### Check is z.variables and x.variables are in inANDout.df ####

try(if(length(targets)!=length(z.variables)) stop("Target values are not equal to the variables set"))

#select the x model param values (model parameters)
simpact.x <- dplyr::select_(simpact.inANDout.df,.dots=x.variables) %>% as.matrix()
#select the z model param values (summary statistics)
simpact.z <- dplyr::select_(simpact.inANDout.df,.dots=z.variables)

#Creating the Latin Hypercube Sample (LHS) for each of the parameters
x.design.name <- as.character(c(2:(length(x.variables)+1)))
x.design <- dplyr::select_(simpact.inANDout.df,.dots=x.design.name)
for( i in 1:length(x.variables)) {
  colnames(x.design)[i] <- paste("x",i,".design", sep="")
}


sum(duplicated(simpact.x)) # Should be zero

##Check the distribution of the summary statistics output
par(mfrow=c(1,1))
multi.hist(simpact.z)
#multi.hist(simpact.x)

##### Need to transform the none normal statistics from the multi.hist(simpact.z)

##Creating a LHS for each summary statistic
x.design.long <- x.design[rep(1:nrow(x.design),length(z.variables)),]
x.design.long <- as.matrix(x.design.long)
#### Computing model output at design points
z_obs <- as.vector(unlist(simpact.z))

#### Before we start the emulation, let's see what the best fit is from the model simulation runs
z.df <- simpact.z
sq.array <- as.data.frame(t(apply(z.df, 1, function(x) (x - t(targets))^2)))
SumSq <- as.numeric(rowSums(sq.array))
which.min(SumSq)
z.df[which.min(SumSq), ]
targets
#### And most importantly, the best estimate for the model parameters:
x.estimate <- as.numeric(simpact.x[which.min(SumSq), ])

#Do principle componet on the summary statistics
z.pc <- princomp(z.df, scores = TRUE, cor = TRUE)
summary(z.pc) # The first 4 components capture 94% of the total variance. That's perhaps sufficient?
comp.number.to.use <- 4 #change this if 5 which captures 98% is better.
plot(z.pc)
biplot(z.pc)
z.pc$loadings
z.pc.df <- data.frame(z.pc$scores)
z.pc.obs <- as.vector(unlist(z.pc.df[ ,1:comp.number.to.use]))

## Decide if you will drop some of the summary statics
x.design.pc.long <- x.design[rep(1:nrow(x.design),comp.number.to.use),]
x.design.pc.long <- as.matrix(x.design.pc.long)


#### Creating the multivator objects for the PCA-based analysis
RS.pc.mdm <- mdm(x.design.pc.long, types = rep(names(z.pc.df)[1:comp.number.to.use], each = dim(simpact.z)[1]))
RS.pc.expt <- experiment(mm = RS.pc.mdm, obs = z.pc.obs)

optima.starttime.pc <- proc.time()
RS.pc.opt.a <- optimal_params(RS.pc.expt, option="a")
optima.endtime.pc <- proc.time() - optima.starttime.pc

comp1.pc.B <- data.frame(t(diag(B(RS.pc.opt.a)[,,1])),"a")
comp2.pc.B <- data.frame(t(diag(B(RS.pc.opt.a)[,,2])),"a")
comp3.pc.B <- data.frame(t(diag(B(RS.pc.opt.a)[,,3])),"a")
comp4.pc.B <- data.frame(t(diag(B(RS.pc.opt.a)[,,4])),"a")

names(comp1.pc.B)[names(comp1.pc.B)=="X.a."] <- "run.type"
names(comp2.pc.B)[names(comp2.pc.B)=="X.a."] <- "run.type"
names(comp3.pc.B)[names(comp3.pc.B)=="X.a."] <- "run.type"
names(comp4.pc.B)[names(comp4.pc.B)=="X.a."] <- "run.type"

RS.opt.pc.var <-  RS.pc.opt.a

optim.check.pc <- proc.time()

## Use the loop to get iterate through different values. So the optimasation process is faster.

for (iter in seq(100,700, 100)){
  print (paste("Working on iteration number: ", iter, sep=" "))
  RS.opt.b.var.iter <- optimal_params(RS.pc.expt, option="b", start_hp = RS.opt.pc.var, control = list(maxit=iter))
  RS.opt.pc.var <- RS.opt.b.var.iter

  comp1.pc.B.var <- data.frame(t(diag(B(RS.opt.b.var.iter)[,,1])),paste("b",iter,sep = ""))
  comp2.pc.B.var <- data.frame(t(diag(B(RS.opt.b.var.iter)[,,2])),paste("b",iter,sep = ""))
  comp3.pc.B.var <- data.frame(t(diag(B(RS.opt.b.var.iter)[,,3])),paste("b",iter,sep = ""))
  comp4.pc.B.var <- data.frame(t(diag(B(RS.opt.b.var.iter)[,,4])),paste("b",iter,sep = ""))

  names(comp1.pc.B.var)[11] <- "run.type"
  names(comp2.pc.B.var)[11] <- "run.type"
  names(comp3.pc.B.var)[11] <- "run.type"
  names(comp4.pc.B.var)[11] <- "run.type"

  comp1.pc.B <- rbind(comp1.pc.B, comp1.pc.B.var)
  comp2.pc.B <- rbind(comp2.pc.B, comp2.pc.B.var)
  comp3.pc.B <- rbind(comp3.pc.B, comp3.pc.B.var)
  comp4.pc.B <- rbind(comp4.pc.B, comp4.pc.B.var)
}
#check how long this took.
RS.pc.opt.b <- RS.opt.pc.var
optim.check.pc.conv <- proc.time() - optim.check.pc

#See the plot of convergency in the B matrix of coefficients.
ggplot(melt(comp1.pc.B,id.vars=c("run.type")), aes(x=run.type,y=value, color=variable)) + geom_point() + facet_grid(variable ~ .)
ggplot(melt(comp2.pc.B,id.vars=c("run.type")), aes(x=run.type,y=value, color=variable)) + geom_point()
ggplot(melt(comp3.pc.B,id.vars=c("run.type")), aes(x=run.type,y=value, color=variable)) + geom_point() + facet_grid(variable ~ .)
ggplot(melt(comp4.pc.B,id.vars=c("run.type")), aes(x=run.type,y=value, color=variable)) + geom_point()


#### Creating the multivator objects for a non-PCA-based analysis
#RS_mdm <- mdm(x.design.long, types = rep(c("growth.rate", "prev.15.50.end"), each = design.points))
RS.mdm <- mdm(x.design.long, types = rep(z.variables, each = dim(simpact.z)[1]))
RS.expt <- experiment(mm = RS.mdm, obs = z_obs)
#RS_opt <- optimal_params(RS_expt, option="a", verbose = TRUE)

optima.starttime <- proc.time()
RS.opt.a <- optimal_params(RS.expt, option="a")
optima.endtime <- proc.time() - optima.starttime



comp1.B <- data.frame(t(diag(B(RS.opt.a)[,,1])),"a")
comp2.B <- data.frame(t(diag(B(RS.opt.a)[,,2])),"a")
comp3.B <- data.frame(t(diag(B(RS.opt.a)[,,3])),"a")
comp4.B <- data.frame(t(diag(B(RS.opt.a)[,,4])),"a")

names(comp1.B)[names(comp1.B)=="X.a."] <- "run.type"
names(comp2.B)[names(comp2.B)=="X.a."] <- "run.type"
names(comp3.B)[names(comp3.B)=="X.a."] <- "run.type"
names(comp4.B)[names(comp4.B)=="X.a."] <- "run.type"

RS.opt.var <-  RS.opt.a

optim.check <- proc.time()

## Use the loop to get iterate through different values. So the optimasation process is faster.

for (iter in seq(100,700, 100)){
  print (paste("Working on iteration number: ", iter, sep=" "))
  RS.opt.b.var.iter <- optimal_params(RS.expt, option="b", start_hp = RS.opt.var, control = list(maxit=iter))
  RS.opt.var <- RS.opt.b.var.iter

  comp1.B.var <- data.frame(t(diag(B(RS.opt.b.var.iter)[,,1])),paste("b",iter,sep = ""))
  comp2.B.var <- data.frame(t(diag(B(RS.opt.b.var.iter)[,,2])),paste("b",iter,sep = ""))
  comp3.B.var <- data.frame(t(diag(B(RS.opt.b.var.iter)[,,3])),paste("b",iter,sep = ""))
  comp4.B.var <- data.frame(t(diag(B(RS.opt.b.var.iter)[,,4])),paste("b",iter,sep = ""))

  names(comp1.B.var)[11] <- "run.type"
  names(comp2.B.var)[11] <- "run.type"
  names(comp3.B.var)[11] <- "run.type"
  names(comp4.B.var)[11] <- "run.type"

  comp1.B <- rbind(comp1.B, comp1.B.var)
  comp2.B <- rbind(comp2.B, comp2.B.var)
  comp3.B <- rbind(comp3.B, comp3.B.var)
  comp4.B <- rbind(comp4.B, comp4.B.var)
}
#check how long this took.
RS.opt.b <- RS.opt.var

optim.check.conv <- proc.time() - optim.check

#See the plot of convergency in the B matrix of coefficients.
ggplot(melt(comp1.B,id.vars=c("run.type")), aes(x=run.type,y=value, color=variable)) + geom_point()
ggplot(melt(comp2.B,id.vars=c("run.type")), aes(x=run.type,y=value, color=variable)) + geom_point()
ggplot(melt(comp3.B,id.vars=c("run.type")), aes(x=run.type,y=value, color=variable)) + geom_point()
ggplot(melt(comp4.B,id.vars=c("run.type")), aes(x=run.type,y=value, color=variable)) + geom_point()


## Plot the cov matrix from the 100, 200, 300 results (just to visualise how they differ)
## could perfom another test ??
##plotcov(M(RS.opt.b.100),M(RS.opt.b.200))
## corr.mat.b.100 <- M(RS.opt.b.100)/sqrt(tcrossprod(diag(M(RS.opt.b.100))))
## cor.plot(corr.mat.b.100, numbers = TRUE)
#### Using the emulator to explore the parameter space on the PCA method
n <- 10000
x.pc.new <- latin.hypercube(n, length(x.variables), names=colnames(x.design))
x.pc.new.long <- x.pc.new[rep(1:nrow(x.pc.new),comp.number.to.use),]
x.pc.new.long <- as.matrix(x.pc.new.long)

RS.pc.new.mdm <- mdm(rbind(x.pc.new.long), types = rep(names(z.pc.df)[1:comp.number.to.use], each = n))

pred.starttime.pc.mult <- proc.time()
RS.pc.prediction.opt.a <- multem(x = RS.pc.new.mdm, expt = RS.pc.expt, hp = RS.pc.opt.a)
pred.endtime.pc.mult <- proc.time() - pred.starttime.pc.mult

pred2.starttime.pc.mult <- proc.time()
RS.pc.prediction.opt.b <- multem(x = RS.pc.new.mdm, expt = RS.pc.expt, hp = RS.pc.opt.b)
pred2.endtime.pc.mult <- proc.time() - pred2.starttime.pc.mult

par(mfrow=c(2,comp.number.to.use))# distribution of pc i
for (i in 1:comp.number.to.use){
  hist(RS.pc.prediction.opt.a[((i-1)*n+1):(i*n)], main = paste("Dist of PC",i, sep = " "), xlab = "opt.pc.a")
}

par(mfrow=c(2,comp.number.to.use))
for (i in 1:comp.number.to.use){
  hist(RS.pc.prediction.opt.b[((i-1)*n+1):(i*n)], main = paste("Dist of PC",i, sep = " "), xlab = "opt.pc.b")
}

#### Using the emulator to explore the parameter space for the none PCA part
x.new <- latin.hypercube(n, length(x.variables), names=colnames(x.design))
x.new.long <- x.new[rep(1:nrow(x.new),length(z.variables)),]
x.new.long <- as.matrix(x.new.long)

RS.new.mdm <- mdm(rbind(x.new.long), types = rep(z.variables, each = n))

pred.starttime <- proc.time()
RS.prediction.opt.a <- multem(x = RS.new.mdm, expt = RS.expt, hp = RS.opt.a)
pred.endtime <- proc.time() - pred.starttime

pred2.starttime <- proc.time()
RS.prediction.opt.b <- multem(x = RS.new.mdm, expt = RS.expt, hp = RS.opt.b)
pred2.endtime <- proc.time() - pred2.starttime

par(mfrow=c(2,length(z.variables))) # distribution of z.variable i with opt.a
for (i in 1:length(z.variables)){
  hist(RS.prediction.opt.a[((i-1)*n+1):(i*n)], main = paste("Dist of z.var",i, sep = " "), xlab = "opt.a")
}

par(mfrow=c(2,length(z.variables))) # distribution of z.variable i with opt.b
for (i in 1:length(z.variables)){
  hist(RS.prediction.opt.b[((i-1)*n+1):(i*n)], main = paste("Dist of z.var",i, sep = " "), xlab = "opt.b")
}


##One way of efficiently comparing emulation output with target statistics is to reshape RS.pc.prediction.* as a dataframe
prediction.pc.a.df <- data.frame(matrix(RS.pc.prediction.opt.a, nrow = n,
                                     dimnames = list(rownames = 1:n, colnames = names(z.pc.df)[1:comp.number.to.use])))
prediction.pc.b.df <- data.frame(matrix(RS.pc.prediction.opt.b, nrow = n,
                                     dimnames = list(rownames = 1:n, colnames = names(z.pc.df)[1:comp.number.to.use])))
##One way of efficiently comparing emulation output with target statistics is to reshape RS.prediction.* as a dataframe
prediction.a.df <- data.frame(matrix(RS.prediction.opt.a, nrow = n,
                                        dimnames = list(rownames = 1:n, colnames = names(z.df)[1:length(z.variables)])))
prediction.b.df <- data.frame(matrix(RS.prediction.opt.b, nrow = n,
                                        dimnames = list(rownames = 1:n, colnames = names(z.df)[1:length(z.variables)])))

#### sum of squared distances between model statistics and target statistics
# Note: we could normalise (and centralise) these statistics to give them more equal weight in the SumSq

predicted.values <- function(prediction.df, pc = FALSE){
  if(pc == TRUE){
    targets <- predict(z.pc, as.data.frame(targets))
    targets <- as.numeric(targets)[1:4]
  }
  #
  sq.a.array <- as.data.frame(t(apply(prediction.df, 1, function(x) (x - t(targets))^2)))
  names(sq.a.array) <- names(prediction.df)
  SumSq <- as.numeric(rowSums(sq.array))
  x.estimate.row <- which.min(SumSq)
  pridicted.output.values <- cbind(x.estimate.row, prediction.df[x.estimate.row, ])
  return(pridicted.output.values)
}

#rbind(targets, predicted.values(prediction.a.df)[2], predicted.values(prediction.pc.a.df)[2])
targets
pred.a <- predicted.values(prediction.a.df)
pred.pc.a <- predicted.values(prediction.pc.a.df, pc = TRUE)
pred.b <- predicted.values(prediction.b.df)
pred.pc.b <- predicted.values(prediction.pc.b.df, pc = TRUE)


#### And most importantly, the best estimate for the model parameters:
x.estimate.a <- as.numeric(x.new[pred.a$x.estimate.row, ])
####x.estimate.pc.a <- as.numeric(x.pc.new[pred.pc.a$x.estimate.row, ])
x.estimate.b <- as.numeric(x.new[pred.b$x.estimate.row, ])
####x.estimate.pc.b <- as.numeric(x.pc.new[pred.pc.b$x.estimate.row, ])
# As an example: the value of the first PC for the target statistics:
# as.numeric(as.numeric(z.pc$loadings[, 1]) %*% ((as.numeric(as.dataframe(targets) - z.pc$center) / z.pc$scale) )


####################################################################

iv <- intervention.introduced(list(27,0,100,2),list(30,200,1.5), list(33,350,1),list(36,500,0.5))

validation.repeats <- 20
simulation.number.count <- 0
x.estimate <- x.estimate.b
rlhs.pred <- matrix(rep(x.estimate, each = validation.repeats), nrow = validation.repeats)

simpact.set.simulation("simpact-cyan")#("maxart") # Is it a standard or a MaxART simulation?
agedist.data.frame <- agedistr.creator(shape = 5, scale = 65)

#### Set input params
##Specifying the initially chosen values for the simulation.
cfg.pred <- input.params.creator(population.simtime = 40, population.numwomen = 500, population.nummen = 500)


#parameter list with input + default values
testoutput.headers <- simpact.run(configParams = cfg.pred, destDir = "temp", agedist = agedist.data.frame,
                                  intervention = iv, seed = 1, parallel = TRUE)
datalist.headers <- readthedata(testoutput.headers)
lhs.pred.df <- as.data.frame(head(datalist.headers$itable,1))
unlink(paste(getwd(),"/temp/*", sep=""))

#rm all the columns that have text + t
lhs.pred.df <- within(lhs.pred.df, rm(t, population.agedistfile, periodiclogging.outfile.logperiodic, logsystem.outfile.logevents,
                            logsystem.outfile.loglocation, logsystem.outfile.logpersons, logsystem.outfile.logrelations,
                            logsystem.outfile.logsettings, logsystem.outfile.logtreatments))

lhs.pred.df <- as.data.frame(lapply(lhs.pred.df, rep, validation.repeats))

x.index <- 0
for (j in x.variables){
  x.index <- x.index + 1

  min.var <- x.variables.boundaries[paste(j,".min",sep = "")][[1]]
  max.var <- x.variables.boundaries[paste(j,".max",sep = "")][[1]]
  col.index <- which(colnames(lhs.pred.df)==j)
  lhs.pred.df[col.index] <- qunif(x.estimate[x.index], min = as.numeric(min.var), max = as.numeric(max.var))
}


lhs.pred.df$formation.hazard.agegapry.numrel_woman <- lhs.pred.df$formation.hazard.agegapry.numrel_man
lhs.pred.df$formation.hazard.agegapry.gap_factor_man_exp <- lhs.pred.df$formation.hazard.agegapry.gap_factor_woman_exp
lhs.pred.df$success.rows <- NA

#Create a dataframe with NA for the summary statistics
summary.stats.pred.df <- data.frame(matrix(NA, nrow = validation.repeats, ncol = length(z.variables)))
names(summary.stats.pred.df) <- z.variables
# Creating a dataframe for input AND output
inANDout.pred.20.df <- cbind.data.frame(sim.id = 1:validation.repeats, rlhs.pred, lhs.pred.df, summary.stats.pred.df)


simpact4emulation.pred <- function(sim.id, lhs.df, cfg, agedist.data.frame, iv){

  for (j in x.variables){
    assign.cfg.value <- lhs.pred.df[sim.id,j]
    cfg.pred[j][[1]] <- assign.cfg.value
  }

  ## Set those parameters that need to be the same
  cfg.pred$person.eagerness.woman.dist.gamma.a <- cfg.pred$person.eagerness.man.dist.gamma.a
  cfg.pred$person.eagerness.woman.dist.gamma.b <- cfg.pred$person.eagerness.man.dist.gamma.b
  cfg.pred$formation.hazard.agegapry.numrel_woman <- cfg.pred$formation.hazard.agegapry.numrel_man
  cfg.pred$formation.hazard.agegapry.gap_factor_man_exp <- cfg.pred$formation.hazard.agegapry.gap_factor_woman_exp

  simpact.seed.id <- sim.id

  testoutput <- simpact.run(configParams = cfg.pred, destDir = "temp", agedist = agedist.data.frame, intervention = iv,
                            seed = simpact.seed.id)


  if (testoutput$simulationtime < cfg.pred$population.simtime)
  {
    if (testoutput$eventsexecuted >= cfg.pred$population.maxevents-1)
    {
      stop("MAXEVENTS: Simulation stopped prematurely, max events reached")
    }
    else
    {
      stop("Simulation stopped prematurely, probably ran out of events")
    }
  }

  datalist.test <- readthedata(testoutput)

}


succInANDOut.pred.df<- function(design.points=10){

  # Running the simpact4emulation.pred function with error catcher in a loop
  for (sim.id in inANDout.pred.20.df$sim.id){

    set.average.number <- 10
    #create a df to collect the repetead runs to be averaged
    outStats.df <- data.frame(matrix(NA, nrow = set.average.number, ncol = length(z.variables)))
    names(outStats.df) <- z.variables


    for (j in 1:set.average.number){
      simulation.number.count <- simulation.number.count + 1
      print(paste("Simulation Count:",simulation.number.count,"Design number:", sim.id,"/",validation.repeats, sep = " "))

      datalist.test <- tryCatch(simpact4emulation.pred(sim.id, lhs.df, cfg.pred, agedist.data.frame, iv), error = errFunction)
      if(length(datalist.test)>1){
        #get the summary statistics for each run
        out.test <- output.summary.maker(datalist.test, growth.rate=list(timewindow.min = 0, timewindow.max = 20),
                     agemix.maker=list(agegroup.min = 15, agegroup.max=30, timepoint =30,
                                       timewindow = 1, start=FALSE, gender = "female"),
                     prev.15.25 = list(age.group.min=15, age.group.max=25, timepoint = 35, gender = "men"),
                     prev.25.50 = list(age.group.min=25, age.group.max=50, timepoint = 35, gender = "men"),
                     art.coverage = list(age.group.min=15, age.group.max=50, timepoint = 35, gender = "men"),
                     inc.15.30 = list(age.group.min=15, age.group.max=30, timewindow.min = 30,
                                      timewindow.max = 40, gender = "women", only.active = "No"),
                     partner.degree = list(age.group.min=15, age.group.max=30, hivstatus = 0, survey.time = 30,
                                           window.width = 1, gender="female", only.new = FALSE))
        out.test <- out.test[,z.variables]
      }else{out.test <- rep(NA,length(z.variables))}

      outStats.df[j,] <- out.test

    }
    #average the ten runs to get the averal summary statistic

    success.pred.complete.df <- dplyr::filter(outStats.df, complete.cases(outStats.df))
    success.pred.rows <- nrow(success.pred.complete.df)
    out.test <- colMeans(outStats.df, na.rm = TRUE)
    # Inserting the output to the inANDout.pred.20.df dataframe

    big.insert <- length(inANDout.pred.20.df) - length(z.variables) + 1

    inANDout.pred.20.df[sim.id, big.insert:length(inANDout.pred.20.df)] <- out.test
    inANDout.pred.20.df$success.pred.rows[sim.id] <- success.pred.rows

    write.csv(inANDout.pred.20.df, file =paste("RowUpdate_pred_20.","-",design.points,"Points",variables,"Par_Partial",Sys.Date(), ".csv", sep=""), row.names = FALSE)
    unlink(paste(getwd(),"/temp/*", sep=""))

    #write the data stating the number of design points and the number of var parameters
  }

  write.csv(inANDout.pred.20.df, file =paste("inANDout.pred.20.df","-",validation.repeats,"Points",
                                             length(z.variables),"Par",Sys.Date(), ".csv", sep=""), row.names = FALSE)

  return(inANDout.pred.20.df)

}

start.time = proc.time()
inANDout.pred.20.df <- succInANDOut.pred.df(validation.repeats)
end.time = proc.time() - start.time


small.insert <- length(inANDout.pred.20.df) - 1

inANDout.pred.20.plusHarling.df <- rbind(inANDout.pred.20.df[, small.insert:length(inANDout.pred.20.df)],
                                         c(0.075, 0.013))
inANDout.pred.20.plusHarling.df$Harling <- c(rep(0, nrow(inANDout.pred.20.df)), 1)


if (nrow(inANDout.pred.20.df) > 0){
  pairs(inANDout.pred.20.plusHarling.df[, 1:2],
        col = 1+inANDout.pred.20.plusHarling.df$Harling)#,
  # xlim = c(0, 1),
  # ylim = c(0, 1))
  2
}


