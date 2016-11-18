#Read the libraries
rm(list=ls())
pacman::p_load(RSimpactCyan, RSimpactHelper, data.table, dplyr, magrittr, exactci,
               nlme, ggplot2,survival, KMsurv, tidyr, expoTree, sna, intergraph,
               igraph,lhs, GGally, emulator, multivator, tidyr, psych)

## mem_used() Check what the memory is. library(pryr)
## gc(TRUE)
## sort(sapply(mget(ls()),object.size)) size of objects that are in the environment

dirname <- "/home/trust/Documents/GIT_Projects/RSimpactHelp/"# getwd()
## Simple Simpact emulator

all.chunks.processed <- c("SummaryOutPut-inANDout.df.chunk-1-1000-2016-11-05-2.csv",
                          "SummaryOutPut-inANDout.df.chunk-1-1000-2016-11-05-3.csv")

file.name.csv <- paste0(dirname, "SummaryOutPut-inANDout.df.chunk-1-1000-2016-11-05-1.csv") # param.varied
# Read the output file from running simpact many times.
inputANDoutput.complete <- read.csv(file = file.name.csv, header = TRUE)
for (chunk.df in all.chunks.processed){
  chuck.df.filename <- paste0(dirname, chunk.df) # param.varied
  chunk.df.data <- read.csv(chuck.df.filename)
  chunk.df.data <- data.frame(chunk.df.data)

  inputANDoutput.complete <- rbind(inputANDoutput.complete, chunk.df.data)
}

## House keeping to know how many of the variables are being tracked.
summary.count <- which(colnames(inputANDoutput.complete)=="sim.id") - 1
xdesign.count <- length(names(dplyr::select(inputANDoutput.complete, contains("xdesign"))))
simpact.par.count <- length(inputANDoutput.complete) - summary.count - xdesign.count - 1

repetition.count <- nrow(inputANDoutput.complete[inputANDoutput.complete$sim.id==1,])


# The x variables (model parameters) that were varied:
x.variables <- varied.simpact.params()[[1]]


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
inputANDoutput.select <- dplyr::filter(inputANDoutput.complete,
                                     complete.cases(inputANDoutput.complete[,z.variables]),
                                     prev.men.15.25 > 0.05)

# remove all rows which do not have sim.id occurring the number of repeat times.
freq.sim.id <- plyr::count(inputANDoutput.select,"sim.id")
#get all the average in all the columns in the selected df
inputANDoutput.select <- aggregate(inputANDoutput.select, by = list(inputANDoutput.select$sim.id), FUN = "mean")
inputANDoutput.select <- left_join(freq.sim.id, inputANDoutput.select, by = "sim.id")
inputANDoutput.select <- inputANDoutput.select[inputANDoutput.select$freq == repetition.count,]

#Select the first 250 testing the optiomal_paras()
#You can also select a fraction of simulated dataset set round(dim(simpact.inANDout.df)[1]*0.10, digits=0)
simpact.inANDout.df <- head(inputANDoutput.select, 250)

### Check is z.variables and x.variables are in inANDout.df ####

try(if(length(targets)!=length(z.variables)) stop("Target values are not equal to the variables set"))

#select the x model param values (model parameters)
simpact.x <- dplyr::select_(inputANDoutput.select,.dots=x.variables) %>% as.matrix()
#select the z model param values (summary statistics)
simpact.z <- dplyr::select_(inputANDoutput.select,.dots=z.variables)

#Creating the Latin Hypercube Sample (LHS) for each of the parameters
x.design.name <- names(dplyr::select(inputANDoutput.select, contains("xdesign")))
x.design <- dplyr::select_(inputANDoutput.select,.dots=x.design.name)

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
#which.min(SumSq)
sim.best.summary <- z.df[which.min(SumSq), ]
#### And most importantly, the best estimate for the model parameters:
x.estimate <- as.numeric(simpact.x[which.min(SumSq), ])

p.stats <- rbind(sim.best.summary, targets)
rownames(p.stats) <- c("Summary Statistics","Targets")

p.stats

#Do principle componet on the summary statistics
z.pc <- princomp(z.df, scores = TRUE, cor = TRUE)
cum.prop.var <- summary(z.pc) # The first 4 components capture 94% of the total variance. That's perhaps sufficient?
comp.number.to.use <- 4 #change this if 5 which captures 98% is better.
plot(z.pc)
#biplot(z.pc)
#z.pc$loadings
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
RS.mdm <- mdm(x.design.long, types = rep(z.variables, each = dim(simpact.z)[1]))
RS.expt <- experiment(mm = RS.mdm, obs = z_obs)

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
as.numeric(as.numeric(z.pc$loadings[, 1])) %*% ((as.numeric(as.data.frame(targets) - z.pc$center) / z.pc$scale) )


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


simpact4emulation.pred <- function(sim.id, rep.sim.id){

  for (j in x.variables){
    assign.cfg.value <- lhs.pred.df[sim.id,j]
    cfg.pred[j][[1]] <- assign.cfg.value
  }

  ## Set those parameters that need to be the same
  cfg.pred$person.eagerness.woman.dist.gamma.a <- cfg.pred$person.eagerness.man.dist.gamma.a
  cfg.pred$person.eagerness.woman.dist.gamma.b <- cfg.pred$person.eagerness.man.dist.gamma.b
  cfg.pred$formation.hazard.agegapry.numrel_woman <- cfg.pred$formation.hazard.agegapry.numrel_man
  cfg.pred$formation.hazard.agegapry.gap_factor_man_exp <- cfg.pred$formation.hazard.agegapry.gap_factor_woman_exp

  simpact.seed.id <- rep.sim.id

  sub.dir.sim.id <- sprintf("%06d",sim.id)
  sub.dir.rep.sim.id <- sprintf("%02d",rep.sim.id)

  sub.dir.rename <- paste("temp/",sub.dir.sim.id,"/",sub.dir.rep.sim.id,sep = "")

  testoutput <- simpact.run(configParams = cfg.pred, destDir = "temp", agedist = agedist.data.frame, intervention = iv,
                            identifierFormat = paste0("%T-%y-%m-%d-%H-%M-%S_%p_%r%r%r%r%r%r%r%r_",sub.dir.sim.id,"_",sub.dir.rep.sim.id,"-"),
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

      datalist.test <- tryCatch(simpact4emulation.pred(sim.id,j), error = errFunction)
      if(length(datalist.test)>1){
        #get the summary statistics for each run
        out.test <- output.summary.maker(datalist.test, growth.rate=list(timewindow =c(0,20)),
                     agemix.maker=list(agegroup=c(15,30), timepoint =30, timewindow = 1, start=FALSE, gender = "female"),
                     prev.15.25 = list(age.group=c(15,25), timepoint = 35, gender = "men"),
                     prev.25.50 = list(age.group=c(25,50), timepoint = 35, gender = "men"),
                     art.coverage = list(age.group=c(15,50), timepoint = 35, gender = "men"),
                     inc.15.30 = list(age.group=c(15,30), timewindow=c(30,40), gender = "women", only.active = "No"),
                     partner.degree = list(age.group=c(15,30), hivstatus = 0, survey.time = 30,
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

    outStats.df$sim.id <- sim.id
    repeat.sum.stats.df <- rbind(repeat.sum.stats.df,outStats.df)

    big.insert <- length(inANDout.pred.20.df) - length(z.variables) + 1

    inANDout.pred.20.df[sim.id, big.insert:length(inANDout.pred.20.df)] <- out.test
    inANDout.pred.20.df$success.pred.rows[sim.id] <- success.pred.rows

    write.csv(inANDout.pred.20.df, file =paste("RowUpdate_pred_20.","-",design.points,"Points",variables,"Par_Partial",Sys.Date(), ".csv", sep=""), row.names = FALSE)
    write.csv(repeat.20.sum.stats.df, file =paste("RepeatAverage.20.df","-",design.points,"Points",variables,"Par_Partial",Sys.Date(), ".csv", sep=""), row.names = FALSE)
    #unlink(paste(getwd(),"/temp/*", sep=""))

    #write the data stating the number of design points and the number of var parameters
  }

  write.csv(inANDout.pred.20.df, file =paste("inANDout.pred.20.df","-",validation.repeats,"Points",
                                             length(z.variables),"Par",Sys.Date(), ".csv", sep=""), row.names = FALSE)
  write.csv(repeat.20.sum.stats.df, file =paste("RepeatAverage.20.df","-",design.points,"Points",variables,"Par",Sys.Date(),
                                             ".csv", sep=""), row.names = FALSE)
  sim.output.result <- list(inANDout.pred.20.df, repeat.20.sum.stats.df)

  return(sim.output.result)


}



start.time = proc.time()
sim.20.output.result <- succInANDOut.pred.df(validation.repeats)
end.time = proc.time() - start.time

inANDout.pred.20.df <- sim.20.output.result[[1]]
repeat.sum.stats.df <- sim.20.output.result[[2]]


small.insert <- length(inANDout.pred.20.df) - 1

inANDout.df <-
repeat.sum.stats.df <- simpact.inANDout.df[[1]][2]


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


