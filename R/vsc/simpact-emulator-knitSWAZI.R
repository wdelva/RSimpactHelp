
pacman::p_load(RSimpactCyan, RSimpactHelper, data.table, dplyr, magrittr, exactci,
               nlme, ggplot2,survival, KMsurv, tidyr, expoTree, sna, intergraph,
               igraph,lhs, GGally, emulator, multivator, tidyr, psych)
# install.packages("readcsvcolumns", repos="http://193.190.10.42/jori/")

dirname <- "~/Documents/GIT_Projects/RSimpactHelp/" #### Change this as needed

file.name.csv <- paste0(dirname, "SummaryOutPut-inANDout.df.chunk-1-100-2016-12-01.csv") # param.varied
# Read the output file from running simpact many times.
inputANDoutput.complete <- data.frame(read.csv(file = file.name.csv, header = TRUE))

summary.count <- which(colnames(inputANDoutput.complete)=="sim.id") - 1
xdesign.count <- length(names(dplyr::select(inputANDoutput.complete, contains("xdesign"))))
simpact.par.count <- length(inputANDoutput.complete) - summary.count - xdesign.count - 1

repetition.count <- nrow(inputANDoutput.complete[inputANDoutput.complete$sim.id==1,])


names(inputANDoutput.complete[,(summary.count+xdesign.count+2):length(inputANDoutput.complete)])

summaryparameters <- names(inputANDoutput.complete[,1:summary.count])
summaryparameters


## We will only consider here those combinations that all produced valid summary statistics (no NA's) and those whose prev.men.15.25 > 0.05.


##The figure below shows the distribution of the summary data to be used. We might need to transform any of the summary statistics should the histogram diverge from normaility.

#summary statistics
z.variables <- c("growth.rate", "prev.men.15.50")

# The x variables (model parameters) that were varied:
x.variables <- c("conception.alpha_base","formation.hazard.agegapry.baseline")

x.variables.boundaries <- list(conception.alpha_base.min = -3.6, conception.alpha_base.max = -1.2,
                               formation.hazard.agegapry.baseline.min = 1.5, formation.hazard.agegapry.baseline.max = 3)

#select the


# remove all rows which do not have sim.id occurring the number of repeat times.
####freq.sim.id <- plyr::count(inputANDoutput.select,"sim.id")
#get all the average in all the columns in the selected df
#inputANDoutput.select <- inputANDoutput.select %>% dplyr::group_by(sim.id) %>% dplyr::summarise_each(funs(mean))

inputANDoutput.select <- aggregate(inputANDoutput.complete, by = list(inputANDoutput.complete$sim.id), FUN = "mean")


## Decide if we want to keep only rows without NA or
#inputANDoutput.select <- nrow(dplyr::filter(inputANDoutput.select,
#                                       complete.cases(inputANDoutput.complete[,z.variables]),
#                                       prev.men.15.50 > 0.25))


#inputANDoutput.select <- left_join(freq.sim.id, inputANDoutput.select, by = "sim.id")
##inputANDoutput.select <- inputANDoutput.select[inputANDoutput.select$freq == repetition.count,]


#Set the targets for the summary statistics.
targets <- c(0.014, 0.32)

#Select the first 250 testing the optiomal_paras()
#You can also select a fraction of simulated dataset set round(dim(simpact.inANDout.df)[1]*0.10, digits=0)
#inputANDoutput.select <- head(inputANDoutput.select, 250)

#select the x model param values (model parameters)
simpact.x <- dplyr::select_(inputANDoutput.select,.dots=x.variables) %>% as.matrix()
#select the z model param values (summary statistics)
simpact.z <- dplyr::select_(inputANDoutput.select,.dots=z.variables)
#select the x.design frame
x.design.name <- names(dplyr::select(inputANDoutput.select, contains("xdesign")))
x.design <- dplyr::select_(inputANDoutput.select,.dots=x.design.name)

### Check is z.variables and x.variables are in inANDout.df ####

try(if(length(targets)!=length(z.variables)) stop("Target values are not equal to the variables set"))


par(mfrow=c(1,1))
multi.hist(simpact.z)


##If we take the simulated data and ask what are the best baseline parameters that give statistics that are close to our target statistics. (We should be able to compare these with the procedural parameters that we would obtain).

##Creating a LHS for each summary statistic
x.design.long <- x.design[rep(1:nrow(x.design),length(z.variables)),]
x.design.long <- as.matrix(x.design.long)
#### Computing model output at design points
z_obs <- as.vector(unlist(simpact.z))

#### Before we start the emulation, let's see what the best fit is from the model simulation runs
z.df <- simpact.z
sq.array <- as.data.frame(t(apply(z.df, 1, function(x) (x - t(targets))^2)))
SumSq <- as.numeric(rowSums(sq.array))
sim.best.summary <- z.df[which.min(SumSq), ]
#### And most importantly, the best estimate for the model parameters:
x.estimate <- as.numeric(simpact.x[which.min(SumSq), ])

p.stats <- rbind(sim.best.summary, targets)
rownames(p.stats) <- c("Summary Statistics","Targets")

p.stats


##And the simpact parameters are:
x.estimate



### Creating the multivator objects for a non-PCA-based analysis
RS.mdm <- mdm(x.design.long, types = rep(z.variables, each = dim(simpact.z)[1]))
RS.expt <- experiment(mm = RS.mdm, obs = z_obs)

optima.starttime <- proc.time()
RS.opt.a <- optimal_params(RS.expt, option="a")
optima.endtime <- proc.time() - optima.starttime

comp1.B <- data.frame(t(diag(B(RS.opt.a)[,,1])),"a")
comp2.B <- data.frame(t(diag(B(RS.opt.a)[,,2])),"a")


names(comp1.B)[names(comp1.B)=="X.a."] <- "run.type"
names(comp2.B)[names(comp2.B)=="X.a."] <- "run.type"


RS.opt.var <-  RS.opt.a

optim.check <- proc.time()

# Use the loop to get iterate through different values. So the optimasation process is faster.

for (iter in seq(100,700, 100)){
  print (paste("Working on iteration number: ", iter, sep=" "))
  RS.opt.b.var.iter <- optimal_params(RS.expt, option="b", start_hp = RS.opt.var, control = list(maxit=iter))
  RS.opt.var <- RS.opt.b.var.iter

  comp1.B.var <- data.frame(t(diag(B(RS.opt.b.var.iter)[,,1])),paste("b",iter,sep = ""))
  comp2.B.var <- data.frame(t(diag(B(RS.opt.b.var.iter)[,,2])),paste("b",iter,sep = ""))


  names(comp1.B.var)[11] <- "run.type"
  names(comp2.B.var)[11] <- "run.type"

  comp1.B <- rbind(comp1.B, comp1.B.var)
  comp2.B <- rbind(comp2.B, comp2.B.var)

}
#check how long this took.
RS.opt.b <- RS.opt.var

optim.check.conv <- proc.time() - optim.check

 #See the plot of convergency in the B matrix of coefficients. -->
ggplot(melt(comp1.B,id.vars=c("run.type")), aes(x=run.type,y=value, color=variable)) + geom_point() + ggtitle("Coeff 1")
ggplot(melt(comp2.B,id.vars=c("run.type")), aes(x=run.type,y=value, color=variable)) + geom_point() + ggtitle("Coeff 2")


## Using the Emulator to Explore the Parameter Space for the None PCA Part

x.new <- latin.hypercube(n, length(x.variables), names=colnames(x.design))
x.new.long <- x.new[rep(1:nrow(x.new),length(z.variables)),]
x.new.long <- as.matrix(x.new.long)

RS.new.mdm <- mdm(rbind(x.new.long), types = rep(z.variables, each = n))

pred.starttime <- proc.time()
RS.prediction.opt.a <- multem(x = RS.new.mdm, expt = RS.expt, hp = RS.opt.a)
red.endtime <- proc.time() - pred.starttime

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


##One way of efficiently comparing emulation output with target statistics is to reshape RS.prediction.* as a dataframe
prediction.a.df <- data.frame(matrix(RS.prediction.opt.a, nrow = n,
                                       dimnames = list(rownames = 1:n, colnames = names(z.df)[1:length(z.variables)])))
prediction.b.df <- data.frame(matrix(RS.prediction.opt.b, nrow = n,
                                        dimnames = list(rownames = 1:n, colnames = names(z.df)[1:length(z.variables)])))

#### sum of squared distances between model statistics and target statistics
# Note: we could normalise (and centralise) these statistics to give them more equal weight in the SumSq

predicted.values <- function(prediction.df, pc = FALSE){

  sq.a.array <- as.data.frame(t(apply(prediction.df, 1, function(x) (x - t(targets))^2)))
  names(sq.a.array) <- names(prediction.df)
  SumSq <- as.numeric(rowSums(sq.array))
  x.estimate.row <- which.min(SumSq)
  pridicted.output.values <- cbind(x.estimate.row, prediction.df[x.estimate.row, ])
  return(pridicted.output.values)
}

rbind(targets, predicted.values(prediction.a.df)[2], predicted.values(prediction.pc.a.df)[2])
targets
pred.a <- predicted.values(prediction.a.df)

pred.b <- predicted.values(prediction.b.df)

# #### And most importantly, the best estimate for the model parameters:
x.estimate.a <- as.numeric(x.new[pred.a$x.estimate.row, ])

x.estimate.b <- as.numeric(x.new[pred.b$x.estimate.row, ])



