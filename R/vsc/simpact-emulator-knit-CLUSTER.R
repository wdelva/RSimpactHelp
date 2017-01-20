#set the packages that we need
pacman::p_load(data.table, dplyr, magrittr, exactci,
               nlme, ggplot2,survival, KMsurv, tidyr, expoTree, sna, intergraph,
               igraph,lhs, GGally, emulator, multivator, tidyr, psych)

pred.starttime.All.time <- proc.time()

file.name.csv <- paste0("/","SummaryOutPut-inANDout.df.chunk-PCA-emu2017-01-19.csv") # param.varied
# Read the output file from running simpact many times.
inputANDoutput.complete <- data.frame(read.csv(file = file.name.csv, header = TRUE))

summary.count <- which(colnames(inputANDoutput.complete)=="sim.id") - 1
xdesign.count <- length(names(dplyr::select(inputANDoutput.complete, contains("xdesign"))))
simpact.par.count <- length(inputANDoutput.complete) - summary.count - xdesign.count - 1

#you can fix this when you know it otherwise let the line below determine
repetition.count <- nrow(inputANDoutput.complete[inputANDoutput.complete$sim.id==1,])

names(inputANDoutput.complete[,(summary.count+xdesign.count+2):length(inputANDoutput.complete)])
summaryparameters <- names(inputANDoutput.complete[,1:summary.count])

#summary statistics
z.variables <- summaryparameters
#Set the targets for the summary statistics.
targets <- c(0.015, 0.016, 0.043, 0.21, 0.47, 0.37, 0.54, 0.33, 0.34, 5)

#target.variables <- c("growth.rate", "inc.men.20.25", "inc.wom.20.25", "prev.men.25.30",
#                      "prev.wom.25.30","prev.men.30.35", "prev.wom.30.35", "ART.cov.men.18.50",
#                      "ART.cov.wom.18.50", "median.wom.18.50.AD")

try(if(length(targets)!=length(z.variables)) stop("Target values are not equal to the variables set"))

# The x variables (model parameters) that were varied:
x.variables <- dplyr::select(inputANDoutput.complete, contains("."))
x.variables <- x.variables[, !(colnames(x.variables) %in% z.variables)]
x.variables <- names(x.variables[2:length(x.variables)])

#get all the average in all the columns in the selected df
inputANDoutput.select <- aggregate(inputANDoutput.complete, by = list(inputANDoutput.complete$sim.id), FUN = "mean")

#### This will generate pairs that with each parameter
inputANDoutput.select$is.complete <- complete.cases(inputANDoutput.select)
sum(inputANDoutput.select$is.complete)
pair.plot <- pairs(inputANDoutput.select[, x.variables[1:8]],
                   col = 1+inputANDoutput.select$is.complete, pch = 16, cex = 2)

## Remove all rows with NA in the summary statistic
inputANDoutput.select <- dplyr::filter(inputANDoutput.select, complete.cases(inputANDoutput.select[,z.variables]))

#impose some conditions on the target statistics as needed, you can adjust as needed
##inputANDoutput.select <- dplyr::filter(inputANDoutput.select, growth.rate > 0)

#Keep the full set of the data to use for validation later
inputANDoutput.selectTTE <- inputANDoutput.select

#Select part of the data for use on multivator - emulator
nrow.sel <- round(nrow(inputANDoutput.select) * 85/100,0) # use 85% of the data always and use the 15% for validation)
inputANDoutput.select <- head(inputANDoutput.select, nrow.sel)

#select the x model param values (model parameters)
simpact.x <- dplyr::select_(inputANDoutput.select,.dots=x.variables) %>% as.matrix()
#select the z model param values (summary statistics)
simpact.z <- dplyr::select_(inputANDoutput.select,.dots=z.variables)
#select the x.design frame
x.design.name <- names(dplyr::select(inputANDoutput.select, contains("xdesign")))
x.design <- dplyr::select_(inputANDoutput.select,.dots=x.design.name)

##The figure below shows the distribution of the summary data to be used.
#If these divert from normalist PCA is the best option as it corrects for normality
#par(mfrow=c(1,1))
simpact.z.plot <- multi.hist(simpact.z)

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


########################################### Creating the multivator objects for a non-PCA-based analysis
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

# Use the loop to iterate through different values. So the optimisation process is faster.
err.function.optim.param <- function(e){
  if (length(grep("computationally singular",e$message)) != 0)
  stop(e)
}

for (iter in seq(100,700, 100)){
  print (paste("Working on option b iteration number: ", iter, sep=" "))
  RS.opt.b.var.iter <- optimal_params(RS.expt, option="b", start_hp = RS.opt.var, control = list(maxit=iter))

  #tryCatch(optimal_params(RS.expt, option="b", start_hp = RS.opt.var, control = list(maxit=iter)),
   #                             error = err.function.optim.param)

  RS.opt.var <- RS.opt.b.var.iter

  comp1.B.var <- data.frame(t(diag(B(RS.opt.b.var.iter)[,,1])),paste("b",iter,sep = ""))
  comp2.B.var <- data.frame(t(diag(B(RS.opt.b.var.iter)[,,2])),paste("b",iter,sep = ""))

  names(comp1.B.var)[length(comp1.B.var)] <- "run.type"
  names(comp2.B.var)[length(comp1.B.var)] <- "run.type"

  comp1.B <- rbind(comp1.B, comp1.B.var)
  comp2.B <- rbind(comp2.B, comp2.B.var)

}
#check how long this took.
RS.opt.b <- RS.opt.var

optim.check.conv <- proc.time() - optim.check

 #See the plot of convergency in the B matrix of coefficients. -->
gg.B1.plot <- ggplot(melt(comp1.B,id.vars=c("run.type")), aes(x=run.type,y=value, color=variable)) + geom_point() + ggtitle("Coeff 1")
gg.B2.plot <- ggplot(melt(comp2.B,id.vars=c("run.type")), aes(x=run.type,y=value, color=variable)) + geom_point() + ggtitle("Coeff 2")

for (iter in seq(100,700, 100)){
  print (paste("Working on option c iteration number: ", iter, sep=" "))
  RS.opt.c.var.iter <- optimal_params(RS.expt, option="c", start_hp = RS.opt.var, control = list(maxit=iter))
  RS.opt.var <- RS.opt.c.var.iter

  comp1.B.var <- data.frame(t(diag(B(RS.opt.c.var.iter)[,,1])),paste("c",iter,sep = ""))
  comp2.B.var <- data.frame(t(diag(B(RS.opt.c.var.iter)[,,2])),paste("c",iter,sep = ""))

  names(comp1.B.var)[length(comp1.B.var)] <- "run.type"
  names(comp2.B.var)[length(comp1.B.var)] <- "run.type"

  comp1.B <- rbind(comp1.B, comp1.B.var)
  comp2.B <- rbind(comp2.B, comp2.B.var)
}
#check how long this took.
RS.opt.c <- RS.opt.var

optim.check.conv <- proc.time() - optim.check

#See the plot of convergency in the B matrix of coefficients. -->
par.conv.plot <- ggplot(melt(comp1.B,id.vars=c("run.type")),
                        aes(x=run.type,y=value, color=variable))
par.conv.plot <- par.conv.plot + geom_point() + ggtitle("Coeff 1") #change this to comp2 as needed.


############### Testing if the prediction of the unused data can be recreated ######################################
n.check <- nrow(inputANDoutput.selectTTE)-nrow.sel
x.new.check <- tail(subset(inputANDoutput.selectTTE, select=x.design.name), n.check)
x.new.long.check <- x.new.check[rep(1:nrow(x.new.check),length(z.variables)),]
x.new.long.check <- as.matrix(x.new.long.check)
RS.new.mdm.check <- mdm(rbind(x.new.long.check), types = rep(z.variables, each = n.check))

RS.opt.a.check <- multem(x = RS.new.mdm.check, expt = RS.expt, hp = RS.opt.a)
RS.a.df.check <- data.frame(matrix(RS.opt.a.check, nrow = n.check,
                                      dimnames = list(rownames = 1:n.check,
                                      colnames = paste0("a",names(z.df))[1:length(z.variables)])))

RS.opt.b.check <- multem(x = RS.new.mdm.check, expt = RS.expt, hp = RS.opt.b)
RS.b.df.check <- data.frame(matrix(RS.opt.b.check, nrow = n.check,
                                      dimnames = list(rownames = 1:n.check,
                                      colnames = paste0("b",names(z.df))[1:length(z.variables)])))

RS.opt.c.check <- multem(x = RS.new.mdm.check, expt = RS.expt, hp = RS.opt.c)
RS.c.df.check <- data.frame(matrix(RS.opt.c.check, nrow = n.check,
                                      dimnames = list(rownames = 1:n.check,
                                      colnames = paste0("c",names(z.df))[1:length(z.variables)])))

model.stats.check <- tail(subset(inputANDoutput.selectTTE, select=z.variables), n.check)
model.stats.check <- cbind(model.stats.check, RS.a.df.check, RS.b.df.check, RS.c.df.check)

### Visualise the results (Choose one of the summary statistics to visualise how thy compare)
#par(mfrow=c(1,1))
stats.compare <- dplyr::select(model.stats.check, contains(z.variables[4]))
est.plot <- matplot(stats.compare, pch = 20, cex = 2)
est.plot <- legend("topleft", colnames(stats.compare),col=seq_len(ncol(stats.compare)),
       cex=0.8,fill=seq_len(ncol(stats.compare)), bty = "n")

############## Using the Emulator to Explore the Parameter Space for the None PCA Part to get the statistics ###################
n <- 10000
set.seed(1)
x.new <- latin.hypercube(n, length(x.variables), names=colnames(x.design))
x.new.long <- x.new[rep(1:nrow(x.new),length(z.variables)),]
x.new.long <- as.matrix(x.new.long)
RS.new.mdm <- mdm(rbind(x.new.long), types = rep(z.variables, each = n))
######################################## Parameter space to interpolate summary statistics #########################################
pred.starttime <- proc.time()
RS.prediction.opt.a <- multem(x = RS.new.mdm, expt = RS.expt, hp = RS.opt.a)
pred.endtime <- proc.time() - pred.starttime

pred2.starttime <- proc.time()
RS.prediction.opt.b <- multem(x = RS.new.mdm, expt = RS.expt, hp = RS.opt.b)
pred2.endtime <- proc.time() - pred2.starttime

pred3.starttime <- proc.time()
RS.prediction.opt.c <- multem(x = RS.new.mdm, expt = RS.expt, hp = RS.opt.c)
pred3.endtime <- proc.time() - pred3.starttime

# par(mfrow=c(3,length(z.variables))) # distribution of z.variable i with opt.a
# for (i in 1:length(z.variables)){
#   hist(RS.prediction.opt.a[((i-1)*n+1):(i*n)], main = paste("Dist of ",z.variables[i], sep = " "), xlab = "opt.a")
# }
#
# for (i in 1:length(z.variables)){
#   hist(RS.prediction.opt.b[((i-1)*n+1):(i*n)], main = paste("Dist of ",z.variables[i], sep = " "), xlab = "opt.b")
# }
#
# for (i in 1:length(z.variables)){
#   hist(RS.prediction.opt.c[((i-1)*n+1):(i*n)], main = paste("Dist of ",z.variables[i], sep = " "), xlab = "opt.c")
# }

##One way of efficiently comparing emulation output with target statistics is to reshape RS.prediction.* as a dataframe
prediction.a.df <- data.frame(matrix(RS.prediction.opt.a, nrow = n,
                                       dimnames = list(rownames = 1:n, colnames = names(z.df)[1:length(z.variables)])))
prediction.b.df <- data.frame(matrix(RS.prediction.opt.b, nrow = n,
                                        dimnames = list(rownames = 1:n, colnames = names(z.df)[1:length(z.variables)])))
prediction.c.df <- data.frame(matrix(RS.prediction.opt.c, nrow = n,
                                     dimnames = list(rownames = 1:n, colnames = names(z.df)[1:length(z.variables)])))

check.me.a <- as.data.frame(t(apply(prediction.a.df, 1, function(x) (x - t(targets))^2)))
names(check.me.a) <- summaryparameters
check.me.b <- as.data.frame(t(apply(prediction.b.df, 1, function(x) (x - t(targets))^2)))
names(check.me.b) <- summaryparameters
check.me.c <- as.data.frame(t(apply(prediction.c.df, 1, function(x) (x - t(targets))^2)))
names(check.me.c) <- summaryparameters

check.plot.multem <- multi.hist(check.me.b)

sum.square.df.a <- transform(check.me.a, sum=rowSums(check.me.a))
sum.square.df.b <- transform(check.me.b, sum=rowSums(check.me.b))
sum.square.df.c <- transform(check.me.c, sum=rowSums(check.me.c))

new.xdesign.ssd.c <- data.frame(cbind(x.new, ssd.c=sum.square.df.c$sum))
new.xdesign.ssd.c <- new.xdesign.ssd.c[order(new.xdesign.ssd.c$ssd.c),]

write.csv(head(new.xdesign.ssd.c, nrow(inputANDoutput.select)), file =paste0("SummaryOutPut-inANDout.df.chunk-NONPCA-BEST-emu",Sys.Date(),".csv"), row.names = FALSE)

# par(mfrow=c(3,1))
# sum.square.df.a.sum <- hist(sum.square.df.a$sum, 100)
# sum.square.df.b.sum <- hist(sum.square.df.b$sum, 100)
# sum.square.df.c.sum <- hist(sum.square.df.c$sum, 100)
#
#
# #### Get the min in sum of squared distances between model statistics and target statistics
# predicted.values <- function(prediction.df, method){
#   sq.a.array <- as.data.frame(t(apply(prediction.df, 1, function(x) (x - t(targets))^2)))
#   names(sq.a.array) <- names(prediction.df)
#   SumSq <- as.numeric(rowSums(sq.array))
#   x.estimate.row <- which.min(SumSq)
#   pridicted.output.values <- cbind(x.estimate.row, prediction.df[x.estimate.row, ], method = method)
#   return(pridicted.output.values)
# }
# targets
# pred.a <- predicted.values(prediction.a.df, "a")
# pred.b <- predicted.values(prediction.b.df, "b")
# pred.c <- predicted.values(prediction.c.df, "c")
#
# targets.row <- data.frame(cbind("NA", t(targets), "target"))
# names(targets.row) <- names(pred.a)
#
# pred.all <- rbind(pred.a, pred.b, pred.c, targets.row)
#
# # #### And most importantly, the best estimate for the model parameters:
# x.estimate.a <- as.numeric(x.new[pred.a$x.estimate.row, ])
# x.estimate.b <- as.numeric(x.new[pred.b$x.estimate.row, ])
# x.estimate.c <- as.numeric(x.new[pred.c$x.estimate.row, ])
#
# ############ Compute the final estimated parameters to use with SIMPACT ##################
# #Select the config parameters that will be varied from the input config
# par.estimated <- data.frame(matrix(NA, nrow = 3, ncol = length(x.variables)+1))
# names(par.estimated) <- c(x.variables,"method")
#
# ###Give the boundaries for the parameters here;
# x.variables.boundaries <- simpact.params.boundaries()
#
# #Create the list of parameters with their min, max vlaue (all will sample from a unif distribution)
# x.index <- 0
# for (j in x.variables){
#   x.index <- x.index + 1
#   col.index <- which(colnames(par.estimated)==j)
#   if(is.null(x.variables.boundaries[[j]])){
#     #if the boundary of the parameter is not set
#     par.estimated[1,col.index] <- "No bound set"
#     par.estimated[2,col.index] <- "No bound set"
#     par.estimated[3,col.index] <- "No bound set"
#   }else{
#     #otherwise compute the estimated parameter value
#     min.var <- x.variables.boundaries[[j]][1]
#     max.var <- x.variables.boundaries[[j]][2]
#     par.estimated[1,col.index] <- qunif(x.estimate.a[x.index], min = as.numeric(min.var), max = as.numeric(max.var))
#     par.estimated[2,col.index] <- qunif(x.estimate.b[x.index], min = as.numeric(min.var), max = as.numeric(max.var))
#     par.estimated[3,col.index] <- qunif(x.estimate.c[x.index], min = as.numeric(min.var), max = as.numeric(max.var))
#   }
# }
#
# par.estimated[1,length(par.estimated)] <- "a"
# par.estimated[2,length(par.estimated)] <- "b"
# par.estimated[3,length(par.estimated)] <- "c"
#
# par.estimated  #The estimated parameters

predpc.starttime.FIN.time <-  proc.time() - pred.starttime.All.time

save.image("~/MaxART/RSimpactHelp/data/NONEPCA-emulator-run2017-01-20-Cluster-WIM.RData")

################################################ END #####################################################







