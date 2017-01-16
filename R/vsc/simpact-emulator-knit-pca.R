#set the packages that we need
pacman::p_load(RSimpactCyan, RSimpactHelper, data.table, dplyr, magrittr, exactci,
               nlme, ggplot2,survival, KMsurv, tidyr, expoTree, sna, intergraph,
               igraph,lhs, GGally, emulator, multivator, tidyr, psych)
# install.packages("readcsvcolumns", repos="http://193.190.10.42/jori/")

comp <- "lin" #lin #mac

if(comp == "win"){
  dirname <- "~/MaxART/RSimpactHelp"
}else if(comp=="lin"){
  dirname <- "~/Documents/GIT_Projects/RSimpactHelp"
}else{
  dirname <- "~/Documents/RSimpactHelp"  #mac directory here
}

file.name.csv <- paste0(dirname, "/","SummaryOutPut-inANDout.df.chunk-1-250-2017-01-14.csv") # param.varied
# Read the output file from running simpact many times.
inputANDoutput.complete <- data.frame(read.csv(file = file.name.csv, header = TRUE))

summary.count <- which(colnames(inputANDoutput.complete)=="sim.id") - 1
xdesign.count <- length(names(dplyr::select(inputANDoutput.complete, contains("xdesign"))))
simpact.par.count <- length(inputANDoutput.complete) - summary.count - xdesign.count - 1

#you can fix this when you know it otherwise let the line below determine
repetition.count <- nrow(inputANDoutput.complete[inputANDoutput.complete$sim.id==1,])

names(inputANDoutput.complete[,(summary.count+xdesign.count+2):length(inputANDoutput.complete)])
summaryparameters <- names(inputANDoutput.complete[,1:summary.count])
summaryparameters

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

#x.variables.boundaries <- list(conception.alpha_base.min = -3.6, conception.alpha_base.max = -1.2,
#                              formation.hazard.agegapry.baseline.min = 1.5, formation.hazard.agegapry.baseline.max = 3)

#get all the average in all the columns in the selected df
inputANDoutput.select <- aggregate(inputANDoutput.complete, by = list(inputANDoutput.complete$sim.id), FUN = "mean")

## Remove all rows with NA in the summary statistic
inputANDoutput.select <- dplyr::filter(inputANDoutput.select,complete.cases(inputANDoutput.select[,z.variables]))

## some condition on the prev.men.15.50 e.g
##inputANDoutput.select <- dplyr::filter(inputANDoutput.select, prev.men.15.50 > 0.25)
inputANDoutput.select <- dplyr::filter(inputANDoutput.select, growth.rate > 0)

#Keep the full set of the data for testing later
inputANDoutput.selectTTE <- inputANDoutput.select

#Select a fraction of simulated dataset
nrow.sel <- floor(nrow(inputANDoutput.select) * 75/100) # use 75% of the data always and use the 25% for validation)
inputANDoutput.select <- head(inputANDoutput.select, nrow.sel)

#select the x model param values (model parameters)
simpact.x <- dplyr::select_(inputANDoutput.select,.dots=x.variables) %>% as.matrix()
#select the z model param values (summary statistics)
simpact.z <- dplyr::select_(inputANDoutput.select,.dots=z.variables)
#select the x.design frame
x.design.name <- names(dplyr::select(inputANDoutput.select, contains("xdesign")))
x.design <- dplyr::select_(inputANDoutput.select,.dots=x.design.name)

##The figure below shows the distribution of the summary data to be used.
##We might need to transform any of the summary statistics should the histogram diverge from normaility.
par(mfrow=c(1,1))
multi.hist(simpact.z)

##Creating a LHS for each summary statistic
x.design.long <- x.design[rep(1:nrow(x.design),length(z.variables)),]
x.design.long <- as.matrix(x.design.long)
#### Computing model output at design points
#z_obs <- as.vector(unlist(simpact.z))

############################ creating the PCA part of the data #####################################
z.df <- simpact.z
z.pc <- princomp(z.df, scores = TRUE, cor = TRUE)
z.pc.summary <- summary(z.pc)
#plot(z.pc)
#biplot(z.pc)
#z.pc$loadings
z.pc.cum.var <- cumsum((z.pc.summary$sdev)^2) / sum(z.pc.summary$sdev^2) #cummulative variance
pc.select.number <- which.min(abs(z.pc.cum.var - 0.96))[[1]] # Want 96% of the variance to be explained by ....
z.pc.df <- data.frame(z.pc$scores)

z.pc.obs <- as.vector(unlist(z.pc.df[ ,1:pc.select.number]))
x.design.pc.long <- x.design[rep(1:nrow(x.design),pc.select.number),]  #You can do length(z.variables) - #PCA not to be used
x.design.pc.long <- as.matrix(x.design.pc.long)

################ Creating the multivator objects for the PCA-based analysis
RS.pc.mdm <- mdm(x.design.pc.long, types = rep(names(z.pc.df)[1:pc.select.number], each = dim(simpact.z)[1])) #You can do names(z.pc.df)[1:2] - #PCA not to be used
RS.pc.expt <- experiment(mm = RS.pc.mdm, obs = z.pc.obs)

optima.starttime.pc <- proc.time()
RS.pc.opt.a <- optimal_params(RS.pc.expt, option="a")
optima.endtime.pc <- proc.time() - optima.starttime.pc

comp1.pc.B <- data.frame(t(diag(B(RS.pc.opt.a)[,,1])),"a")
comp2.pc.B <- data.frame(t(diag(B(RS.pc.opt.a)[,,2])),"a")

names(comp1.pc.B)[names(comp1.pc.B)=="X.a."] <- "run.type"
names(comp2.pc.B)[names(comp2.pc.B)=="X.a."] <- "run.type"

RS.pc.opt.var <-  RS.pc.opt.a

optim.pc.check <- proc.time()

# Use the loop to iterate through different values. So the optimisation process is faster.

for (iter in seq(100,700, 100)){
  print (paste("Working on pc iteration number: ", iter, sep=" "))
  RS.pc.opt.b.var.iter <- optimal_params(RS.pc.expt, option="b", start_hp = RS.pc.opt.var, control = list(maxit=iter))
  RS.pc.opt.var <- RS.pc.opt.b.var.iter

  comp1.pc.B.var <- data.frame(t(diag(B(RS.pc.opt.b.var.iter)[,,1])),paste("b",iter,sep = ""))
  comp2.pc.B.var <- data.frame(t(diag(B(RS.pc.opt.b.var.iter)[,,2])),paste("b",iter,sep = ""))

  names(comp1.pc.B.var)[length(comp1.pc.B.var)] <- "run.type"
  names(comp2.pc.B.var)[length(comp2.pc.B.var)] <- "run.type"

  comp1.pc.B <- rbind(comp1.pc.B, comp1.pc.B.var)
  comp2.pc.B <- rbind(comp2.pc.B, comp2.pc.B.var)
}
#check how long this took.
RS.pc.opt.b <- RS.pc.opt.var

optim.pc.check.conv <- proc.time() - optim.pc.check

#See the plot of convergency in the B matrix of coefficients. -->
ggplot(melt(comp1.pc.B,id.vars=c("run.type")), aes(x=run.type,y=value, color=variable)) + geom_point() + ggtitle("Coeff pc 1")
ggplot(melt(comp2.pc.B,id.vars=c("run.type")), aes(x=run.type,y=value, color=variable)) + geom_point() + ggtitle("Coeff pc 2")

for (iter in seq(100,700, 100)){
  print (paste("Working on pc iteration number: ", iter, sep=" "))
  RS.pc.opt.c.var.iter <- optimal_params(RS.pc.expt, option="c", start_hp = RS.pc.opt.var, control = list(maxit=iter))
  RS.pc.opt.var <- RS.pc.opt.c.var.iter

  comp1.pc.B.var <- data.frame(t(diag(B(RS.pc.opt.c.var.iter)[,,1])),paste("c",iter,sep = ""))
  comp2.pc.B.var <- data.frame(t(diag(B(RS.pc.opt.c.var.iter)[,,2])),paste("c",iter,sep = ""))

  names(comp1.pc.B.var)[length(comp1.pc.B.var)] <- "run.type"
  names(comp2.pc.B.var)[length(comp2.pc.B.var)] <- "run.type"

  comp1.pc.B <- rbind(comp1.pc.B, comp1.pc.B.var)
  comp2.pc.B <- rbind(comp2.pc.B, comp2.pc.B.var)
}
#check how long this took.
RS.pc.opt.c <- RS.pc.opt.var

optim.pc.check.conv <- proc.time() - optim.pc.check

#See the plot of convergency in the B matrix of coefficients. -->
ggplot(melt(comp1.pc.B,id.vars=c("run.type")), aes(x=run.type,y=value, color=variable)) + geom_point() + ggtitle("Coeff pc 1")
ggplot(melt(comp2.pc.B,id.vars=c("run.type")), aes(x=run.type,y=value, color=variable)) + geom_point() + ggtitle("Coeff pc 2")

############### Testing if the prediction of the unused data can be recreated ######################################
n.check <- nrow(inputANDoutput.selectTTE)-nrow.sel
x.new.check <- tail(subset(inputANDoutput.selectTTE, select=x.design.name), n.check)
x.new.long.check <- x.new.check[rep(1:nrow(x.new.check),pc.select.number),]
x.new.long.check <- as.matrix(x.new.long.check)
RS.new.mdm.check <- mdm(rbind(x.new.long.check), types = rep(names(z.pc.df)[1:pc.select.number], each = n.check))


RS.pc.opt.a.check <- multem(x = RS.new.mdm.check, expt = RS.pc.expt, hp = RS.pc.opt.a)
RS.pc.a.df.check <- data.frame(matrix(RS.pc.opt.a.check, nrow = n.check,
                                   dimnames = list(rownames = 1:n.check, colnames = paste0("a",names(z.pc.df))[1:pc.select.number])))

RS.pc.opt.b.check <- multem(x = RS.new.mdm.check, expt = RS.pc.expt, hp = RS.pc.opt.b)
RS.pc.b.df.check <- data.frame(matrix(RS.pc.opt.b.check, nrow = n.check,
                                   dimnames = list(rownames = 1:n.check, colnames = paste0("b",names(z.pc.df))[1:pc.select.number])))

RS.pc.opt.c.check <- multem(x = RS.new.mdm.check, expt = RS.pc.expt, hp = RS.pc.opt.c)
RS.pc.c.df.check <- data.frame(matrix(RS.pc.opt.c.check, nrow = n.check,
                                   dimnames = list(rownames = 1:n.check, colnames = paste0("c",names(z.pc.df))[1:pc.select.number])))

model.stats.check.pc <- tail(subset(inputANDoutput.selectTTE, select=z.variables), n.check)
model.stats.check.pc <- predict(z.pc, model.stats.check.pc)
model.stats.check.pc <- cbind(model.stats.check.pc[,1:pc.select.number], RS.pc.a.df.check, RS.pc.b.df.check, RS.pc.c.df.check)

### Visualise the results (Choose one of the summary statistics to visualise how they compare)
stats.compare.pc <- dplyr::select(model.stats.check.pc, contains("Comp.4"))
matplot(stats.compare.pc, pch = 20, cex = 2)
legend("topleft", colnames(stats.compare.pc),col=seq_len(ncol(stats.compare.pc)),cex=0.8,fill=seq_len(ncol(stats.compare.pc)), bty = "n")

############################ Using the Emulator to Explore the Parameter Space for the PCA Part to get the statistics
n<-10000
set.seed(1)
x.new <- latin.hypercube(n, length(x.variables), names=colnames(x.design))
x.new.long <- x.new[rep(1:nrow(x.new),pc.select.number),]
x.new.long <- as.matrix(x.new.long)
RS.new.mdm <- mdm(rbind(x.new.long), types = rep(names(z.pc.df)[1:pc.select.number], each = n))
#################################### Parameter space to interpolate summary statistics  ########################
predpc.starttime <- proc.time()
RS.prediction.pc.opt.a <- multem(x = RS.new.mdm, expt = RS.pc.expt, hp = RS.pc.opt.a)
predpc.endtime <- proc.time() - predpc.starttime

pred2pc.starttime <- proc.time()
RS.prediction.pc.opt.b <- multem(x = RS.new.mdm, expt = RS.pc.expt, hp = RS.pc.opt.b)
pred2pc.endtime <- proc.time() - pred2pc.starttime

pred3pc.starttime <- proc.time()
RS.prediction.pc.opt.c <- multem(x = RS.new.mdm, expt = RS.pc.expt, hp = RS.pc.opt.c)
pred3pc.endtime <- proc.time() - pred3pc.starttime

par(mfrow=c(3,pc.select.number)) # distribution of z.variable i with pc.opt.a then with pc.opt.b
for (i in 1:pc.select.number){
  hist(RS.prediction.pc.opt.a[((i-1)*n+1):(i*n)], main = paste("Dist of ",names(z.pc.df)[i], sep = " "), xlab = "opt.a")
}

for (i in 1:pc.select.number){
  hist(RS.prediction.pc.opt.b[((i-1)*n+1):(i*n)], main = paste("Dist of ",names(z.pc.df)[i], sep = " "), xlab = "opt.b")
}

for (i in 1:pc.select.number){
  hist(RS.prediction.pc.opt.c[((i-1)*n+1):(i*n)], main = paste("Dist of ",names(z.pc.df)[i], sep = " "), xlab = "opt.c")
}

##One way of efficiently comparing emulation output with target statistics is to reshape RS.prediction.* as a dataframe
prediction.pc.a.df <- data.frame(matrix(RS.prediction.pc.opt.a, nrow = n,
                                     dimnames = list(rownames = 1:n, colnames = names(z.pc.df)[1:pc.select.number])))
prediction.pc.b.df <- data.frame(matrix(RS.prediction.pc.opt.b, nrow = n,
                                     dimnames = list(rownames = 1:n, colnames = names(z.pc.df)[1:pc.select.number])))
prediction.pc.c.df <- data.frame(matrix(RS.prediction.pc.opt.c, nrow = n,
                                        dimnames = list(rownames = 1:n, colnames = names(z.pc.df)[1:pc.select.number])))

## Predicting the PC values against the targets
######as.numeric(as.numeric(z.pc$loadings[, 1]) %*% ((as.numeric(targets.df) - z.pc$center) / z.pc$scale) )
# All PCs for the target statistics:
targets.df <- data.frame(t(targets))
names(targets.df) <- z.variables
targets.pc <- predict(z.pc, targets.df)[,1:pc.select.number]
targets.pc.vector <- as.numeric(targets.pc)

predicted.values <- function(prediction.df, method){
  sq.a.array <- as.data.frame(t(apply(prediction.df, 1, function(x) (x - t(targets.pc.vector))^2)))
  names(sq.a.array) <- names(prediction.df)
  SumSq <- as.numeric(rowSums(sq.a.array))
  x.estimate.row <- which.min(SumSq)
  pridicted.output.values <- data.frame(cbind(x.estimate.row, prediction.df[x.estimate.row, ], method = method))
  return(pridicted.output.values)
}

pred.pc.a <- data.frame(predicted.values(prediction.pc.a.df, "a"))
pred.pc.b <- data.frame(predicted.values(prediction.pc.b.df, "b"))
pred.pc.c <- data.frame(predicted.values(prediction.pc.c.df, "c"))
targets.pc.vector.row <- data.frame(cbind("NA", t(targets.pc.vector), "target.pc"))
names(targets.pc.vector.row) <- names(pred.pc.a)

##pred.pc.all <- rbind(targets.pc.vector.row, pred.pc.a, pred.pc.b, pred.pc.c)

# #### And most importantly, the best estimate for the model parameters from the none PC part
x.estimate.pc.a <- as.numeric(x.new[pred.pc.a$x.estimate.row, ])
x.estimate.pc.b <- as.numeric(x.new[pred.pc.b$x.estimate.row, ])
x.estimate.pc.c <- as.numeric(x.new[pred.pc.c$x.estimate.row, ])

############ Compute the final estimated parameters to use with SIMPACT ##################
#Select the config parameters that will be varied from the input config
par.estimated.pc <- data.frame(matrix(NA, nrow = 3, ncol = length(x.variables)+1))
names(par.estimated.pc) <- c(x.variables,"method")

###Give the boundaries for the parameters here;
x.variables.boundaries <- simpact.params.boundaries()

#Create the list of parameters with their min, max vlaue (all will sample from a unif distribution)
x.index <- 0
for (j in x.variables){
  x.index <- x.index + 1
  col.index <- which(colnames(par.estimated.pc)==j)

  if(is.null(x.variables.boundaries[[j]])){
    #if the boundary of the parameter is not set
    par.estimated[1,col.index] <- NA
    par.estimated[2,col.index] <- NA
    par.estimated[3,col.index] <- NA
  }else{
    min.var <- x.variables.boundaries[[j]][1]
    max.var <- x.variables.boundaries[[j]][2]
    par.estimated.pc[1,col.index] <- qunif(x.estimate.pc.a[x.index], min = as.numeric(min.var), max = as.numeric(max.var))
    par.estimated.pc[2,col.index] <- qunif(x.estimate.pc.b[x.index], min = as.numeric(min.var), max = as.numeric(max.var))
    par.estimated.pc[3,col.index] <- qunif(x.estimate.pc.c[x.index], min = as.numeric(min.var), max = as.numeric(max.var))
  }
}

par.estimated.pc[1,length(par.estimated.pc)] <- "a"
par.estimated.pc[2,length(par.estimated.pc)] <- "b"
par.estimated.pc[3,length(par.estimated.pc)] <- "c"

par.estimated.pc  #The estimated parameters

################################################ END #####################################################



