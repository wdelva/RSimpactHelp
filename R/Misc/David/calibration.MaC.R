# Calibration with MaC

## the model has two parameters and outputs two summary statistics.
## defining a simple toy model:

toy_model<-function(x){ c( x[1] + x[2] + rnorm(1,0,0.1) , x[1] * x[2] + rnorm(1,0,0.1) ) }

## define prior information
toy_prior=list(c("unif",0,1),c("normal",1,2))

# a uniform prior distribution between 0 and 1 for parameter 1, and a normal distribution
# of mean 1 and standard deviation of 2 for parameter 2.

## define the targeted summary statistics
sum_stat_obs=c(1.5,0.5)

## ABC
library(EasyABC)
## to perform the Marjoram et al. (2003)'s method:
##
ABC_Marjoram_original<-ABC_mcmc(method="Marjoram_original", model=toy_model, prior=toy_prior,
                                summary_stat_target=sum_stat_obs)
ABC_Marjoram_original

### MaC


# My version with error

MaC.fit <- MaC(targets.empirical = sum_stat_obs, # dummy.targets.empirical,
               RMSD.tol.max = 2,
               min.givetomice = 2, # 12, # 64,
               n.experiments = 2, # 122, #256,
               lls = c(0,1),
               uls = c(1,2),
               model = toy_model, #VEME.wrapper, # simpact.wrapper,
               strict.positive.params = c(1,2),
               probability.params = c(1),
               method = "norm",
               predictorMatrix = "complete",
               maxit = 50,
               maxwaves = 4,
               n_cluster = 4) # n_cluster)


# Version from supervisor

# me adding packages
library(data.table)
library(mice)
library(pcaPP)
library(randtoolbox)
library(readcsvcolumns)
library(gsubfn)
library(glmnet)

# BEGIN

library(RSimpactHelper)
# THIS IS ESSENTIAL: So that the imputation functions are in the global environment, and can be found.
mice.impute.norm <- mice::mice.impute.norm
mice.impute.rf <- mice::mice.impute.rf

toy_model <- function(x){
  set.seed(x[1])   # THIS IS ESSENTIAL
  c( x[2] + x[3] + rnorm(1,0,0.1) ,   # x[1] and x[2] now become x[2] and x[3]
     x[2] * x[3] + rnorm(1,0,0.1) )
}

sum_stat_obs <- c(1.5, 0.5)
lls = c(0,1)
uls = c(1,2)

help(MaC)

MaC.toy <- MaC(targets.empirical = sum_stat_obs,
               RMSD.tol.max = 2,
               min.givetomice = 200,
               n.experiments = 2000,
               lls = lls,
               uls = uls,
               model = toy_model,
               strict.positive.params = 0,
               probability.params = 0,
               method = "norm",
               predictorMatrix = "complete",
               maxit = 20,
               maxwaves = 6,
               n_cluster = 8)
# Plotting the input parameters of the calibrated model
plot(MaC.toy$selected.experiments[[6]][, 1],
     MaC.toy$selected.experiments[[6]][, 2])
# Plotting the summary statistics of the calibrated model
plot(MaC.toy$selected.experiments[[6]][, 3],
     MaC.toy$selected.experiments[[6]][, 4])

## END

targets.stat <- read.csv("~/Desktop/mastermodeltest/features.matrix.csv")

median.targets.stat <-  colMedians(as.matrix(targets.stat))  # Ok # library(robustbase)

classic.target <- colMedians(as.matrix(targets.stat[,12:27]))

dummy.targets.empirical <- as.numeric(classic.target)

MaC.fit <- MaC(targets.empirical = dummy.targets.empirical, # dummy.targets.empirical,
               RMSD.tol.max = 2,
               min.givetomice =64,
               n.experiments = 256,
               lls = c(0.5, 0, -0.1, -0.5, 2, 0, -1, -0.9, 0, -2, -100, 0, 0, -0.5, 3, 5, 10, -3.5),
               uls = c(2, 0, 3, 0.5, 4, 1, 0, 0, 1, 0, -80, 1, 1, 1, 7, 9, 14, -1.7),
               model = wrapper.test.study.1, #VEME.wrapper, # simpact.wrapper,
               strict.positive.params = c(3, 5, 6, 9, 12, 13, 15, 16, 17),
               probability.params = c(6, 9, 12, 13),
               method = "norm",
               predictorMatrix = "complete",
               maxit = 50,
               maxwaves = 4,
               n_cluster = 4) # n_cluster)



