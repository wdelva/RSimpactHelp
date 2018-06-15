library(devtools)
# install_github("wdelva/RSimpactHelp")
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
lls = c(0, 0)
uls = c(1.5, 2)

# To get more info about the arguments of the MaC function:
# help(MaC)

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


# Let's compare this to accept-reject ABC
library(EasyABC)
# help("ABC_rejection")
toy_prior <- list(c("unif", lls[1], uls[1]),
                  c("unif", lls[2], uls[2]))
Rej.toy <- ABC_rejection(model = toy_model,
                         prior = toy_prior,
                         summary_stat_target = sum_stat_obs,
                         nb_simul = 12000,
                         use_seed = TRUE,
                         seed_count = 1,
                         n_cluster = 8,
                         tol = 200/12000)

# Let's compare this to sequential ABC
# help("ABC_sequential")
Seq.toy <- ABC_sequential(model = toy_model,
                          method = "Lenormand",
                          prior = toy_prior,
                          summary_stat_target = sum_stat_obs,
                          nb_simul = 2000,
                          alpha = 0.1,
                          p_acc_min = 0.03,
                          use_seed = TRUE,
                          seed_count = 1,
                          n_cluster = 8,
                          inside_prior = FALSE)
# To see how many waves were done:
1 + (Seq.toy$nsim - 2000) / 1800

# Plotting the input parameters of the calibrated model
plot(Rej.toy$param[, 1],
     Rej.toy$param[, 2],
     pch = 16,
     col = "black",
     xlab = "parameter 1",
     ylab = "parameter 2",
     xlim = c(0,1.5),
     ylim = c(0, 2))
points(Seq.toy$param[, 1],
       Seq.toy$param[, 2],
       pch = 16,
       col = "blue2")
points(MaC.toy$selected.experiments[[6]][, 1],
       MaC.toy$selected.experiments[[6]][, 2],
       pch = 16,
       col = "orange")
# Plotting the summary statistics of the calibrated model
plot(Rej.toy$stats[, 1],
     Rej.toy$stats[, 2],
     pch = 16,
     col = "black",
     xlab = "summary statistic 1",
     ylab = "summary statistic 2",
     xlim = c(1.4, 1.6),
     ylim = c(0.3, 0.7))
points(Seq.toy$stats[, 1],
       Seq.toy$stats[, 2],
       pch = 16,
       col = "blue2")
points(MaC.toy$selected.experiments[[6]][, 3],
       MaC.toy$selected.experiments[[6]][, 4],
       pch = 16,
       col = "orange")

points(sum_stat_obs[1],
       sum_stat_obs[2],
       pch = "+",
       cex = 3,
       col = "red3")

# Comparing simulation time:
MaC.toy$secondspassed[3]
Rej.toy$computime
Seq.toy$computime
