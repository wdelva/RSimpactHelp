library(magrittr)
library(dplyr)
library(emulator)
library(multivator)
## Simple Simpact emulator

simpact.inANDout.df <- read.csv(file = "/Users/delvaw/Documents/RSimpactHelper/data/inANDout.df2016-09-27.csv",
         header = TRUE,
         sep = ",")

# The x variables (model parameters) that were varied:
# concept.base
# eagerness.a
# eagerness.b
# formation.hazard.agegapry.numrel_man = formation.hazard.agegapry.numrel_woman
# formation.hazard.agegapry.eagerness_sum
# formation.hazard.agegapry.eagerness_diff
simpact.x <- dplyr::select(simpact.inANDout.df,
              concept.base,
              eagerness.a,
              eagerness.b,
              formation.hazard.agegapry.numrel_man,
              formation.hazard.agegapry.eagerness_sum,
              formation.hazard.agegapry.eagerness_diff) %>% as.matrix()

# The "z" variables (model summary statistics), for which we will define targets:
# growth.rate
# prev.15.50.end
simpact.z <- dplyr::select(simpact.inANDout.df,
                           growth.rate,
                           prev.15.50.end)

# The simulations ran for 40 years, with a starting population of 1000 people.


names(simpact.inANDout.df)
hist(simpact.inANDout.df$concept.base)
hist(simpact.inANDout.df$formation.hazard.agegapry.eagerness_sum)


#library(lhs)

# Defining the target statistics for which we will seek the most likely parameters.

# z1 (growth.rate) ~ 0
# z2 (prev.15.50.end) ~ 0.1


targets <- c(0, 0.1)

#### Setting up design. Creating the LHS over the 0-1 uniform parameter space
design.points <- 1000 #40 #20
variables <- 6

# On VSC, when each core or node gets its own series of simulations to run, the seed must be set to a different value for each core or node.
# set.seed(1)
# rlhs <- randomLHS(design.points, variables) # We should check if objects created with randomLHS also work. They should work.
x1.design <- punif(simpact.x[, 1],
                   min = round(min(simpact.x[, 1])),
                   max = round(max(simpact.x[, 1])))
x2.design <- punif(simpact.x[, 2],
                   min = round(min(simpact.x[, 2])),
                   max = round(max(simpact.x[, 2])))
x3.design <- punif(simpact.x[, 3],
                   min = round(min(simpact.x[, 3])),
                   max = round(max(simpact.x[, 3])))
x4.design <- punif(simpact.x[, 4],
                   min = round(min(simpact.x[, 4])),
                   max = round(max(simpact.x[, 4])))
x5.design <- punif(simpact.x[, 5],
                   min = round(min(simpact.x[, 5])),
                   max = round(max(simpact.x[, 5])))
x6.design <- punif(simpact.x[, 6],
                   min = round(min(simpact.x[, 6])),
                   max = round(max(simpact.x[, 6])))

# x.design <- latin.hypercube(design.points, variables, names=c("x1","x2", "x3", "x4", "x5"))
x.design <- cbind(x1.design, x2.design, x3.design, x4.design, x5.design, x6.design)
x.design.long <- rbind(x.design, x.design) # one for each z component

#### Computing model output at design points
z_obs <- as.vector(unlist(simpact.z))

#### Before we start the emulation, let's see what the best fit is from the model simulation runs
z.df <- simpact.z
hist(z.df$growth.rate)
hist(z.df$prev.15.50.end)
#names(z.df) <- c("z1", "z2", "z3", "z4")
sq.array <- as.data.frame(t(apply(z.df, 1, function(x) (x - t(c(0, 0.1)))^2)))
SumSq <- as.numeric(rowSums(sq.array))
which.min(SumSq)
z.df[which.min(SumSq), ] #                   growth.rate prev.15.50.end
                                          #  -0.02872134     0.07633588
# And the targets:                            0              0.10

#### And most importantly, the best estimate for the model parameters:
                                                         #  concept.base eagerness.a eagerness.b
x.estimate <- as.numeric(simpact.x[which.min(SumSq), ])  #  -1.56448040  0.70100469  68.50956153


                                                        #  formation.hazard.agegapry.numrel_man  formation.hazard.agegapry.eagerness_sum formation.hazard.agegapry.eagerness_diff
                                                        #   -0.01584842                          -1.87933009                                  -0.45373292


names(simpact.x)


#### Creating the multivator objects
RS_mdm <- mdm(x.design.long, types = rep(c("growth.rate", "prev.15.50.end"), each = design.points))
RS_expt <- experiment(mm = RS_mdm, obs = z_obs)
#RS_opt <- optimal_params(RS_expt, option="a")
RS_opt_c <- optimal_params(RS_expt, option="c")

#### Using the emulator to explore the parameter space
n <- 5000
x.new <- latin.hypercube(n, variables, names=colnames(x.design))
RS_new_mdm <- mdm(rbind(x.new, x.new), types = rep(c("growth.rate", "prev.15.50.end"), each = n))
RS_prediction_c <- multem(x = RS_new_mdm, expt = RS_expt, hp = RS_opt_c)

hist(RS_prediction_c[1:n]) # distribution of z1
hist(RS_prediction_c[(n+1):(2*n)]) # distribution of z2


#### One way of efficiently comparing emulation output with target statistics is to reshape RS_prediction as a dataframe
prediction.df <- data.frame(matrix(RS_prediction_c, nrow = n, dimnames = list(rownames = 1:n, colnames = c("growth.rate", "prev.15.50.end"))))

#### sum of squared distances between model statistics and target statistics
# Note: we could normalise (and centralise) these statistics to give them more equal weight in the SumSq
sq.array <- as.data.frame(t(apply(prediction.df, 1, function(x) (x - t(c(0, 0.1)))^2)))
names(sq.array) <- names(prediction.df)
SumSq <- as.numeric(rowSums(sq.array))
which.min(SumSq)
prediction.df[which.min(SumSq), ] #            growth.rate prev.15.50.end
#                                              -0.02847503     0.03123356
# very bad, in light of targets:                0              0.1

#### And most importantly, the best estimate for the model parameters:
x.estimate <- as.numeric(simpact.x[which.min(SumSq), ])

# The x variables (model parameters) that were varied:
# concept.base -4.8868513
# eagerness.a 0.6228266
# eagerness.b 144.5345675
# formation.hazard.agegapry.numrel_man = formation.hazard.agegapry.numrel_woman -0.6530806
# formation.hazard.agegapry.eagerness_sum -1.4261494
# formation.hazard.agegapry.eagerness_diff -0.2541594




#### From here, we can improve on this method in at least 3 non-mutually exclusive ways:
# 1. We run more initial model simulations (more design.points) so that the emulator has more data to work with (it worked with 20 and 40)
# 2. Try different argument values for optimal_params()
# 3. We expore more narrowly around the parameter estimates obtained from the first step with a second step of simulation+emulation
# 4. We apply Principal Component Analysis to reduce the dimensions of the model output (see McNeal example page 16 of multivator tutorial)

# Principal component analysis
pairs(z.df)

z.pc <- princomp(z.df, scores = TRUE, cor = TRUE)
summary(z.pc)
plot(z.pc)
biplot(z.pc)
z.pc$loadings
z.pc$scores
targets.df = data.frame(growth.rate = targets[1], prev.15.50.end = targets[2])#, z3 = targets[3], z4 = targets[4])
# As an example: the value of the first PC for the target statistics:
as.numeric(as.numeric(z.pc$loadings[, 1]) %*% ((as.numeric(targets.df) - z.pc$center) / z.pc$scale) )
# All PCs for the target statistics:
targets.pc <- predict(z.pc, targets.df)[1:2]
# So the proposal is to only use the first 3 components
# because they jointly capture 96% of the total (normalised, standardised) variance

# As an example: recovering the value of the first statistics (z1), based on the loadings and scores
as.numeric(as.numeric(z.pc$loadings[1, ]) %*% as.numeric(z.pc$scores[1, ])) * as.numeric(z.pc$scale[1]) + as.numeric(z.pc$center[1])
# And the second statistic:
as.numeric(as.numeric(z.pc$loadings[2, ]) %*% as.numeric(z.pc$scores[1, ])) * as.numeric(z.pc$scale[2]) + as.numeric(z.pc$center[2])
z.df[1, ]

# If we, however, discard the third PC, because it captures less than 5% of the total variance,
# we can only recover statistics to some approximate level of accuracy:
as.numeric(as.numeric(z.pc$loadings[1, 1:3]) %*% as.numeric(z.pc$scores[1, 1:3])) * as.numeric(z.pc$scale[1]) + as.numeric(z.pc$center[1])
# And the second statistic:
as.numeric(as.numeric(z.pc$loadings[2, 1:3]) %*% as.numeric(z.pc$scores[1, 1:3])) * as.numeric(z.pc$scale[2]) + as.numeric(z.pc$center[2])
z.df[1, ]

# Let's build the emulator now, using only the first 3 PCs.
z.pc$scores

#### Creating the multivator objects for the PC-based emulator
x.design.long.pc <- rbind(x.design, x.design, x.design) # one for each of the 3 PC components
pc_obs <- as.vector(z.pc$scores)[1:(3*design.points)]

RS_mdm.pc <- mdm(x.design.long.pc, types = rep(c("pc1", "pc2", "pc3"), each = design.points))
RS_expt.pc <- experiment(mm = RS_mdm.pc, obs = pc_obs)
RS_opt.pc <- optimal_params(RS_expt.pc, option="b")

#### Using the emulator to explore the parameter space
n <- 5000
x.new <- latin.hypercube(n, variables, names=c("x1","x2", "x3", "x4", "x5"))
RS_new_mdm.pc <- mdm(rbind(x.new, x.new, x.new), types = rep(c("pc1", "pc2", "pc3"), each = n))
RS_prediction.pc <- multem(x = RS_new_mdm.pc, expt = RS_expt.pc, hp = RS_opt.pc)
hist(RS_prediction.pc[1:n]) # distribution of pc1
hist(RS_prediction.pc[(n+1):(2*n)]) # distribution of pc2

#### The targets, transformed to the first 3 PCs:
targets.pc   # 1.8825113 -0.4188274 -0.2740838


#### One way of efficiently comparing emulation output with target statistics is to reshape RS_prediction as a dataframe
prediction.pc.df <- data.frame(matrix(RS_prediction.pc, nrow = n, dimnames = list(rownames = 1:n, colnames = c("pc1", "pc2", "pc3"))))

#### sum of squared distances between model statistics and target statistics
# Note: we could normalise (and centralise) these statistics to give them more equal weight in the SumSq
sq.pc.array <- as.data.frame(t(apply(prediction.pc.df, 1, function(x) (x - t(targets.pc))^2)))
names(sq.pc.array) <- names(prediction.pc.df)
SumSq.pc <- as.numeric(rowSums(sq.pc.array))
which.min(SumSq.pc)
prediction.pc.df[which.min(SumSq.pc), ] #            pc1       pc2       pc3
#                                                1.951745 -0.4689429 -0.2492109
# Not bad, in light of targets:                  1.882511 -0.4188274 -0.2740838

#### And most importantly, the best estimate for the model parameters:
x.standardised.pc <- as.numeric(x.new[which.min(SumSq.pc), ])
x1 <- qunif(x.standardised.pc[1], min = x1.min.max[1], max = x1.min.max[2])
x2 <- qunif(x.standardised.pc[2], min = x2.min.max[1], max = x2.min.max[2])
x3 <- qunif(x.standardised.pc[3], min = x3.min.max[1], max = x3.min.max[2])
x4 <- qunif(x.standardised.pc[4], min = x4.min.max[1], max = x4.min.max[2])
x5 <- qunif(x.standardised.pc[5], min = x5.min.max[1], max = x5.min.max[2])
x.estimate <- c(x1, x2, x3, x4, x5)
#                       -0.08530    1.38561   22.93771   -1.96600  1861
# But the truth was:     0          1.5       10         -1.5      1000

#### The average model output for the best estimate of the model parameters:
reps <- 10000
df <- data.frame(matrix(rep(NA, time = reps, each = 4), ncol = 4))
for(rep in 1:reps){
  df[rep, ] <- truth(x = x.estimate)
}

colMeans(df)
  #old values before PCA        # -0.1894938  2.3169974 31.3353000  4.8122847

                                # -0.08514   1.38520   162.373      5.510
# Much closer to targets:          0         1.5       158.7        5.5














# # Let's highlight the x-y coordinates that fit (produce z values within some target range)
# a.fit <- a > 0.7 & a < 0.9
# b.fit <- b > 4 & b < 6
# ab.fit <- a.fit & b.fit # 6 coordinates give a good fit
# x.fit <- x_a[ab.fit]
# y.fit <- y_a[ab.fit]
# points(x.fit, y.fit, col = "white", pch = 16)
#
# # Lastly, let's look at the "true" a and b for the x and y values that were identified by the emulator as good fits.
# xy.calibrated <- cbind(x.fit, y.fit)
# a_validated <- apply(xy.calibrated,1,fa)
# b_validated <- apply(xy.calibrated,1,fb)
#
#
#
#
# # Creating the LHS dataframe
# lhs.df <- data.frame(x1 = rep(NA, design.points),
#                      x2 = rep(NA, design.points),
#                      x3 = rep(NA, design.points),
#                      x4 = rep(NA, design.points),
#                      x5 = rep(NA, design.points))
# lhs.df$x1 <- qunif(rlhs[ , 1], min = -5, max = -1)
# lhs.df$x2 <- qunif(rlhs[ , 2], min = 0.1, max = 6)
# lhs.df$x3 <- qunif(rlhs[ , 3], min = 0.01, max = 100)
# lhs.df$x4 <- qunif(rlhs[ , 4], min = -10, max = 10)
# lhs.df$x5 <- qunif(rlhs[ , 5], min = 10, max = 10000)
#
#
# # Creating a dataframe for input AND output
# inANDout.df <- cbind.data.frame(sim.id = 1:design.points,
#                                 lhs.df,
#                                 z1 = rep(NA, design.points),
#                                 z2 = rep(NA, design.points),
#                                 z3 = rep(NA, design.points),
#                                 z4 = rep(NA, design.points))
#
# model <- function(X){
#   x1 <- X[1]; x2 <- X[2]; x3 <- X[3]; x4 <- X[4]; x5 <- X[5];
#   y1 <- rnorm(n = x5, mean = x1, sd = x2)
#   y2 <- y1[y1 < x4]
#
#   z1 <- mean(y1)
#   z2 <- sd(y1)
#   z3 <- length(y2)
#   z4 <- log(z2 * z3 + x3)
#   output <- c(x1, x2, x3, x4, x5, z1, z2, z3, z4)
#   return(output)
# }
#
# simpact4emulation <- function(sim.id, lhs.df){
#   X <- as.numeric(lhs.df[sim.id, ])
#   output <- model(X)
#   return(output)
# }
#
# # First wave of naive simulations
# for (sim.id in inANDout.df$sim.id){
#   inANDout.df[sim.id, 2:10] <- simpact4emulation(sim.id, lhs.df)
# }
#
# targets <- c(0, 1.5, 158.7, 5.5)
#
# distance <- (inANDout.df[, 7] - targets[1])^2 +
#   (inANDout.df[, 8] - targets[2])^2 +
#   (inANDout.df[, 8] - targets[3])^2 +
#   (inANDout.df[, 10] - targets[4])^2
#
# # Looking for a clever and efficient emulation scheme
# # We start with an emulator for x1
# emul.x1 <- lm(x1 ~ z1*z2*z3*z4, data = inANDout.df)
# newdata0 <- data.frame(z1 = targets[1],
#                        z2 = targets[2],
#                        z3 = targets[3],
#                        z4 = targets[4])
# p1 <- as.numeric(predict(emul.x1, newdata = newdata0, interval = "prediction"))
# p1vect <- seq(from = p1[2], to = p1[3], length.out = design.points)
#
# emul.x2 <- lm(x2 ~ z1*z2*z3*z4*x1, data = inANDout.df)
# newdata1 <- data.frame(z1 = targets[1],
#                        z2 = targets[2],
#                        z3 = targets[3],
#                        z4 = targets[4],
#                        x1 = p1vect)
# p2vect <- as.numeric(predict(emul.x2, newdata = newdata1))
#
# emul.x3 <- lm(x3 ~ z1*z2*z3*z4*x1*x2, data = inANDout.df)
# newdata2 <- data.frame(z1 = targets[1],
#                        z2 = targets[2],
#                        z3 = targets[3],
#                        z4 = targets[4],
#                        x1 = p1vect,
#                        x2 = p2vect)
# p3vect <- as.numeric(predict(emul.x3, newdata = newdata2))
#
# emul.x4 <- lm(x4 ~ z1*z2*z3*z4*x1*x2*x3, data = inANDout.df)
# newdata3 <- data.frame(z1 = targets[1],
#                        z2 = targets[2],
#                        z3 = targets[3],
#                        z4 = targets[4],
#                        x1 = p1vect,
#                        x2 = p2vect,
#                        x3 = p3vect)
# p4vect <- as.numeric(predict(emul.x4, newdata = newdata3))
#
# emul.x5 <- lm(x5 ~ z1*z2*z3*z4*x1*x2*x3*x4, data = inANDout.df)
# newdata4 <- data.frame(z1 = targets[1],
#                        z2 = targets[2],
#                        z3 = targets[3],
#                        z4 = targets[4],
#                        x1 = p1vect,
#                        x2 = p2vect,
#                        x3 = p3vect,
#                        x4 = p4vect)
# p5vect <- as.numeric(predict(emul.x5, newdata = newdata4))
#
#
# corr.df <- data.frame(x1 = p1vect,
#                       x2 = p2vect,
#                       x3 = p3vect,
#                       x4 = p4vect,
#                       x5 = round(p5vect))
#
# # Creating a dataframe for input AND output
# inANDout2.df <- cbind.data.frame(sim.id = 1:design.points,
#                                 corr.df,
#                                 z1 = rep(NA, design.points),
#                                 z2 = rep(NA, design.points),
#                                 z3 = rep(NA, design.points),
#                                 z4 = rep(NA, design.points))
#
# inANDout2.df <- filter(inANDout2.df, x5>0)
#
# # Second wave of simulations
# for (sim.id in inANDout2.df$sim.id){
#   inANDout2.df[sim.id, 2:10] <- simpact4emulation(sim.id, corr.df)
# }
#
#
# distance <- (inANDout2.df[, 7] - targets[1])^2 +
#   (inANDout2.df[, 8] - targets[2])^2 +
#   (inANDout2.df[, 8] - targets[3])^2 +
#   (inANDout2.df[, 10] - targets[4])^2
#
# hist(distance)
#
#
#
#
#
#
#
#
#
# p1 <- predict(emul.x1,
#               newdata = newdata0,
#               interval = "prediction")
# newdata1 <- data.frame(z1 = targets[1],
#                        z2 = targets[2],
#                        z3 = targets[3],
#                        z4 = targets[4],
#                        x1 = seq(from = p1[2], to = p1[3], length.out = design.points))
#
#
#
#
#
# emul.eagerness.a <- glm(eagerness.a ~ growth.rate + prev.15.50.end,
#                         family = poisson(link = "log"),
#                         data = inANDout.df)
# emul.eagerness.b <- glm(eagerness.b ~ growth.rate + prev.15.50.end,
#                         family = poisson(link = "log"),
#                         data = inANDout.df)
#
#
#
#
#
# na <- 33
# nb <- 9
# xa <- latin.hypercube(na,2)
# xb <- xa[seq_len(nb),]
#
# fa <- function(xy){
#   sin(5*(xy[1]+xy[2]))
# }
#
# fb <- function(xy){
#   7 * sin(5*(xy[1]+xy[2])) + sin(20*(xy[1]-xy[2]))
# }
#
# a_obs <- apply(xa,1,fa)
# b_obs <- apply(xb,1,fb)
#
# RS_mdm <- mdm(rbind(xa,xb),types=c(rep("a",na),rep("b",nb)))
# RS_expt <- experiment(mm=RS_mdm, obs= c(a_obs,b_obs))
# RS_opt <- optimal_params(RS_expt, option="b")
#
