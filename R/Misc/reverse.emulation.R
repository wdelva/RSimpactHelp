# Simple test case for model fitting with "reverse emulation"

# Calling the libraries we will need
library(RSimpactCyan)
library(RSimpactHelper)
library(data.table)
library(dplyr)
library(magrittr)
library(exactci)
library(nlme)
library(ggplot2)
install.packages("readcsvcolumns", repos="http://193.190.10.42/jori/")
library(readcsvcolumns)
library(survival)
library(KMsurv)
library(tidyr)
library(lhs)

# Is it a standard or a MaxART simulation?
simpact.set.simulation("simpact-cyan")#("maxart")

agedist.data.frame <- agedistr.creator(shape = 5, scale = 65)

testinput <- list()
testinput$mortality.normal.weibull.shape <- 5
testinput$mortality.normal.weibull.scale <- 65
testinput$mortality.normal.weibull.genderdiff <- 0
testinput$periodiclogging.interval <- 1
testinput$syncrefyear.interval <- 1
testinput$formation.hazard.type <- "agegapry"
testinput$person.eagerness.dist.type <- "gamma"
testinput$person.eagerness.dist.type <- "gamma"
testinput$person.eagerness.dist.gamma.a <- 0.231989836885#0.15 #0.425#3.4#1.7#0.85 #0.1 ~ kappa
testinput$person.eagerness.dist.gamma.b <- 45#70#100 #3.5#5#10#20 #170 ~ theta
testinput$person.agegap.man.dist.type <- "normal"
testinput$person.agegap.woman.dist.type <- "normal"
testinput$person.agegap.man.dist.normal.mu <- -5
testinput$person.agegap.woman.dist.normal.mu <- 2.5
testinput$person.agegap.man.dist.normal.sigma <- 1
testinput$person.agegap.woman.dist.normal.sigma <- 1
testinput$formation.hazard.agegapry.numrel_man <- -0.5
testinput$formation.hazard.agegapry.numrel_woman <- -0.5
testinput$formation.hazard.agegapry.gap_factor_man_exp <- -0.15# -0.5
testinput$formation.hazard.agegapry.gap_factor_woman_exp <- -0.15# -0.5
testinput$formation.hazard.agegapry.gap_factor_man_const <- 0
testinput$formation.hazard.agegapry.gap_factor_woman_const <- 0
testinput$formation.hazard.agegapry.gap_agescale_man <- 0.23
testinput$formation.hazard.agegapry.gap_agescale_woman <- 0.23
testinput$formation.hazard.agegapry.eagerness_sum <- 0.1
testinput$formation.hazard.agegapry.eagerness_diff <- -0.048#-0.110975
testinput$dissolution.alpha_0 <- -0.52#-0.1 # baseline
testinput$dissolution.alpha_4 <- -0.05
testinput$debut.debutage <- 14
testinput$population.simtime <- 40
testinput$population.nummen <- 500
testinput$population.numwomen <- 500
testinput$population.maxevents <- testinput$population.simtime * testinput$population.nummen * 4 # If 4 events happen per person per year, something's wrong.
testinput$population.eyecap.fraction <- 0.2
testinput$hivseed.type <- "amount"
testinput$hivseed.amount <- 20
testinput$hivseed.age.min <- 15
testinput$hivseed.age.max <- 25
testinput$hivseed.time <- 10
testinput$transmission.param.a <- -1.0352239
testinput$transmission.param.b <- -89.339994
testinput$transmission.param.c <- 0.4948478
testinput$transmission.param.f1 <- log(5) # ~1.6 such that the hazard is x 5 in 15 yo
testinput$transmission.param.f2 <- log(log(2.5) / log(5)) / 5
testinput$conception.alpha_base <- -3


# Creating the LHS over the 0-1 uniform parameter space
design.points <- 100
variables <- 3
# On VSC, when each core or node gets its own series of simulations to run, the seed must be set to a different value for each core or node.
set.seed(1)
rlhs <- randomLHS(design.points, variables)

# Creating the LHS dataframe
lhs.df <- data.frame(concept.base = rep(NA, design.points), # tunes conception rate
                     eagerness.a = rep(NA, design.points), # tunes relationship rate and distribution
                     eagerness.b = rep(NA, design.points)) # tunes relationship rate and distribution
lhs.df$concept.base <- qunif(rlhs[ , 1], min = -5, max = -1)
lhs.df$eagerness.a <- qunif(rlhs[ , 2], min = 0.1, max = 1)
lhs.df$eagerness.b <- qunif(rlhs[ , 3], min = 10, max = 200)

# Creating a dataframe for input AND output
inANDout.df <- cbind.data.frame(sim.id = 1:design.points,
                                lhs.df,
                                growth.rate = rep(NA, design.points),
                                prev.15.50.end = rep(NA, design.points))

# Creating a new Simpact4emulation function
simpact4emulation <- function(sim.id, lhs.df, testinput, agedist.data.frame){

  testinput$conception.alpha_base <- lhs.df$concept.base[sim.id]
  testinput$person.eagerness.dist.gamma.a <- lhs.df$eagerness.a[sim.id]
  testinput$person.eagerness.dist.gamma.b <- lhs.df$eagerness.b[sim.id]
  simpact.seed.id <- sim.id

  testoutput <- simpact.run(configParams = testinput,
                            destDir = "temp",
                            agedist = agedist.data.frame,
                            seed = simpact.seed.id)

  if (testoutput$simulationtime < testinput$population.simtime)
  {
    # Ik kan op dit moment twee redenen bedenken waarom de simulatie te vroeg zou stoppen
    #  - maximaal aantal events is bereikt, kan op gecheckt worden dmv ret["eventsexecuted"]
    #  - geen events meer, gebeurt bvb als populatie uitsterft.
    if (testoutput$eventsexecuted >= testinput$population.maxevents-1)
    {
      # Ik doe hier een -1 omdat in R getallen standaard voorgesteld worden als floating
      # point getallen, en echte gelijkheden daarmee nogal gevaarlijk zijn.
      stop("MAXEVENTS: Simulation stopped prematurely, max events reached")
    }
    else
    {
      # Misschien moet er in dit geval niet gestopt worden
      stop("Simulation stopped prematurely, probably ran out of events")
    }
  }


  datalist.test <- readthedata(testoutput)
  growth.rate <- pop.growth.calculator(datalist = datalist.test, timewindow = c(0, testinput$population.simtime))
  prev.15.50.end <- prevalence.calculator(datalist = datalist.test, agegroup = c(15, 50), timepoint = testinput$population.simtime)

  output.stats <- c(growth.rate, prev.15.50.end$pointprevalence[3])
  return(output.stats)
}


# Testing simpact4emulation
#out.test <- simpact4emulation(1, lhs.df, testinput, agedist.data.frame)
# This worked (before the IF statements were inserted), but growth.rate was NA, because simulation was stopped prematurely, due to population.maxevents being reached

# Creating an error function to catch the case when population.maxevents is reached before population.simtime is reached
errFunction <- function(e){
  if (length(grep("MAXEVENTS",e$message)) != 0)
    return(inANDout.test = rep(NA, 2))
  # Een andere foutmelding dan MAXEVENTS, zorg dat we toch stoppen
  stop(e)
}



# Running the simpact4emulation function with error catcher in a loop
for (sim.id in inANDout.df$sim.id){
  out.test <- tryCatch(simpact4emulation(sim.id, lhs.df, testinput, agedist.data.frame), error = errFunction)
  # Inserting the output to the inANDout dataframe
  inANDout.df[sim.id, c(5, 6)] <- out.test
}

# As a reference point, let's calculate what fraction of the simulations was close enough to the
# target values of the output statistics
target.growth.rate <- 0
target.prev.15.50.end <- 0.1
# Let's accept growth rates between -0.02 and 0.02
# Let's accept HIV prevalences between 0.08 and 0.12
tolerance <- 0.02
success <- filter(inANDout.df, growth.rate >= target.growth.rate - tolerance &
                    growth.rate <= target.growth.rate + tolerance &
                    prev.15.50.end >= target.prev.15.50.end - tolerance &
                    prev.15.50.end <= target.prev.15.50.end + tolerance)

success.rate <- nrow(success) / nrow(inANDout.df) # success.rate is defined as the fraction of simulations that produced output within the tolerance band


# Now that we have a dataframe with n input parameters (3 in this case) and k output statistics (2 in this case),
# We fit n regression models (one for each input parameter): input_i = linkfunction(output1, output2) + random error function

emul.concept.base <- lm(concept.base ~ growth.rate + prev.15.50.end, data = inANDout.df)
emul.eagerness.a <- lm(eagerness.a ~ growth.rate + prev.15.50.end, data = inANDout.df)
emul.eagerness.b <- lm(eagerness.b ~ growth.rate + prev.15.50.end, data = inANDout.df)

# Let's now look at the prediction interval for the parameters, evaluating the emulator at the
# target values of the output statistics
newdata <- data.frame(growth.rate = target.growth.rate,
                      prev.15.50.end =target.prev.15.50.end)

p1 <- predict(emul.concept.base,
        newdata = newdata,
        interval = "prediction")

p2 <- predict(emul.eagerness.a,
        newdata = newdata,
        interval = "prediction")

p3 <- predict(emul.eagerness.b,
        newdata = newdata,
        interval = "prediction")

# To determine which of these prediction intervals is the narrowest, in relative terms,
# let's compare the width to the expected value, and also to the width of the initial lhs.
rel.widths <- abs(c((p1[3] - p1[2]) / p1[1],
                    (p2[3] - p2[2]) / p2[1],
                    (p3[3] - p3[2]) / p3[1]))
# rel.widths are 0.3061137 2.7696663 3.1030063 so the first prediction interval is much narrower.

# The initial range for p1 was -5 to -1. Prediction: -3.76 to -2.76
# So the range shrank from 4 to 1

# The initial range for p2 was 0.1 to 1. Prediction: -0.11 to 0.66 (but we KNOW that p2 must be >0)
# So the range shrank from 0.9 to 0.8 (or adjusted: shrank to 0.66)

# The initial range for p3 was 10 to 200. Prediction: -38.65 to 178 (but we KNOW that p3 must be >0)
# So the range INCREASED from 190 to 217 (or adjusted: shrank to 178)

# Let's do a second wave of simulations, this time with updated the ranges for the input parameters:
# For simplicity, we set the new ranges equal to the prediction intervals.


# Creating the LHS over the 0-1 uniform parameter space
design.points <- 100
variables <- 3
# On VSC, when each core or node gets its own series of simulations to run, the seed must be set to a different value for each core or node.
set.seed(2) # Because we want a new wave of simulations
rlhs2 <- randomLHS(design.points, variables)

# Creating the LHS dataframe
lhs2.df <- data.frame(concept.base = rep(NA, design.points), # tunes conception rate
                     eagerness.a = rep(NA, design.points), # tunes relationship rate and distribution
                     eagerness.b = rep(NA, design.points)) # tunes relationship rate and distribution
lhs2.df$concept.base <- qunif(rlhs2[ , 1], min = p1[2], max = p1[3])
lhs2.df$eagerness.a <- qunif(rlhs2[ , 2], min = 0, max = p2[3])
lhs2.df$eagerness.b <- qunif(rlhs2[ , 3], min = 0, max = p3[3])

# Creating a dataframe for input AND output
inANDout2.df <- cbind.data.frame(sim.id = 1:design.points,
                                lhs2.df,
                                growth.rate = rep(NA, design.points),
                                prev.15.50.end = rep(NA, design.points))

# Running the simpact4emulation function with error catcher in a loop
for (sim.id in inANDout2.df$sim.id){
  out.test <- tryCatch(simpact4emulation(sim.id, lhs2.df, testinput, agedist.data.frame), error = errFunction)
  # Inserting the output to the inANDout dataframe
  inANDout2.df[sim.id, c(5, 6)] <- out.test
}

success2 <- filter(inANDout2.df, growth.rate >= target.growth.rate - tolerance &
                    growth.rate <= target.growth.rate + tolerance &
                    prev.15.50.end >= target.prev.15.50.end - tolerance &
                    prev.15.50.end <= target.prev.15.50.end + tolerance)

success2.rate <- nrow(success2) / nrow(inANDout2.df) # success.rate is defined as the fraction of simulations that produced output within the tolerance band


### Instead of doing a second wave with a completely independent sampling of all parameters,
### We could use the output of the first wave to model parameter dependencies.
# Build an emulator that also has p1 as an independent ("predictor") variable for p2
# But we only use that part of the dataset for which p1 is within the prediction range.
# Then we create a new dataset, with datapoints for p1 spanning its prediction range,
# And we look at the range of expected values for p2.
# This second emulator with p1 on the right-hand side can also be used to derive
# a functional dependency between p1 and p2.

inANDout.narrow.df <- filter(inANDout.df, concept.base >= p1[2] & concept.base <= p1[3]) # That leaves only 25 rows

emul.eagerness.a <- lm(eagerness.a ~ growth.rate + prev.15.50.end + concept.base, data = inANDout.narrow.df)#inANDout.df)
emul.eagerness.b <- lm(eagerness.b ~ growth.rate + prev.15.50.end + concept.base, data = inANDout.narrow.df)#inANDout.df)

newdata.narrow <- data.frame(growth.rate = target.growth.rate,
                      prev.15.50.end =target.prev.15.50.end,
                      concept.base = inANDout.narrow.df$concept.base)

p2.new <- predict(emul.eagerness.a,
              newdata = newdata.narrow,
              interval = "prediction")

plot(inANDout.narrow.df$concept.base, p2.new[, 1], type = "l", ylim = c(-0.5, 0.5))
lines(inANDout.narrow.df$concept.base, p2.new[, 2], lty = "dashed")
lines(inANDout.narrow.df$concept.base, p2.new[, 3], lty = "dashed")

# Ideas to take further:
# 1. sample p2 and p3 in the second wave, conditional on the value of p1,
# With the sampling function given by the emulator that includes p1 as independent variable.

# 2. Choose more natural link functions (e.g. Poisson) for eagerness.a and eagerness.b because they must be > 0

# Key Question: do we want to identify "ALL" parameter sets that give good output?
# Or is our goal to find ONE? If the latter, we only care about the predicted values of p2 and p3 ("fit")

# Creating a non-random, but correlated dataframe
corr.df <- data.frame(concept.base = rep(NA, design.points), # tunes conception rate
                      eagerness.a = rep(NA, design.points), # tunes relationship rate and distribution
                      eagerness.b = rep(NA, design.points)) # tunes relationship rate and distribution

# For p1, we still use a naively uniform, but narrower sampling window
corr.df$concept.base <- qunif(rlhs2[ , 1], min = p1[2], max = p1[3])

# For p2 and p3 we use the emulator output
newdata.p1 <- data.frame(growth.rate = target.growth.rate,
                             prev.15.50.end =target.prev.15.50.end,
                             concept.base = corr.df$concept.base)

corr.df$eagerness.a <- predict(emul.eagerness.a, newdata = newdata.p1)

newdata.p1.p2 <- data.frame(growth.rate = target.growth.rate,
                            prev.15.50.end =target.prev.15.50.end,
                            concept.base = corr.df$concept.base,
                            eagerness.a = corr.df$eagerness.a)

corr.df$eagerness.b <- predict(emul.eagerness.b, newdata = newdata.p1.p2)


# Creating a dataframe for input AND output
inANDout2B.df <- cbind.data.frame(sim.id = 1:design.points,
                                  corr.df,
                                  growth.rate = rep(NA, design.points),
                                  prev.15.50.end = rep(NA, design.points))

# Running the simpact4emulation function with error catcher in a loop
for (sim.id in inANDout2B.df$sim.id){
  out.test <- tryCatch(simpact4emulation(sim.id, corr.df, testinput, agedist.data.frame), error = errFunction)
  # Inserting the output to the inANDout dataframe
  inANDout2B.df[sim.id, c(5, 6)] <- out.test
}

success2B <- filter(inANDout2B.df, growth.rate >= target.growth.rate - tolerance &
                     growth.rate <= target.growth.rate + tolerance &
                     prev.15.50.end >= target.prev.15.50.end - tolerance &
                     prev.15.50.end <= target.prev.15.50.end + tolerance)

success2B.rate <- nrow(success2B) / nrow(inANDout2B.df) # success.rate is defined as the fraction of simulations that produced output within the tolerance band


# From 2% success rate initially, to 6% with the first, naive method, and 13% with wave 2B.

# Let's take the output of wave 2B, and add a third wave:
# As before, we first construct an emulator for concept.base (p1):
emul.concept.base <- lm(concept.base ~ growth.rate + prev.15.50.end, data = inANDout2B.df)
emul.eagerness.a <- lm(eagerness.a ~ growth.rate + prev.15.50.end, data = inANDout2B.df)
emul.eagerness.b <- lm(eagerness.b ~ growth.rate + prev.15.50.end, data = inANDout2B.df)

# Let's now look at the prediction interval for the parameters, evaluating the emulator at the
# target values of the output statistics
newdata <- data.frame(growth.rate = target.growth.rate,
                      prev.15.50.end = target.prev.15.50.end)

p1 <- predict(emul.concept.base,
              newdata = newdata,
              interval = "prediction")

p2 <- predict(emul.eagerness.a,
              newdata = newdata,
              interval = "prediction")

p3 <- predict(emul.eagerness.b,
              newdata = newdata,
              interval = "prediction")

# To determine which of these prediction intervals is the narrowest, in relative terms,
# let's compare the width to the expected value, and also to the width of the initial lhs.
rel.widths <- abs(c((p1[3] - p1[2]) / p1[1],
                    (p2[3] - p2[2]) / p2[1],
                    (p3[3] - p3[2]) / p3[1]))
# rel.widths are 0.0617698 0.3193454 0.1841983 so the first prediction interval is narrower, but the others have narrowed a lot, compared to the first emulation round.

# The range for p1 was -3.76 to -2.76. Prediction: -3.18 to -3.28
# So the range shrank from 1 to 0.1

# The range for p2 was 0.011 to 0.30. Prediction: 0.15 to 0.21
# So the range shrank from 0.29 to 0.06

# The range for p3 was 69.8 to 175.8. Prediction: 103.7 to 124.7
# So the range shrank from 106 to 21



inANDout2B.narrow.df <- filter(inANDout2B.df, concept.base >= p1[2] & concept.base <= p1[3]) # That leaves only 20 rows

emul.eagerness.a <- lm(eagerness.a ~ growth.rate + prev.15.50.end + concept.base, data = inANDout2B.narrow.df)#inANDout.df)

# An open question whether this third emulator should be trained, using only the range of eagerness.a that came out of:
# corr2.df$eagerness.a <- predict(emul.eagerness.a, newdata = newdata.p1)
# For now we keep it a bit broader and more naive.
emul.eagerness.b <- lm(eagerness.b ~ growth.rate + prev.15.50.end + concept.base, data = inANDout2B.narrow.df)#inANDout.df)
# If we go with this more restrictive approach:
emul.eagerness.b <- lm(eagerness.b ~ growth.rate + prev.15.50.end + concept.base, data = inANDout2B.narrow.df)#inANDout.df)


# Creating a non-random, but correlated dataframe for the second time (which will be used in simulation wave 3)
corr2.df <- data.frame(concept.base = rep(NA, design.points), # tunes conception rate
                      eagerness.a = rep(NA, design.points), # tunes relationship rate and distribution
                      eagerness.b = rep(NA, design.points)) # tunes relationship rate and distribution

design.points <- 100
variables <- 3
# On VSC, when each core or node gets its own series of simulations to run, the seed must be set to a different value for each core or node.
set.seed(3) # Because we want a new wave of simulations
rlhs3 <- randomLHS(design.points, variables)

# For p1, we still use a naively uniform, but narrower sampling window
corr2.df$concept.base <- qunif(rlhs3[ , 1], min = p1[2], max = p1[3])

# For p2 and p3 we use the emulator output
newdata.p1 <- data.frame(growth.rate = target.growth.rate,
                         prev.15.50.end =target.prev.15.50.end,
                         concept.base = corr2.df$concept.base)

corr2.df$eagerness.a <- predict(emul.eagerness.a, newdata = newdata.p1)

newdata.p1.p2 <- data.frame(growth.rate = target.growth.rate,
                            prev.15.50.end =target.prev.15.50.end,
                            concept.base = corr2.df$concept.base,
                            eagerness.a = corr2.df$eagerness.a)

corr2.df$eagerness.b <- predict(emul.eagerness.b, newdata = newdata.p1.p2)


# Creating a dataframe for input AND output
inANDout3.df <- cbind.data.frame(sim.id = 1:design.points,
                                  corr2.df,
                                  growth.rate = rep(NA, design.points),
                                  prev.15.50.end = rep(NA, design.points))

# Running the simpact4emulation function with error catcher in a loop
for (sim.id in inANDout3.df$sim.id){
  out.test <- tryCatch(simpact4emulation(sim.id, corr2.df, testinput, agedist.data.frame), error = errFunction)
  # Inserting the output to the inANDout dataframe
  inANDout3.df[sim.id, c(5, 6)] <- out.test
}

success3 <- filter(inANDout3.df, growth.rate >= target.growth.rate - tolerance &
                      growth.rate <= target.growth.rate + tolerance &
                      prev.15.50.end >= target.prev.15.50.end - tolerance &
                      prev.15.50.end <= target.prev.15.50.end + tolerance)

success3.rate <- nrow(success3) / nrow(inANDout3.df) # success.rate is defined as the fraction of simulations that produced output within the tolerance band

# Succes rate in wave 3 is 10%, which is lower than the 13% in wave 2.

# Maybe the emulators in preparation of the third wave should fix concept.base at its highest expectation value,
# and should still allow for a wider range of eagerness.a and eagerness.be

# OR
# We should restrict the dataset for the third emulator, as suggested on line 376 (An open question...)

