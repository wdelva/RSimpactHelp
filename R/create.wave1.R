#lls = c(0.1, -3, 0.1)
#uls = c(0.9, 0, 5.0)
library(devtools)
#install_github("wdelva/RSimpactHelp")

library(RSimpactCyan)
library(RSimpactHelper)
library(randtoolbox)
source(file = "/user/data/gent/vsc400/vsc40070/agemixing/scripts/wave1.creator.R") 
wave1.creator(file.name = "/user/data/gent/vsc400/vsc40070/agemixing/scripts/wave1.csv",
              n.experiments = 48,
              limits <- list(prior.x.1 = c(0.1, 0.9),
                             prior.x.2 = c(-3, 0),
                             prior.x.3 = c(0.1, 5)))
