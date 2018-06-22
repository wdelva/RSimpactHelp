wave1.creator <- function(file.name = "/user/data/gent/vsc400/vsc40070/agemixing/scripts/wave1.csv",
                          n.experiments,
                          limits){ # limits <- list(prior.x.1 = c(0,1), prior.x.2 = c(2,4), prior.x.3 = c(-1,1))
  lls <- as.numeric(as.data.frame(limits)[1, ])
  uls <- as.numeric(as.data.frame(limits)[2, ])
  range.width <- uls - lls
  ll.mat <- matrix(rep(lls, n.experiments), nrow = n.experiments, byrow = TRUE)
  range.width.mat <- matrix(rep(range.width, n.experiments), nrow = n.experiments, byrow = TRUE)
  sobol.seq.0.1 <- sobol(n = n.experiments, dim = length(lls), init = TRUE, scrambling = 1, seed = 1, normal = FALSE)
  experiments <- cbind(1:n.experiments, ll.mat + sobol.seq.0.1 * range.width.mat)
  experiments.col.names <- c("seed", paste0("x", seq(1:length(lls)),"stop"))
  write.table(experiments,
              sep = ",",
              dec = ".",
              quote = FALSE,
            row.names = FALSE,
            col.names = experiments.col.names,
            file = file.name)
}
