#' Plot age-mixing pattern.
#'
#' Visualizes the age-mixing pattern from a Simpact-simulated population. This
#' function takes as input a datalist produced by \code{\link{pattern.modeller}}
#' and outputs a scatterplot of the age of the individual in the population
#' versus the age of their partner. The plot contains a dashed line that
#' represents a scenario where people choose partners that are the same age as
#' themselves. Juxtaposed with this is a line that represents the
#' population-average predicted partner ages based upon the generalised linear
#' mixed effects model in \code{\link{pattern.modeller}}.
#'
#'
#' @param dl The datalist that is produced by \code{\link{pattern.modeller()}}
#'
#' @return a large ggplot object of the age-mixing pattern
#'
#' @examples
#' load(obs)
#' amp <- pattern.plotter(dl = pat)
#' amp
#'
#' @importFrom magrittr %>%
#' @import dplyr
#' @import ggplot2


pattern.plotter <- function(dl) {

  df <- dl[[1]]

  theme <- theme(axis.title.x = element_text(face = "bold", size = 12),
                 axis.text.x  = element_text(angle = 90,
                                             face = "bold",
                                             colour = "black"),
                 axis.title.y = element_text(face = "bold", size = 12),
                 axis.text.y = element_text(angle = 90,
                                            face = "bold",
                                            colour = "black"),
                 plot.title = element_text(lineheight = .8, face = "bold"),
                 panel.grid.major = element_line(colour = "black"),
                 panel.grid.minor = element_line(colour = NA),
                 panel.background = element_rect(fill = "white"),
                 strip.background = element_rect(fill = "white"),
                 legend.title = element_text(colour = "black",
                                             size = 10,
                                             face = "bold"),
                 legend.text = element_text(size = 8, face = "bold"),
                 plot.margin = unit(c(1.5, 1, 1, 1), "lines"))

  fig <- df %>%
    ggplot(aes(x = agerelform, y = pred)) +
    geom_point(aes(x = agerelform, y = pagerelform),
               position = position_jitter(width = 0.75, height = 0.75),
               alpha = 0.5) +
    geom_abline(size = 1,
                aes(intercept = 0, slope = 1, linetype = "Same age"),
                show.legend = FALSE) +
    geom_line(aes(linetype = "Population mean"),
              size = 1) +
    facet_grid(. ~ Gender) +
    scale_y_continuous(name = "Partner's ages") +
    scale_linetype_manual('Lines',
                          values = c("Population mean" = 1, "Same age" = 2)) +
    xlab("Individual in population") +
    guides(linetype = guide_legend(keywidth = 2, keyheight = 1)) +
    theme

  return(fig)
}
