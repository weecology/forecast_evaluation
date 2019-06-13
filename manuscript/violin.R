
# type options
# l, b, p, o, h, r (new is "rectangle histogram")


source("violin_functions.R")

set.seed(123)
xx <- rnorm(1e5, 0, 1)
boxplot(xx, xlim = c(-1, 2), ylim = c(-5, 5))
violin(xx, 1)
violin(round(xx), 1, col = rgb(0, 0.4, 0.9), wex = 10)


xx <- rpois(1e5, 1)
boxplot(xx, xlim = c(-1, 2), ylim = c(-5, 7))
violin(xx, 1, col = rgb(0, 0.4, 0.9), lty = 3, lwd = 4)


class(density(seq(1:20)))