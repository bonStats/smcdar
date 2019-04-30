# lotka-volterra models
library(markovchain)
library(smcdar)

# simulate
y1_max <- 100
y2_max <- 100
true_theta <- runif(3)

G <- generator_matrix_lotka_volterra_ctmc(theta = true_theta, y1_max = y1_max, y2_max = y2_max)

lotka_volt_ctmc <- new("ctmc", states = 1:nrow(G), generator = G)
