# lotka-volterra models
library(markovchain)
library(smcdar)

# simulate
y1_max <- 100
y2_max <- 100
true_theta <- runif(3)

G <- generator_matrix_lotka_volterra_ctmc(theta = true_theta, y1_max = y1_max, y2_max = y2_max)

lotka_volt_ctmc <- new("ctmc", states = as.character(1:nrow(G)), generator = as.matrix(G))

initial <- rep(0, times = 10202)
initial[6000] <- 1

rd <- rctmc(n=Inf, lotka_volt_ctmc, T = 7)

system.time({
  replicate(n = 100,
            rctmc(n=Inf, lotka_volt_ctmc, T = 7)
  )
})

system.time({
  replicate(n = 100,
  simulate_ctmc(generator_matrix = G, total_time = 7, initial_state  = 6000)
  )
})

