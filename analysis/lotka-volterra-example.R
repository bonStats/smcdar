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

theta <- c(0.5,0.1,0.9)

sol <- deSolve::ode(y = c(10,9,0,0,0),
             times = times,
             func = d_dt_lotka_volterra_lna,
             parms = theta,
             method = "ode45", rtol = 1e-06, atol = 1e-06)


y1 <- round(sol[,2]) + rpois(10,1)
y2 <- round(sol[,3]) + rpois(10,1)

smcdar::log_lhood_lotka_volterra_ctmc(theta = 2*theta+0.7, y1 = y1, y2 = y2, times = times, y1_max = 30, y2_max = 30)
smcdar::log_lhood_lotka_volterra_lna(theta = theta+0.7, y1 = y1, y2 = y2, times = times, rtol = 1e-02, atol = 1e-02)


system.time({
  replicate(n = 100,
     smcdar::log_lhood_lotka_volterra_ctmc(theta = theta, y1 = y1, y2 = y2, times = times, y1_max = 30, y2_max = 30)
  )
})

system.time({
  replicate(n = 500,
     smcdar::log_lhood_lotka_volterra_lna(theta = theta*2, y1 = y1, y2 = y2, times = times, tol = list(rtol = 1e-01, atol = 1e-01))
  )
})



system.time({
  replicate(n = 100,
            simulate_ctmc(generator_matrix = G, total_time = 7, initial_state  = 6000)
  )
})
