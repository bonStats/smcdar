library(smcdar)

# simulate data

true_theta <- c(0.2,0.015,0.4)

sim <- simulate_lotka_volterra_ctmc(true_theta, y1_min = 0, y2_min = 0, y1_max =100, y2_max=100, initial_y1 = 20, initial_y2 =15, total_time=50)


Gmat <- generator_matrix_lotka_volterra_ctmc(theta = true_theta, y1_min = 0, y2_min = 0, y1_max =100, y2_max=100)

simulate_ctmc(Gmat, initial_state = prod(20,15), total_time = 10)

smcdar:::flattened_state_lotka_volterra(18,0, 0,100,0,100)


plot(sim[,"time"],sim[,"y1"], type = "l", col = "blue", ylim = c(0,40))

lines(sim[,"time"],sim[,"y2"], type = "l", col = "red")
unique(diff(sim[,"y2"]) - diff(sim[,"y1"]))
