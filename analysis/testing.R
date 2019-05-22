## change to unit testing

library(Matrix)
library(expm)
theta <- runif(3)
y1_min <- 0
y2_min <- 0
y1_max <- 10
y2_max <- 10

G <- generator_matrix_lotka_volterra_ctmc(theta = theta, y1_min = y1_min, y1_max = y1_max, y2_min = y2_min, y2_max = y2_max)

# check probs
f_init_s <- 11
f_final_s <- 10
time <- 2
P <- expm(time * G, method = "R_Eigen")

rowSums(P)

P[1:3,1:3]

v <- rep(0, times = nrow(P))
v[f_final_s] <- 1

exp_tG_v1 <- expAtv(A = G, t = time, v = v)
exp_tG_v2 <- P %*% v

max(abs(exp_tG_v1$eAtv - exp_tG_v2))

smcdar:::prob_lotka_volterra_ctmc(flat_init_state = f_init_s, flat_final_state = f_final_s,elapsed_time = time, generator_matrix = G)

exp_tG_v1$eAtv[f_init_s]

P[f_init_s,f_final_s]


#### LL ####

y1 <- c(4,5,6,2)
y2 <- c(3,4,6,1)
times <- 1:4

log_lhood_lotka_volterra_ctmc(theta = theta, y1 = y1, y2 = y2, times = times, y1_max = y1_max, y2_max = y2_max)
