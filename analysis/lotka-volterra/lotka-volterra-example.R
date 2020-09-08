
true_theta <- c(0.17,0.01,0.2)

sim <- simulate_lotka_volterra_ctmc(
  true_theta,
  y1_min = 0,
  y2_min = 0,
  y1_max =100,
  y2_max=100,
  initial_y1 = 40,
  initial_y2 =25,
  total_time=50
)

tidy_sim <- sim %>% as.data.frame() %>%
  gather(pop,level,-time) %>%
  mutate(pop_name = ifelse(pop == "y1", "prey", "pred"))

tidy_sim %>% ggplot() +
  geom_line(aes(x = time, y = level, colour = pop_name)) +
  theme_bw()

sim <- sim[1:100,]

#### Priors, likelihoods ####

#prior: theta ~ logNormal(-1,0.5) or log(theta) ~ Normal(-1,0.5)
log_prior <- function(log_theta){

  sum( dnorm(log_theta, mean = c(-2,-5,-2), sd = rep(0.5, 3), log = T) )

}

log_like <- function(log_theta){
  log_lhood_lotka_volterra_ctmc(
    theta = exp(log_theta),
    y1 = sim[,"y1"], y2 = sim[,"y2"],
    y1_max = 100, y2_max = 100,
    times = sim[,"time"]
  ) #* exp(log_theta) # Don't need change of var since this is the likelihood?
}

log_like_approx <- function(log_theta){
  log_lhood_lotka_volterra_lna(
    theta = exp(log_theta),
    y1 = sim[1:g_pars$N_approx,"y1"], y2 = sim[1:g_pars$N_approx,"y2"],
    times = sim[1:g_pars$N_approx,"time"]
  )
}

####

source("analysis/lotka-volterra-ctmc-smc-da.R")

f_pars <-  list(
  num_p = 100,
  step_scale_set = c(0.01, 0.05, 0.1, 0.2),
  b_s_start = c(0,0,0,0,0,0),
  refresh_ejd_threshold = 0.01
)

# list2env(f_pars, envir = globalenv()); use_da = T; use_approx = F

smc_da <- with(f_pars, smc_lotka_volterra_da(
  use_da = T, use_approx = F,
  num_p = num_p,
  step_scale_set = step_scale_set,
  b_s_start = b_s_start,
  refresh_ejd_threshold = refresh_ejd_threshold
)
)

smc_full <- with(f_pars, smc_lotka_volterra_da(
  use_da = F, use_approx = F,
  num_p = num_p,
  step_scale_set = step_scale_set,
  b_s_start = b_s_start,
  refresh_ejd_threshold = refresh_ejd_threshold
)
)

smc_approx <- with(f_pars, smc_lotka_volterra_da(
  use_da = F, use_approx = T,
  num_p = num_p,
  step_scale_set = step_scale_set,
  b_s_start = b_s_start,
  refresh_ejd_threshold = refresh_ejd_threshold
)
)

sapply(list(smc_da,smc_full,smc_approx), FUN = getElement, name = "total_time")





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
