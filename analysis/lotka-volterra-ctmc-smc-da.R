#### load packages ####
  library(smcdar)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(extraDistr)

#### simulate data ####
  # (birth prey, death prey/ birth pred, death pred)
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

#### SMC ####

  #prior: theta ~ logNormal(-1,1) or log(theta) ~ Normal(-1,1)
  log_prior <- function(theta){

    sum( dlnorm(theta, meanlog = -2, sdlog = 0.5, log = T) )

  }

  log_like <- function(theta){
    log_lhood_lotka_volterra_ctmc(
      theta = theta,
      y1 = sim[,"y1"], y2 = sim[,"y2"],
      y1_max = 100, y2_max = 100,
      times = sim[,"time"]
    )
  }

  log_like_approx <- function(theta){
    log_lhood_lotka_volterra_lna(
      theta = theta,
      y1 = sim[,"y1"], y2 = sim[,"y2"],
      times = sim[,"time"]
    )
  }

  log_ann_post_ctmc <-
    log_likelihood_anneal_func(
      log_likelihood = log_like,
      log_prior = log_prior
      )

  draw_prior <- function(n){
    matrix(
      rlnorm(n*3, meanlog = -2, sdlog = 0.5), ncol = 3
    )
  }

  ## setup
  num_p <- 100

  # particles
  log_theta_start <- draw_prior(num_p)
  partl <- particles(theta = log_theta_start)
  log_post <- papply(partl, fun = log_prior)

  ## start

  # update weights


