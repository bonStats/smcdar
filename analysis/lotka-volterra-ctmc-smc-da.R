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

  sim <- sim[sim[,"time"] < 30,]

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
      log_prior = log_prior,
      use_memoise = T
      )

  draw_prior <- function(n){
    matrix(
      rlnorm(n*3, meanlog = c(-2,-5,-2), sdlog = 0.5), ncol = 3
    )
  }

  ## setup
  num_p <- 20
  curr_temp <- 0

  # particles (must have support on )
  log_theta_start <- draw_prior(num_p)


  partl <- particles(theta = log_theta_start)
  curr_log_post <- papply(partl, fun = log_prior) #temp = 0

  for(i in 2:length(temps)){

    ## new temp
    # use ESS to find temp
    half_ess_temp_obj <- function(temp) {
      log_post <- papply(partl, fun = log_ann_post_ctmc, temp = temp)
      log_post <- replace(x = log_post, is.na(log_post), -Inf)

      partl_temp <- partl
      weights(partl_temp, log = T) <- weights(partl_temp, log = T) + log_post - curr_log_post

      res <- (ess(partl_temp) - num_p/2)^2

      ifelse(is.finite(res), res, 1e08)

    }

    temp_optim <- optim(par = curr_temp + 0.1, fn = half_ess_temp_obj, lower = 0, upper = 1, method = "Brent")


    new_temp <- temp_optim$par

    # new log posterior
    curr_log_post <- papply(partl, fun = log_ann_post_ctmc, temp = new_temp)
    curr_log_post <- replace(x = curr_log_post, is.na(curr_log_post), -Inf)

    # weight update
    weights(partl, log = T) <- weights(partl, log = T) + curr_log_post - prev_log_post

    # resample
    logw <- weights(partl, log = T)
    groups_logw <- cut(logw, breaks = quantile(logw, probs = seq(0,1, length.out = 5)), include.lowest = T)


  }

  #log_post <-

  ## start

  # update weights


