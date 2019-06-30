#### load packages ####
  library(smcdar)
  library(dplyr)
  library(tidyr)
  library(ggplot2)

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

  sim <- sim[sim[,"time"] < 4,]

  vec_summary <- function(x){

    outv <- round(c(mean(x),range(x)),2)

    paste0(outv[1], " (",outv[2],", ", outv[3],")")

  }

#### SMC ####

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
      y1 = sim[,"y1"], y2 = sim[,"y2"],
      times = sim[,"time"]
    ) #* exp(log_theta) # Don't need change of var since this is the likelihood?
  }

  log_ann_post_ctmc <-
    log_likelihood_anneal_func(
      log_likelihood = log_like,
      log_prior = log_prior,
      use_memoise = T
      )

  draw_prior <- function(n){
    matrix(
      rnorm(n*3, mean = c(-2,-5,-2), sd = rep(0.5, 3)), ncol = 3, byrow = T
    )
  }

  mh_step <- function(new_particles, old_particles, var, temp){
    # using MVN (symmetric kernel)
    accept <- exp(
      papply(new_particles, fun = log_ann_post_ctmc, temp = temp) -
      papply(old_particles, fun = log_ann_post_ctmc, temp = temp)
    ) > runif(num_particles(new_particles))

    maha_dist <- papply(new_particles - old_particles, function(x){t(x) %*% solve(var, x)})

    dist <- sqrt( maha_dist ) * accept

    return(list(
      accept = accept,
      dist = dist
    ))

  }

  best_step_scale_dist <- function(step_scale, dist){

    tibble(step_scale = step_scale, dist = dist) %>%
      group_by(step_scale) %>%
      summarise(m = median(dist)) %>%
      {list(
        step_scale = .$step_scale[which.max(.$m)],
        dist = max(.$m)
      )}

  }

  ## setup
  num_p <- 20
  temps <- 0
  step_scale_set <- c(0.3, 0.5, 0.7, 0.8, 0.9)
  ttime <- 0 # total time
  i <- 1

  # particles (must have support on )
  log_theta_start <- draw_prior(num_p)

  curr_partl <- particles(log_theta = log_theta_start)
  curr_log_post <- papply(curr_partl, fun = log_prior) #temp = 0

  while(tail(temps,1) < 1){
    stime <- Sys.time()
    ## new temp, update posterior
    # use ESS to find temp
    half_ess_temp_obj <- function(temp) {
      log_post <- papply(curr_partl, fun = log_ann_post_ctmc, temp = temp)
      log_post <- replace(x = log_post, is.na(log_post), -Inf)

      partl_temp <- curr_partl
      weights(partl_temp, log = T) <- weights(partl_temp, log = T) + log_post - curr_log_post

      res <- (ess(partl_temp) - num_p/2)^2

      ifelse(is.finite(res), res, 1e08)

    }

    temp_optim <- optim(par = (tail(temps,1) + 1)/2, fn = half_ess_temp_obj, lower = tail(temps,1), upper = 1, method = "Brent")

    temps <- c( temps, ifelse(temp_optim$par > 0.95, 1, temp_optim$par))

    # new log posterior
    prev_log_post <- curr_log_post
    curr_log_post <- papply(curr_partl, fun = log_ann_post_ctmc, temp = tail(temps,1) )
    curr_log_post <- replace(x = curr_log_post, is.na(curr_log_post), -Inf)

    # weight update
    weights(curr_partl, log = T) <- weights(curr_partl, log = T) + curr_log_post - prev_log_post

    mvn_var <- cov.wt(curr_partl, wt =  weights(curr_partl), method = "unbiased")

    ## resample (index also returned, to see duplicates.)
    rs_obj <- resample_stratified(curr_partl, num_strata = 10)
    resampled_partl <- rs_obj$particles # resampled particles

    ## refresh
    sample_step_scale <- step_scale_set[sample.int(num_particles(resampled_partl)) %% length(step_scale_set) + 1]
    proposed_partl <- mvn_jitter(particles = resampled_partl, step_scale = sample_step_scale, var = mvn_var$cov) # jittered particles
    mh_res <- mh_step(new_particles = proposed_partl, old_particles = resampled_partl, var = mvn_var$cov, temp = tail(temps,1) )
    # update particles for 1 iteration
    curr_partl <- replace_particles(new_particles = proposed_partl, old_particles = resampled_partl, index = mh_res$accept)
    best_step_scale <- best_step_scale_dist(step_scale = sample_step_scale, dist = mh_res$dist)
    accept_prop <- mean(mh_res$accept)

    # update additional times with best_step_scale
    for(j in 1:2){
      proposed_partl <- mvn_jitter(particles = curr_partl, step_scale = best_step_scale$step_scale, var = mvn_var$cov)
      mh_res <- mh_step(new_particles = proposed_partl, old_particles = curr_partl, var = mvn_var$cov, temp = tail(temps,1) )
      curr_partl <- replace_particles(new_particles = proposed_partl, old_particles = curr_partl, index = mh_res$accept)
      accept_prop <- c(accept_prop, mean(mh_res$accept))
    }

    curr_log_post <- papply(curr_partl, fun = log_ann_post_ctmc, temp = tail(temps,1) )
    target_log_post <-  papply(curr_partl, fun = log_ann_post_ctmc, temp = 1)

    etime <- Sys.time()

    ttime <- ttime + etime-stime

    cat("*iter:", i,
        "*temp:", round(tail(temps,1),4), "\n\t",
        "*llike-target:", vec_summary(target_log_post),
        "*step-scale:", best_step_scale$step_scale,
        "accept-prop:", mean(accept_prop),
        "*time:", round(difftime(etime, stime, units = "secs"),1),
        "\n")

    i <- i + 1
  }


