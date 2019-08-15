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

  sim <- sim[sim[,"time"] < 0.5,]

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

  #thin_i <- 1:floor(nrow(sim) / 4)
  thin_i <- 1:nrow(sim)

  log_like_approx <- function(log_theta){
    log_lhood_lotka_volterra_lna(
      theta = exp(log_theta),
      y1 = sim[thin_i,"y1"], y2 = sim[thin_i,"y2"],
      times = sim[thin_i,"time"]
    ) #* exp(log_theta) # Don't need change of var since this is the likelihood?
  }

  # log_ann_post_ctmc <-
  #   log_likelihood_anneal_func(
  #     log_likelihood = log_like,
  #     log_prior = log_prior,
  #     use_memoise = T
  #     )

  # unified interface to posterior, likelihoods, so that memoising is efficient
  log_likelihood_anneal_func_da <- function(log_likelihood, log_like_approx, log_prior){

    mem_log_likelihood <- memoise::memoise(log_likelihood)
    mem_log_likelihood_approx <-  memoise::memoise(log_like_approx)
    mem_log_prior <- memoise::memoise(log_prior)

    unique_x <- NULL

    f <- function(x, temp = 1, lh_trans = identity, type = "full_posterior", ...){

      if(type %in% c("full_posterior", "full_approx_lhood_ratio", "full_likelihood", "full_posterior")){ # if will evaulate true likelihood
        new_unique_x <- rbind(unique_x, x) # matrix
        unique_x <<- unique(new_unique_x) # unique rows of matrix
      }

      switch(type,
        approx_posterior = temp * mem_log_likelihood_approx(lh_trans(x), ...) + mem_log_prior(x, ...),
        full_approx_lhood_ratio = temp * (  mem_log_likelihood(x, ...) - mem_log_likelihood_approx(lh_trans(x), ...) ),
        approx_likelihood = temp * mem_log_likelihood_approx(lh_trans(x), ...),
        full_likelihood =  temp * mem_log_likelihood(x, ...),
        full_posterior = temp * mem_log_likelihood(x, ...) + mem_log_prior(x, ...),
        prior =  mem_log_prior(x, ...)
      )

    }

  }

  draw_prior <- function(n){
    matrix(
      rnorm(n*3, mean = c(-2,-5,-2), sd = rep(0.5, 3)), ncol = 3, byrow = T
    )
  }

  mh_step <- function(new_particles, old_particles, var, temp, loglike, type){
    # using MVN (symmetric kernel)

    new_loglike_type <- papply(new_particles, fun = loglike, temp = temp, type = type, comp_time = T)
    old_loglike_type <- papply(old_particles, fun = loglike, temp = temp, type = type, comp_time = F)

    accept <- exp(
      new_loglike_type - old_loglike_type
    ) > runif(num_particles(new_particles))

    maha_dist <- papply(new_particles - old_particles, function(x){t(x) %*% solve(var, x)})

    dist <- sqrt( maha_dist ) * accept

    return(list(
      pre_accept = NA,
      accept = accept,
      dist = dist,
      comp_time = attr(new_loglike_type, "comptime")
    ))

  }

  mh_da_step <- function(new_particles, old_particles, var, temp, loglike, pre_trans = identity){
    # using MVN (symmetric kernel)

    trans_new_particles <- t(papply(new_particles, pre_trans))
    trans_old_particles <- t(papply(old_particles, pre_trans))

    accept_1 <- exp(

      papply(trans_new_particles, fun = loglike, type = "approx_likelihood", temp = temp) +
        papply(new_particles, fun = loglike, type = "prior")  -
        (
        papply(trans_old_particles, fun = loglike, type = "approx_likelihood", temp = temp) +
        papply(old_particles, fun = loglike, type = "prior")
        )

    ) > runif(num_particles(new_particles))

    accept_2 <- exp(
      papply(new_particles[accept_1,,drop = F], fun = loglike, type = "full_likelihood", temp = temp) -
        papply(old_particles[accept_1,,drop = F], fun = loglike, type = "full_likelihood", temp = temp) -
        (
          papply(trans_new_particles[accept_1,,drop = F], fun = loglike, type = "approx_likelihood", temp = temp) -
            papply(trans_old_particles[accept_1,,drop = F], fun = loglike, type = "approx_likelihood", temp = temp)
        )
    ) > runif(sum(accept_1))

    accept <- accept_1
    accept[accept_1] <- accept_2

    maha_dist <- papply(new_particles - old_particles, function(x){t(x) %*% solve(var, x)})

    dist <- sqrt( maha_dist ) * accept

    return(list(
      pre_accept = accept_1,
      accept = accept,
      dist = dist
    ))

  }


  mh_da_step_bglr <- function(new_particles, old_particles, var, temp, loglike, pre_trans = identity){
    # using MVN (symmetric kernel)

    # uses http://aimsciences.org//article/doi/10.3934/fods.2019005
    # to give every particle some small chance of progressing
    c_const <- 0.01 # how to choose?
    d_const <- 2
    log_b_const <- ( 1 / ( d_const - 1) ) * log(c_const)

    #trans_new_particles <- t(papply(new_particles, pre_trans))
    #trans_old_particles <- t(papply(old_particles, pre_trans))

    approx_trans_post_new_particles <- papply(new_particles, fun = loglike, lh_trans = pre_trans, type = "approx_posterior", temp = temp, comp_time = T)
    approx_trans_post_old_particles <- papply(old_particles, fun = loglike, lh_trans = pre_trans, type = "approx_posterior", temp = temp, comp_time = T)

    # approximate likelihood threshold (surrogate model)
    log_rho_tilde_1 <- approx_trans_post_new_particles - approx_trans_post_old_particles

    # conservative version
    log_rho_1 <- pmin( -log_b_const, pmax( log_b_const, log_rho_tilde_1 ) )

    accept_1 <- exp(log_rho_1) > runif(num_particles(new_particles))

    full_post_new_particles <- papply(new_particles[accept_1,,drop = F], fun = loglike, type = "full_likelihood", temp = temp, comp_time = T)
    full_post_old_particles <- papply(old_particles[accept_1,,drop = F], fun = loglike, type = "full_likelihood", temp = temp, comp_time = F)
    approx_trans_llh_new_particles <- papply(new_particles[accept_1,,drop = F], fun = loglike, lh_trans = pre_trans, type = "approx_likelihood", temp = temp)
    approx_trans_llh_old_particles <- papply(old_particles[accept_1,,drop = F], fun = loglike, lh_trans = pre_trans, type = "approx_likelihood", temp = temp)

    log_rho_tilde_2 <-
      ( full_post_new_particles - full_post_old_particles ) -
      ( approx_trans_llh_new_particles - approx_trans_llh_old_particles)

    rho_2 <- exp(log_rho_tilde_1[accept_1] + log_rho_tilde_2 - log_rho_1[accept_1]) # only for accepted particles

    accept_2 <- rho_2 > runif(sum(accept_1))

    accept <- accept_1
    accept[accept_1] <- accept_2

    maha_dist <- papply(new_particles - old_particles, function(x){t(x) %*% solve(var, x)})

    comp_time <- attr(approx_trans_post_new_particles, "comptime") + attr(approx_trans_post_old_particles, "comptime")

    comp_time[accept_1] <- comp_time[accept_1] + attr(full_post_new_particles, "comptime")

    dist <- sqrt( maha_dist ) * accept

    return(list(
      pre_accept = accept_1,
      accept = accept,
      dist = dist,
      comp_time = comp_time
    ))

  }

  optimise_pre_approx_llhood_transformation <- function(b_s_start, log_theta, loglike){

    #log_theta is matrix
    # should only send values of log_like that have been cached already
    true_llhood <- apply(log_theta, MARGIN = 1, loglike, type = "full_likelihood")

    if(missing(b_s_start)){

      approx_llhood <-apply(log_theta, MARGIN = 1, loglike, type = "approx_likelihood")
      mode_true <- which.max(true_llhood)
      mode_approx <- which.max(approx_llhood)
      b_s_start <- c(log_theta[mode_approx,] - log_theta[mode_true,],0,0,0)

    }

    which_upper <- which( true_llhood > quantile(true_llhood, probs = 0.1) )

    use_index <- sample(which_upper, size = 20, replace = FALSE)

    topt <- function(b_s){

      b_mat <- matrix(b_s[1:3], ncol = 3, nrow = length(use_index), byrow = T)
      s_mat <- matrix(b_s[4:6], ncol = 3, nrow = length(use_index), byrow = T)

      llhood_approx <- apply((log_theta[use_index,] - b_mat)/exp(s_mat), MARGIN = 1, loglike, type = "approx_likelihood")

      keep_ll <- is.finite(llhood_approx)

      # chi-square discrepancy
      sum( # p * (p/q - 1)^2
        true_llhood[use_index[keep_ll]] * (exp(llhood_approx[keep_ll] - true_llhood[use_index[keep_ll]]) - 1)^2
      )

    }

    optim(par = b_s_start, fn = topt, control = list(maxit = 5))

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

  best_step_scale_ejd_v_time <- function(step_scale, dist, comptime){

    tibble(step_scale = step_scale, dist = dist, comptime = comptime) %>%
      group_by(step_scale) %>%
      summarise(m = median(dist/comptime)) %>%
      {list(
        step_scale = .$step_scale[which.max(.$m)],
        dist = max(.$m)
      )}


  }

  ## setup
  # num_p <- 50
  # step_scale_set <- c(0.3, 0.5, 0.7, 0.8, 0.9)
  # use_da <- T
  # use_approx <- F
  # b_s_start <- c(0,0,0,0,0,0)

smc_lotka_volterra_da <- function(num_p, step_scale_set, use_da, use_approx = F, refresh_ejd_threshold, b_s_start = c(0,0,0,0,0,0) ){

  stopifnot(! (use_da & use_approx) )

  log_ann_post_ctmc_da <-
    log_likelihood_anneal_func_da(
      log_likelihood = log_like,
      log_like_approx = log_like_approx,
      log_prior = log_prior
    )

  i <- 1
  ttime <- 0 # total time
  temps <- 0
  optima_matrix <- NULL

  # particles (must have support on )
  log_theta_start <- draw_prior(num_p)

  curr_partl <- particles(log_theta = log_theta_start)

  log_z <- log(num_p)

  while(tail(temps,1) < 1){
    stime <- Sys.time()

    curr_log_post <- papply(curr_partl, fun = log_ann_post_ctmc_da, temp = tail(temps,1), type = ifelse(use_approx, "approx_posterior", "full_posterior"))

    ## new temp, update posterior
    # use ESS to find temp
    half_ess_temp_obj <- function(temp) {
      log_post <- papply(curr_partl, fun = log_ann_post_ctmc_da, temp = temp, type = ifelse(use_approx, "approx_posterior", "full_posterior"))
      log_post <- replace(x = log_post, is.na(log_post), -Inf)

      partl_temp <- curr_partl
      weights(partl_temp, log = T) <- weights(partl_temp, log = T) + log_post - curr_log_post

      res <- (ess(partl_temp) - num_p/2)^2

      res <- ifelse(is.finite(res), res, 1e08)
      res
    }

    temp_optim <- optim(par = (tail(temps,1) + 1)/2, fn = half_ess_temp_obj, lower = tail(temps,1), upper = 1, method = "Brent")

    temps <- c( temps, ifelse(temp_optim$par > 0.98, 1, temp_optim$par))

    # new log posterior
    prev_log_post <- curr_log_post
    curr_log_post <- papply(curr_partl, fun = log_ann_post_ctmc_da, temp = tail(temps,1), type = ifelse(use_approx, "approx_posterior", "full_posterior") )
    curr_log_post <- replace(x = curr_log_post, is.na(curr_log_post), -Inf)

    # weight update
    weights(curr_partl, log = T) <- weights(curr_partl, log = T) + curr_log_post - prev_log_post

    log_z <- log_z + log(sum(weights(curr_partl, log = F, normalise = F))) - log(num_p)

    mvn_var <- cov.wt(curr_partl, wt =  weights(curr_partl), method = "unbiased")

    ## resample (index also returned, to see duplicates.)
    rs_obj <- resample_stratified(curr_partl, num_strata = 10)
    resampled_partl <- rs_obj$particles # resampled particles

    ## refresh
    sample_step_scale <- step_scale_set[sample.int(num_particles(resampled_partl)) %% length(step_scale_set) + 1]
    proposed_partl <- mvn_jitter(particles = resampled_partl, step_scale = sample_step_scale, var = mvn_var$cov) # jittered particles

    if(use_da){
      # pre transformation
      if(tail(temps,1) < 0.1){
        partial_optim <- optimise_pre_approx_llhood_transformation(log_theta = curr_partl, loglike = log_ann_post_ctmc_da)
      } else {
        partial_optim <- optimise_pre_approx_llhood_transformation(b_s_start = b_s_start, log_theta = curr_partl, loglike = log_ann_post_ctmc_da)
      }
      #optima_matrix <- rbind(optima_matrix, partial_optim$par)
      #pre_trans_pars <- colMeans(optima_matrix)
      pre_trans <- function(x){ (x - partial_optim$par[1:3]) / exp(partial_optim$par[4:6]) }
      b_s_start <- partial_optim$par
      mh_res <- mh_da_step_bglr(new_particles = proposed_partl, old_particles = resampled_partl, loglike = log_ann_post_ctmc_da, var = mvn_var$cov, temp = tail(temps,1),  pre_trans = pre_trans)
    } else {
      mh_res <- mh_step(new_particles = proposed_partl, old_particles = resampled_partl, loglike = log_ann_post_ctmc_da, var = mvn_var$cov, temp = tail(temps,1), type = ifelse(use_approx, "approx_posterior", "full_posterior"))
    }

    # update particles for 1 iteration
    curr_partl <- replace_particles(new_particles = proposed_partl, old_particles = resampled_partl, index = mh_res$accept)

    # optimise mh step
    best_step_scale <- best_step_scale_ejd_v_time(step_scale = sample_step_scale, dist = mh_res$dist, comptime = mh_res$comp_time)
    accept_prop <- mean(mh_res$accept)
    pre_accept_prop <- mean(mh_res$pre_accept)

    # distance threshold
    total_dist <- mh_res$dist

    mh_step_count <- 1

    # update additional times with best_step_scale
    while(median(total_dist) < refresh_ejd_threshold){
      proposed_partl <- mvn_jitter(particles = curr_partl, step_scale = best_step_scale$step_scale, var = mvn_var$cov)
      #mh_res <- mh_func(new_particles = proposed_partl, old_particles = curr_partl, var = mvn_var$cov, temp = tail(temps,1) )
      if(use_da){
        mh_res <- mh_da_step_bglr(new_particles = proposed_partl, old_particles = resampled_partl, loglike = log_ann_post_ctmc_da, var = mvn_var$cov, temp = tail(temps,1),  pre_trans = pre_trans)
      } else {
        mh_res <- mh_step(new_particles = proposed_partl, old_particles = resampled_partl, loglike = log_ann_post_ctmc_da, var = mvn_var$cov, temp = tail(temps,1), type = ifelse(use_approx, "approx_posterior", "full_posterior"))
      }
      curr_partl <- replace_particles(new_particles = proposed_partl, old_particles = curr_partl, index = mh_res$accept)
      accept_prop <- c(accept_prop, mean(mh_res$accept))
      pre_accept_prop <- c(pre_accept_prop, mean(mh_res$pre_accept))

      total_dist <- total_dist + mh_res$dist
      mh_step_count <- mh_step_count + 1
    }

    target_log_post <-  papply(curr_partl, fun = log_ann_post_ctmc_da, temp = 1, type = ifelse(use_approx, "approx_posterior", "full_posterior") )

    etime <- Sys.time()

    ttime <- ttime + etime - stime

    cat("*iter:", i,
        "*temp:", round(tail(temps,1),3), "\n\t",
        "*llh-target:", vec_summary(target_log_post),
        "*step-scale:", best_step_scale$step_scale,
        "*step-scale-dist:", round(best_step_scale$dist,3),
        "*mh-steps:", mh_step_count,
        "*pre-accept-pr:", round(mean(pre_accept_prop),3),
        "*accept-pr:", round(mean(accept_prop),3),
        "*time:", round(difftime(etime, stime, units = "secs"),1),
        "\n")

    i <- i + 1
  }

  # reweight:
  prev_log_post <- curr_log_post
  curr_log_post <- papply(curr_partl, fun = log_ann_post_ctmc_da, temp = tail(temps,1), type = ifelse(use_approx, "approx_posterior", "full_posterior") )
  curr_log_post <- replace(x = curr_log_post, is.na(curr_log_post), -Inf)
  weights(curr_partl, log = T) <- weights(curr_partl, log = T) + curr_log_post - prev_log_post

  # make function to return approx posterior...

  return(list(
    particles = curr_partl,
    log_z = log_z,
    eve_var_est = eve_var_est(curr_partl, log_z = log_z, num_iter = i),
    total_time = ttime,
    temps = temps,
    log_ann_post_ctmc_da = log_ann_post_ctmc_da
  ))

}

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

# plot how well log-likelihood approximation is doing each interation.

# CHRIS:
# double annealing?
# check and see if subsampling the lotka-volterra approximation is a big speed up
# other models with high gains from approx likelihood
# ensemble particle filter paper + DA
# SMC squared speed up?

# implement/check eve estimator
# refactor for next apps

#https://xianblog.wordpress.com/2019/08/11/delayed-acceptance-ada-boosted/

# UPDATES:
# could make optimisation for approx likelihood only run when approx vs full posterior discrepancy is beyond some cutoff?
# check and see if subsampling the lotka-volterra approximation is a big speed up

# IDEAS:
# data subsampling with non-iid. Find good subset to base off. Hierarchal logistic or something.
