#### Helper ####

vec_summary <- function(x){

  outv <- round(c(mean(x),range(x)),2)

  paste0(outv[1], " (",outv[2],", ", outv[3],")")

}

logit <- function(x) exp(x) / (1 + exp(x))
ilogit <- function(x) log(x) - log(1 - x)

get_hash_memoise <- function(f, x){

  called_args <- list(x = x)
  default_args <- structure(list(), .Names = character(0))
  args <- c(lapply(called_args, eval, parent.frame()),
            lapply(default_args, eval, envir = environment()))

  get("_cache", envir = environment(f))$digest(
    c(
      body(get("_f", envir = environment(f))),
      args,
      list()
      )
    )
}

has_hash_value_memoise <- function(f, hash){

  get("_cache", envir = environment(f))$has_key(hash)

}

####

#### SMC ####

# unified interface to posterior, likelihoods, so that memoising is efficient
log_likelihood_anneal_func_da <- function(log_likelihood, log_like_approx, log_prior){

  mem_log_likelihood <- memoise::memoise(log_likelihood)
  mem_log_likelihood_approx <-  memoise::memoise(log_like_approx)
  mem_log_prior <- memoise::memoise(log_prior)

  unique_x <- NULL

  f <- function(x, temp = 1, lh_trans = identity, type = "full_posterior", ...){

    # make matrix of x values memoised...
    if(type %in% c("full_posterior", "full_approx_lhood_ratio", "full_likelihood")){ # if will evaulate true likelihood

      hash_key <- get_hash_memoise(f = mem_log_likelihood, x = x)

      if( !has_hash_value_memoise(f = mem_log_likelihood, hash = hash_key) ){
        unique_x <<- rbind(unique_x, x) # matrix
      }

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

  return(f)

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

  avg_full_like_cost <- mean(attr(new_loglike_type, "comptime"))

  return(list(
    pre_accept = NA,
    accept = accept,
    dist = dist,
    comp_time = attr(new_loglike_type, "comptime"),
    avg_full_like_cost = avg_full_like_cost,
    avg_surr_like_cost = NA
  ))

}

mh_da_step_bglr <- function(new_particles, old_particles, var, temp, loglike, pre_trans = identity, timing_groups){
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
  full_post_old_particles <- papply(old_particles[accept_1,,drop = F], fun = loglike, type = "full_likelihood", temp = temp, comp_time = F) # should be memoised.
  approx_trans_llh_new_particles <- papply(new_particles[accept_1,,drop = F], fun = loglike, lh_trans = pre_trans, type = "approx_likelihood", temp = temp, comp_time = T)
  approx_trans_llh_old_particles <- papply(old_particles[accept_1,,drop = F], fun = loglike, lh_trans = pre_trans, type = "approx_likelihood", temp = temp, comp_time = T)

  log_rho_tilde_2 <-
    ( full_post_new_particles - full_post_old_particles ) -
    ( approx_trans_llh_new_particles - approx_trans_llh_old_particles)

  rho_2 <- exp(log_rho_tilde_1[accept_1] + log_rho_tilde_2 - log_rho_1[accept_1]) # only for accepted particles

  accept_2 <- rho_2 > runif(sum(accept_1))

  accept <- accept_1
  accept[accept_1] <- accept_2

  maha_dist <- papply(new_particles - old_particles, function(x){t(x) %*% solve(var, x)})

  comp_time <- attr(approx_trans_post_new_particles, "comptime") + attr(approx_trans_post_old_particles, "comptime")

  comp_time[accept_1] <- comp_time[accept_1] +
    attr(full_post_new_particles, "comptime") +
    attr(approx_trans_llh_new_particles, "comptime")

  avg_full_like_cost <- mean( attr(full_post_new_particles, "comptime") )
  avg_surr_like_cost <-
    mean(attr(approx_trans_post_new_particles, "comptime")) +
    mean(attr(approx_trans_post_old_particles, "comptime")) +
    mean(attr(approx_trans_llh_new_particles, "comptime"))

  dist <- sqrt( maha_dist ) * accept

  return(list(
    pre_accept = accept_1,
    accept = accept,
    dist = dist,
    comp_time = comp_time,
    avg_full_like_cost = avg_full_like_cost,
    avg_surr_like_cost = avg_surr_like_cost
  ))

}


time_steps_to_min_quantile_dist_emp <- function(dist, D, rho, max_T = 10){

  if( all(abs(dist) < 1e-04) ) return(list(prob = 0, iter = Inf, sufficient_iter = F))

  nz_index <- abs(dist) > 1e-04
  acceptance_rate <- mean( nz_index )
  ecdf_nz_dist <- ecdf(dist[nz_index])

  for(tT in 1:max_T){

    pr_D <- sum( ( (1 - ecdf_nz_dist(D/1:tT))^(1:tT) ) * dbinom(1:tT, size = tT, prob = acceptance_rate) )
    if(pr_D > rho) break;

  }

  res <- list(prob = pr_D, iter = tT, sufficient_iter = T)

  if(pr_D < rho) res$sufficient_iter <- F

  return(res)

}

time_steps_to_min_quantile_dist_normal <- function(dist, D, rho, max_T = 10){

  if( all(abs(dist) < 1e-04) ) return(list(prob = 0, iter = Inf, sufficient_iter = F))

  nz_index <- abs(dist) > 1e-04
  acceptance_rate <- mean( nz_index )
  n_mean <- mean(dist[nz_index])
  n_sd <- sqrt(var(dist[nz_index]))

  for(tT in 1:max_T){

    pr_D <- sum( (1 - pnorm(D, mean = (1:tT)*n_mean, sd = sqrt(1:tT)*n_sd)) * dbinom(1:tT, size = tT, prob = acceptance_rate) )
    if(pr_D > rho) break;

  }

  res <- list(prob = pr_D, iter = tT, sufficient_iter = T)

  if(pr_D < rho) res$sufficient_iter <- F

  return(res)

}

time_steps_to_min_quantile_dist_gamma <- function(dist, D, rho, max_T = 10){

  if( all(abs(dist) < 1e-04) ) return(list(prob = 0, iter = Inf, sufficient_iter = F))

  nz_index <- abs(dist) > 1e-04
  acceptance_rate <- mean( nz_index )
  gam_fit <- MASS::fitdistr(dist[nz_index], "gamma", lower = c(0,0))
  g_rate <- gam_fit$estimate["rate"]
  g_shape <- gam_fit$estimate["shape"]

  for(tT in 1:max_T){

    pr_D <- sum( (1 - pgamma(D, shape = (1:tT)*g_shape, rate = g_rate)) * dbinom(1:tT, size = tT, prob = acceptance_rate) )
    if(pr_D > rho) break;

  }

  res <- list(prob = pr_D, iter = tT, sufficient_iter = T)

  if(pr_D < rho) res$sufficient_iter <- F

  return(res)

}

time_steps_to_min_quantile_dist_bootstrap <- function(dist, D, rho, max_T = 10, boot_samples = 100){

  if( all(abs(dist) < 1e-04) ) return(list(prob = 0, iter = Inf, sufficient_iter = F))

  nz_index <- abs(dist) > 1e-04
  acceptance_rate <- mean( nz_index )

  boot_prob_gr_D <- vector(mode = "numeric", length = max_T)

  for(tT in 1:max_T){

    boot_prob_gr_D[tT] <- mean(rowSums(matrix(sample(dist[nz_index], size = boot_samples * tT, replace = T), ncol = tT)) > D)
    pr_D <- sum( ( 1 - boot_prob_gr_D[1:tT] ) * dbinom(1:tT, size = tT, prob = acceptance_rate) )
    if(pr_D > rho) break;

  }

  res <- list(prob = pr_D, iter = tT, sufficient_iter = T)

  if(pr_D < rho) res$sufficient_iter <- F

  return(res)

}

min_da_mh_cost <- function(min_T, surrogate_acceptance, surrogate_cost, full_cost){

  stopifnot(length(min_T) == length(surrogate_acceptance),
            length(surrogate_cost) == 1,
            length(full_cost) == 1
            )

  total_costs <- min_T * (surrogate_acceptance * full_cost + surrogate_cost)

  which.min(total_costs)

}

min_mh_cost <- function(min_T){

  which.min(min_T)

}

run_smc_da <- function(num_p, step_scale_set, use_da, use_approx = F, refresh_ejd_threshold, b_s_start,
                       log_prior, log_like, log_like_approx, draw_prior,
                       optimise_pre_approx_llhood_transformation,
                       find_best_step_scale,
                       verbose = F
){

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
  b_s_start <- NULL

  # need to generalise (or take out of function and have pre-defined?) generate_intial?
  curr_partl <- particles(beta = draw_prior(num_p))

  log_z <- 0
  iter_summary <- list()
  particle_list <- list(curr_partl)

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

      res <- ess(partl_temp) - num_p/2

      res <- ifelse(is.finite(res), res, 1e08)
      res
    }

    #temp_optim <- optim(par = tail(temps)+0.0001, fn = half_ess_temp_obj, lower = tail(temps,1), upper = 1, method = "L-BFGS-B")

    if(half_ess_temp_obj(1) > 0){

      temp_optim <- list(root = 1)

    } else {

      temp_optim <- uniroot(f = half_ess_temp_obj, interval = c(tail(temps,1),1), tol = 1e-08)

    }

    temps <- c( temps, ifelse(temp_optim$root > 0.98, 1, temp_optim$root))

    # new log posterior
    prev_log_post <- curr_log_post
    curr_log_post <- papply(curr_partl, fun = log_ann_post_ctmc_da, temp = tail(temps,1), type = ifelse(use_approx, "approx_posterior", "full_posterior") )
    curr_log_post <- replace(x = curr_log_post, is.na(curr_log_post), -Inf)

    # weight update
    weights(curr_partl, log = T) <- weights(curr_partl, log = T) + curr_log_post - prev_log_post

    log_z <- log_z + log(sum(weights(curr_partl, log = F, normalise = F)))

    mvn_var <- cov.wt(curr_partl, wt =  weights(curr_partl), method = "unbiased")

    ## resample (index also returned, to see duplicates.)
    rs_obj <- resample_stratified(curr_partl, num_strata = 10)
    resampled_partl <- rs_obj$particles # resampled particles

    ## refresh
    sample_step_scale <- step_scale_set[sample.int(num_particles(resampled_partl)) %% length(step_scale_set) + 1]
    proposed_partl <- mvn_jitter(particles = resampled_partl, step_scale = sample_step_scale, var = mvn_var$cov) # jittered particles

    if(use_da){

      # pre-tranformation for approx likelihood
      if(tail(temps,1) > 0.05){
        partial_optim <- optimise_pre_approx_llhood_transformation(particles = curr_partl, loglike = log_ann_post_ctmc_da, temp = tail(temps, 1), max_iter = ceiling(tail(temps, 1) * 100), b_s_start = NULL)
        pre_trans <- function(x){ (x - partial_optim$par[1:5]) / exp(partial_optim$par[6:10]) }
        b_s_start <- partial_optim$par
      } else {
        pre_trans <- identity
      }

      mh_res <- mh_da_step_bglr(new_particles = proposed_partl, old_particles = resampled_partl, loglike = log_ann_post_ctmc_da, var = mvn_var$cov, temp = tail(temps,1),  pre_trans = pre_trans)

    } else {
      mh_res <- mh_step(new_particles = proposed_partl, old_particles = resampled_partl, loglike = log_ann_post_ctmc_da, var = mvn_var$cov, temp = tail(temps,1), type = ifelse(use_approx, "approx_posterior", "full_posterior"))
    }

    # update particles for 1 iteration
    curr_partl <- replace_particles(new_particles = proposed_partl, old_particles = resampled_partl, index = mh_res$accept)

    # optimise mh step
   # best_step_scale <- find_best_step_scale(step_scale = sample_step_scale, dist = mh_res$dist, comptime = mh_res$comp_time)
    best_ss <- find_best_step_scale(eta = sample_step_scale,
                                            dist = mh_res$dist,
                                            surrogate_acceptance = mh_res$pre_accept,
                                            surrogate_cost = mh_res$avg_surr_like_cost,
                                            full_cost = mh_res$avg_full_like_cost,
                                            da = use_da)

    accept_prop <- mean(mh_res$accept)
    pre_accept_prop <- mean(mh_res$pre_accept)

    # distance threshold
    total_dist <- mh_res$dist

    mh_step_count <- 1
##### TURN TIMING OFF
    # update additional times with best_step_scale
    while(median(total_dist) < refresh_ejd_threshold){
      proposed_partl <- mvn_jitter(particles = curr_partl, step_scale = best_ss$step_scale, var = mvn_var$cov)
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

    if(verbose){
      cat("*iter:", i,
          "*temp:", round(tail(temps,1),3), "\n\t",
          "*llh-target:", vec_summary(target_log_post),
          "*step-scale:", best_ss$step_scale,
          "*E-mh-steps:", best_ss$expected_iter,
          "*mh-steps:", mh_step_count,"\n\t",
          "*pre-accept-pr:", round(mean(pre_accept_prop[-1]),3),
          "*accept-pr:", round(mean(accept_prop[-1]),3),
          "*time:", round(difftime(etime, stime, units = "secs"),1),
          "\n")
    }

    iter_summary[[i]] <-
                      list(
                           temp = tail(temps,1),
                           llh_target = target_log_post,
                           approx_llh_pre_trans = if(use_da){ pre_trans } else {NULL},
                           step_scale = best_ss$step_scale,
                           expected_mh_steps = best_ss$expected_iter,
                           mh_steps = mh_step_count,
                           pre_accept_pr = pre_accept_prop,
                           accept_pr = accept_prop,
                           ll_opt_par = b_s_start,
                           time = difftime(etime, stime, units = "secs")
                      )


    i <- i + 1

    particle_list <- c(particle_list,
                       list(curr_partl)
    )

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
    log_ann_post_ctmc_da = log_ann_post_ctmc_da,
    iter_summary = iter_summary,
    particle_list = particle_list
  ))

}
