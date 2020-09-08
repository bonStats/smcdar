#' Run SMC sampler with DA (or possibly not DA methods)
#'
#' @param num_p Number of particles.
#' @param step_scale_set Set  of scales used in MVN kernel.
#' @param use_da Logical. Use Delayed-acceptance?
#' @param use_approx Logical. Use approximate likelihood?
#' @param start_from_approx Logical. Start from approximate likelihood the transition to full with delayed-acceptance?
#' @param refresh_ejd_threshold Distance required to travel for expected jumping distance.
#' @param par_start Approximate likelihood optimisation parameter starting value.
#' @param log_prior Function for log-prior.
#' @param log_like  Function for (full) log-likelihood.
#' @param log_like_approx Function for approximate log-likelihood.
#' @param log_like_approx_comps Function for approximate log-likelihood that returns components.
#' @param draw_prior Function to draw random sample from prior.
#' @param start_from_approx_fit SMC object (list) with approximate fit if start_from_approx is TRUE.
#' @param calibrate_approx_likelihood Function to optimise perform calibration of approximate likelihood within DA.
#' @param find_best_step_scale Function to find best step scale.
#' @param max_anneal_temp Temperature to anneal likelihood (posterior) to. Max 1.
#' @param max_mh_steps Maximum allowable cycles per mutation step.
#' @param use_robust_cov Should the covariance be calculated by \code{robust_mean_cov}. Defaults to \code{FALSE} or equivalently \code{cov.wt}.
#' @param save_post_interface Logical. Should the memoised interface to the likelihood functions be returned (very large).
#' @param cores Number of cores.
#' @param verbose Logical. Should the SMC iterations update be printed to console?
#' @param ... Arguments to be passed to functions.
#'
#' @return A list
#' @export
#'
run_smc_da <- function(num_p, step_scale_set, use_da, use_approx = F, start_from_approx = F, refresh_ejd_threshold, par_start,
                       log_prior, log_like, log_like_approx, log_like_approx_comps = NULL, draw_prior,
                       start_from_approx_fit,
                       calibrate_approx_likelihood,
                       find_best_step_scale,
                       max_anneal_temp = 1,
                       max_mh_steps = 50,
                       use_robust_cov = F,
                       save_post_interface,
                       cores = 1L,
                       verbose = F, ...
){

  stopifnot(
    ! (use_da & use_approx),
    #!start_from_approx | use_da,
    is.numeric(max_anneal_temp),
    max_anneal_temp <= 1,
    max_anneal_temp > 0
  )

  dots <- list(...)

  if( !start_from_approx ){

    log_post_llh_interface <-
      log_likelihood_anneal_func_da(
        log_likelihood = log_like,
        log_like_approx = log_like_approx,
        cwise_log_like_approx = log_like_approx_comps,
        log_prior = log_prior,
        cores = cores
      )

    i <- 1
    ttime <- as.difftime(0, units = "secs") # total time
    temps <- 0

    curr_partl <- particles(draw_prior(num_p))
    log_z <- 0 # Z_1 = 1 due to annealing

  } else {

    log_post_llh_interface <-
      log_likelihood_anneal_func_da(
        log_likelihood = log_like,
        log_like_approx = log_like_approx,
        cwise_log_like_approx = log_like_approx_comps,
        log_prior = log_prior,
        max_approx_anneal = start_from_approx_fit$max_anneal_temp,
        cores = cores,
        start_from_approx = TRUE,
        count_start = start_from_approx_fit$evaluation_counts
      )

    i <- 1
    ttime <- start_from_approx_fit$total_time # total time
    temps <- 0

    curr_partl <- start_from_approx_fit$particles
    log_z <- start_from_approx_fit$log_z

  }

  iter_summary <- list()
  particle_list <- list(curr_partl)

  while(tail(temps,1) < max_anneal_temp){
    stime <- Sys.time()

    curr_log_post <- log_post_llh_interface(
      curr_partl,
      temp = tail(temps, 1),
      type = ifelse(use_approx, "approx_posterior", "full_posterior")
      )

    ## new temp, update posterior
    # use ESS to find temp
    half_ess_temp_obj <- function(temp) {
      log_post <- log_post_llh_interface(curr_partl, temp = temp, type = ifelse(use_approx, "approx_posterior", "full_posterior"))
      log_post <- replace(x = log_post, is.na(log_post), -Inf)

      partl_temp <- curr_partl
      weights(partl_temp, log = T) <- weights(partl_temp, log = T) + log_post - curr_log_post

      res <- ess(partl_temp) - num_p/2

      res <- ifelse(is.finite(res), res, 1e08)
      res
    }

    if(half_ess_temp_obj(max_anneal_temp) > 0){

      temp_optim <- list(root = max_anneal_temp)

    } else {

      temp_optim <- stats::uniroot(f = half_ess_temp_obj,
                                   interval = c(tail(temps,1), max_anneal_temp),
                                   tol = 1e-08
                                   )

    }

    temps <- c( temps, min(max_anneal_temp, temp_optim$root) )

    # new log posterior
    prev_log_post <- curr_log_post
    curr_log_post <- log_post_llh_interface(
      curr_partl,
      temp = tail(temps,1),
      type = ifelse(use_approx, "approx_posterior", "full_posterior")
      )
    curr_log_post <- replace(x = curr_log_post, is.na(curr_log_post), -Inf)

    # weight update
    weights(curr_partl, log = T) <- weights(curr_partl, log = T, normalise = T) + curr_log_post - prev_log_post

    log_z <- log_z + log(sum(weights(curr_partl, log = F, normalise = F)))

    if(use_robust_cov){
      mvn_var <- robust_mean_cov(curr_partl, wt =  weights(curr_partl, normalise = T), method = "unbiased")
    } else {
      mvn_var <- stats::cov.wt(curr_partl, wt =  weights(curr_partl, normalise = T), method = "unbiased")
    }

    ## resample (index also returned, to see duplicates.)
    rs_obj <- resample_stratified(curr_partl, num_strata = 10)
    resampled_partl <- rs_obj$particles # resampled particles

    ## refresh
    sample_step_scale <- step_scale_set[sample.int(num_particles(resampled_partl)) %% length(step_scale_set) + 1]
    proposed_partl <- mvn_jitter(particles = resampled_partl, step_scale = sample_step_scale, var = mvn_var$cov) # jittered particles

    total_artifical_time <- as.difftime(0, units = "secs")

    if(use_da){

      if(i == 1) initial_par_start <- par_start

      # pre-tranformation for approx likelihood
      if(i > 1 & !is.null(calibrate_approx_likelihood) ){

        optim_time <- Sys.time()

        approx_ll_tune_optim <- calibrate_approx_likelihood(
          particles = curr_partl,
          loglike = log_post_llh_interface,
          log_like_approx_comps = log_like_approx_comps,
          par_start = par_start * 0.5,
          initial_par_start = initial_par_start)

        optim_time <- Sys.time() - optim_time
        par_start <- approx_ll_tune_optim$par

      } else {

        optim_time <- NULL
        approx_ll_tune_optim <- list(par = NULL,
                                     trans = identity,
                                     weights = NULL)

      }

      mh_res <- mh_da_step_bglr(new_particles = proposed_partl, old_particles = resampled_partl, loglike = log_post_llh_interface, var = mvn_var$cov, temp = tail(temps,1),  pre_trans = approx_ll_tune_optim$trans, approx_ll_weights = approx_ll_tune_optim$weights, mvn_step_scale = sample_step_scale)

    } else {
      mh_res <- mh_step(new_particles = proposed_partl, old_particles = resampled_partl, loglike = log_post_llh_interface, var = mvn_var$cov, temp = tail(temps,1), type = ifelse(use_approx, "approx_posterior", "full_posterior"))
      optim_time <- NULL
    }

    # if likelihood has artifical timing, need to record total time for stime and etime.
    total_artifical_time <- total_artifical_time + mh_res$extra_artifical_time

    # update particles for 1 iteration
    curr_partl <- replace_particles(new_particles = proposed_partl, old_particles = resampled_partl, index = mh_res$accept)

    # distance threshold start
    expected_sq_dist <- if(use_da) {mh_res$est_expected_dist^2} else {mh_res$expected_dist^2}

    # optimise mh step
    # best_step_scale <- find_best_step_scale(step_scale = sample_step_scale, dist = mh_res$dist, comptime = mh_res$comp_time)
    best_ss <- find_best_step_scale(eta = sample_step_scale,
                                    dist = mh_res$proposal_dist,
                                    adjust_D = stats::median(expected_sq_dist),
                                    prob_accept = if(use_da) {mh_res$est_prob_accept} else {mh_res$prob_accept},
                                    surrogate_expected_acceptance = mh_res$expected_pre_accept,
                                    surrogate_cost = mh_res$avg_surr_like_cost,
                                    full_cost = mh_res$avg_full_like_cost,
                                    da = use_da)

    accept_prop <- mean(mh_res$accept)
    pre_accept_prop <- mean(mh_res$pre_accept)

    mh_step_count <- 1

    # update additional times with best_step_scale
    while( (stats::median(expected_sq_dist) < refresh_ejd_threshold) & (mh_step_count < max_mh_steps) ){

      proposed_partl <- mvn_jitter(particles = curr_partl, step_scale = best_ss$step_scale, var = mvn_var$cov)
      #mh_res <- mh_func(new_particles = proposed_partl, old_particles = curr_partl, var = mvn_var$cov, temp = tail(temps,1) )
      if(use_da){
        mh_res <- mh_da_step_bglr(new_particles = proposed_partl, old_particles = resampled_partl, loglike = log_post_llh_interface, var = mvn_var$cov, temp = tail(temps,1),  pre_trans = approx_ll_tune_optim$trans, time_on = T, approx_ll_weights = approx_ll_tune_optim$weights, reg_warnings = F)
      } else {
        mh_res <- mh_step(new_particles = proposed_partl, old_particles = resampled_partl, loglike = log_post_llh_interface, var = mvn_var$cov, temp = tail(temps,1), type = ifelse(use_approx, "approx_posterior", "full_posterior"), time_on = T)
      }
      curr_partl <- replace_particles(new_particles = proposed_partl, old_particles = curr_partl, index = mh_res$accept)
      accept_prop <- c(accept_prop, mean(mh_res$accept))
      pre_accept_prop <- c(pre_accept_prop, mean(mh_res$pre_accept))

      expected_sq_dist <- if(use_da) {expected_sq_dist + (mh_res$est_expected_dist^2)} else {expected_sq_dist + (mh_res$expected_dist^2)}

      mh_step_count <- mh_step_count + 1
      total_artifical_time <- total_artifical_time + mh_res$extra_artifical_time

    }

    target_log_post <- log_post_llh_interface(curr_partl, temp = 1, type = ifelse(use_approx, "approx_posterior", "full_posterior") )

    etime <- Sys.time()

    ttime <- ttime + difftime(etime, stime, units = "secs") + total_artifical_time

    if(verbose){
      cat("*iter:", i,
          "*temp:", round(tail(temps,1),3), "\n\t",
          "*llh-target:", vec_summary(target_log_post),
          "*step-scale:", best_ss$step_scale,
          "*E-mh-steps:", best_ss$expected_iter,
          "*mh-steps:", mh_step_count,"\n\t",
          "*pre-accept-pr:", round(mean(pre_accept_prop[-1]),3),
          "*accept-pr:", round(mean(accept_prop[-1]),3),
          "*time:", round(difftime(etime, stime, units = "secs") + total_artifical_time, 1),
          "\n")
    }

    iter_summary[[i]] <-
      list(
        temp = tail(temps,1),
        llh_target = target_log_post,
        step_scale = best_ss$step_scale,
        expected_mh_steps = best_ss$expected_iter,
        mh_steps = mh_step_count,
        pre_accept_pr = pre_accept_prop,
        accept_pr = accept_prop,
        approx_ll_calib = if(exists("approx_ll_tune_optim")) list(par = approx_ll_tune_optim$par, weights = approx_ll_tune_optim$weights) else NULL,
        optim_time = optim_time,
        true_time = difftime(etime, stime, units = "secs"),
        time = difftime(etime, stime, units = "secs") + total_artifical_time
      )


    i <- i + 1

    particle_list <- c(particle_list,
                       list(curr_partl)
    )

  }

  # # reweight: #DONT NEED TO
  # prev_log_post <- curr_log_post
  # curr_log_post <- papply(curr_partl, fun = log_post_llh_interface, temp = tail(temps,1), type = ifelse(use_approx, "approx_posterior", "full_posterior") )
  # curr_log_post <- replace(x = curr_log_post, is.na(curr_log_post), -Inf)
  # weights(curr_partl, log = T) <- weights(curr_partl, log = T) + curr_log_post - prev_log_post

  # make function to return approx posterior...

  return(
    list(
      particles = curr_partl,
      log_z = log_z,
      eve_var_est = eve_var_est(curr_partl, log_z = log_z, num_iter = i),
      total_time = ttime,
      evaluation_counts = environment(log_post_llh_interface)$evaluation_counts,
      temps = temps,
      log_post_llh_interface = if(save_post_interface){ log_post_llh_interface} else {NULL},
      iter_summary = iter_summary,
      particle_list = particle_list,
      max_anneal_temp = max_anneal_temp
    )
  )

}

vec_summary <- function(x){

  outv <- round(c(mean(x),range(x)),2)

  paste0(outv[1], " (",outv[2],", ", outv[3],")")

}
