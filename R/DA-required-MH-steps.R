#' Calculate number of MH steps required until threshold is met.
#'
#' Emperical version.
#'
#' @param dist Proposed distances travelled by each particle.
#' @param prob_accept Vector (or mean) of MH acceptance probabilities for proposal.
#' @param D Threshold distance required to travel in total.
#' @param rho Probability level of exceeding trheshold distance.
#' @param max_T Maximum number of iterations to compute before "failing".
#'
#' @return List with results for each parameter.
#' @export
time_steps_to_min_quantile_dist_emp <- function(dist, prob_accept, D, rho, max_T = 10){

  acceptance_rate <- mean( prob_accept )
  ecdf_nz_dist <- stats::ecdf(dist^2)

  for(tT in 1:max_T){

    pr_D <- sum( ( (1 - ecdf_nz_dist(D/1:tT))^(1:tT) ) * stats::dbinom(1:tT, size = tT, prob = acceptance_rate) )
    if(pr_D > rho) break;

  }

  res <- list(prob = pr_D, iter = tT, sufficient_iter = T)

  if(pr_D < rho) res$sufficient_iter <- F

  return(res)

}

#' Calculate number of MH steps required until threshold is met.
#'
#' Normal version.
#'
#' @param dist Proposed distances travelled by each particle.
#' @param prob_accept Vector (or mean) of MH acceptance probabilities for proposal.
#' @param D Threshold distance required to travel in total.
#' @param rho Probability level of exceeding trheshold distance.
#' @param max_T Maximum number of iterations to compute before "failing".
#'
#' @return List with results for each parameter.
#' @export
time_steps_to_min_quantile_dist_normal <- function(dist, prob_accept, D, rho, max_T = 10){

  expected_dist <- prob_accept * (dist^2)
  n_mean <- mean(expected_dist)
  n_sd <- sqrt(stats::var(expected_dist))

  pr_D <- stats::pnorm(D, mean = (1:max_T)*n_mean, sd = sqrt(1:max_T)*n_sd, lower.tail = F)

  suff_Ts <- which(pr_D > rho)

  if(length(suff_Ts) == 0){

    return(
      list(prob = NA, iter = max_T, sufficient_iter = F)
      )

  } else {

    min_suff_Ts <- min(suff_Ts)

    return(
      list(prob = pr_D[min_suff_Ts], iter = min_suff_Ts, sufficient_iter = T)
    )

  }

}

#' Calculate number of MH steps required until threshold is met.
#'
#' Gamma version.
#'
#' @param dist Proposed distances travelled by each particle.
#' @param prob_accept Vector (or mean) of MH acceptance probabilities for proposal.
#' @param D Threshold distance required to travel in total.
#' @param rho Probability level of exceeding trheshold distance.
#' @param max_T Maximum number of iterations to compute before "failing".
#'
#' @return List with results for each parameter.
#' @export
time_steps_to_min_quantile_dist_gamma <- function(dist, prob_accept, D, rho, max_T = 10){

  expected_dist <- prob_accept * (dist^2)

  # use (small bias) mixed type log-moment estimators (some MLE implementations can fail with infinite spike at zero)
  N <- length(expected_dist)
  d_sum <- sum(expected_dist)
  logd_sum <- sum(log(expected_dist))
  d_logd_sum <- sum(expected_dist * log(expected_dist))

  g_shape <- ( d_sum ) / ( d_logd_sum - ( logd_sum * d_sum / N)  )
  g_scale <- (d_logd_sum / N) - ( logd_sum * d_sum / (N^2) )

  pr_D <- stats::pgamma(D, shape = (1:max_T)*g_shape, scale = g_scale, lower.tail = F)

  suff_Ts <- which(pr_D > rho)

  if(length(suff_Ts) == 0){

    return(
      list(prob = NA, iter = max_T, sufficient_iter = F)
    )

  } else {

    min_suff_Ts <- min(suff_Ts)

    return(
      list(prob = pr_D[min_suff_Ts], iter = min_suff_Ts, sufficient_iter = T)
    )

  }

}

#' Calculate number of MH steps required until threshold is met.
#'
#' Bootstrap version.
#'
#' @param dist Proposed distances travelled by each particle.
#' @param prob_accept Vector (or mean) of MH acceptance probabilities for proposal.
#' @param D Threshold distance required to travel in total.
#' @param rho Probability level of exceeding threshold distance.
#' @param max_T Maximum number of iterations to compute before "failing".
#' @param boot_samples Number of bootstrap samples to use.
#'
#' @return List with results for each parameter.
#' @export
time_steps_to_min_quantile_dist_bootstrap <- function(dist, prob_accept, D, rho, max_T = 10, boot_samples = 100){

  expected_dist <- prob_accept * (dist^2)

  pr_D <- vector(mode = "numeric", length = max_T)

  for(tT in 1:max_T){
    pr_D[tT] <- mean(rowSums(matrix(sample(expected_dist, size = boot_samples * tT, replace = T), ncol = tT)) > D)
  }

  suff_Ts <- which(pr_D > rho)

  if(length(suff_Ts) == 0){

    return(
      list(prob = NA, iter = max_T, sufficient_iter = F)
    )

  } else {

    min_suff_Ts <- min(suff_Ts)

    return(
      list(prob = pr_D[min_suff_Ts], iter = min_suff_Ts, sufficient_iter = T)
    )

  }

}

#' Calculate number of MH steps required until threshold is met.
#'
#' Median version.
#'
#' @param dist Proposed distances travelled by each particle.
#' @param prob_accept Vector (or mean) of MH acceptance probabilities for proposal.
#' @param D Threshold distance required to travel in total.
#' @param rho Probability level of exceeding trheshold distance.
#' @param max_T Maximum number of iterations to compute before "failing".
#'
#' @return List with results for each parameter.
#' @export
time_steps_to_min_quantile_dist_median <- function(dist, prob_accept, D, rho, max_T = 10){

  expected_dist <- prob_accept * (dist^2)

  # use the median of the expected_dist
  iter <- as.integer( ceiling( D / median(expected_dist) ) )

  if( !is.finite(iter) ){

    return(
      list(prob = NA, iter = max_T, sufficient_iter = F)
    )

  } else {

    return(
      list(prob = 0.5, iter = iter, sufficient_iter = T)
    )

  }

}
