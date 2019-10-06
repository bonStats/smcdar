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
  ecdf_nz_dist <- stats::ecdf(dist)

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

  acceptance_rate <- mean( prob_accept )
  n_mean <- mean(dist)
  n_sd <- sqrt(stats::var(dist))

  for(tT in 1:max_T){

    pr_D <- sum( (1 - stats::pnorm(D, mean = (1:tT)*n_mean, sd = sqrt(1:tT)*n_sd)) * stats::dbinom(1:tT, size = tT, prob = acceptance_rate) )
    if(pr_D > rho) break;

  }

  res <- list(prob = pr_D, iter = tT, sufficient_iter = T)

  if(pr_D < rho) res$sufficient_iter <- F

  return(res)

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

  acceptance_rate <- mean( prob_accept )
  gam_fit <- MASS::fitdistr(dist, "gamma", lower = c(0,0))
  g_rate <- gam_fit$estimate["rate"]
  g_shape <- gam_fit$estimate["shape"]

  for(tT in 1:max_T){

    pr_D <- sum( (1 - stats::pgamma(D, shape = (1:tT)*g_shape, rate = g_rate)) * stats::dbinom(1:tT, size = tT, prob = acceptance_rate) )
    if(pr_D > rho) break;

  }

  res <- list(prob = pr_D, iter = tT, sufficient_iter = T)

  if(pr_D < rho) res$sufficient_iter <- F

  return(res)

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

  acceptance_rate <- mean( prob_accept )

  boot_prob_gr_D <- vector(mode = "numeric", length = max_T)

  for(tT in 1:max_T){

    boot_prob_gr_D[tT] <- mean(rowSums(matrix(sample(dist, size = boot_samples * tT, replace = T), ncol = tT)) > D)
    pr_D <- sum( ( 1 - boot_prob_gr_D[1:tT] ) * stats::dbinom(1:tT, size = tT, prob = acceptance_rate) )
    if(pr_D > rho) break;

  }

  res <- list(prob = pr_D, iter = tT, sufficient_iter = T)

  if(pr_D < rho) res$sufficient_iter <- F

  return(res)

}
