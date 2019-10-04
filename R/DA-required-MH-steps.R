#' Calculate number of MH steps required until threshold is met.
#'
#' Emperical version.
#'
#' @param dist (Expected) distances travelled by each particle.
#' @param D Threshold distance required to travel in total.
#' @param rho Probability level of exceeding trheshold distance.
#' @param max_T Maximum number of iterations to compute before "failing".
#'
#' @return List with results for each parameter.
#' @export
time_steps_to_min_quantile_dist_emp <- function(dist, D, rho, max_T = 10){

  nz_index <- abs(dist) > 1e-04

  # not enough non-zeros to estimate
  if( sum(nz_index) < 3 )  return(list(prob = 0, iter = Inf, sufficient_iter = F))

  acceptance_rate <- mean( nz_index )
  ecdf_nz_dist <- stats::ecdf(dist[nz_index])

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
#' @param dist (Expected) distances travelled by each particle.
#' @param D Threshold distance required to travel in total.
#' @param rho Probability level of exceeding trheshold distance.
#' @param max_T Maximum number of iterations to compute before "failing".
#'
#' @return List with results for each parameter.
#' @export
time_steps_to_min_quantile_dist_normal <- function(dist, D, rho, max_T = 10){

  nz_index <- abs(dist) > 1e-04

  # not enough non-zeros to estimate
  if( sum(nz_index) < 4 )  return(list(prob = 0, iter = Inf, sufficient_iter = F))

  acceptance_rate <- mean( nz_index )
  n_mean <- mean(dist[nz_index])
  n_sd <- sqrt(stats::var(dist[nz_index]))

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
#' @param dist (Expected) distances travelled by each particle.
#' @param D Threshold distance required to travel in total.
#' @param rho Probability level of exceeding trheshold distance.
#' @param max_T Maximum number of iterations to compute before "failing".
#'
#' @return List with results for each parameter.
#' @export
time_steps_to_min_quantile_dist_gamma <- function(dist, D, rho, max_T = 10){

  nz_index <- abs(dist) > 1e-04

  # not enough non-zeros to estimate
  if( sum(nz_index) < 4 )  return(list(prob = 0, iter = Inf, sufficient_iter = F))

  acceptance_rate <- mean( nz_index )
  gam_fit <- MASS::fitdistr(dist[nz_index], "gamma", lower = c(0,0))
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
#' @param dist (Expected) distances travelled by each particle.
#' @param D Threshold distance required to travel in total.
#' @param rho Probability level of exceeding threshold distance.
#' @param max_T Maximum number of iterations to compute before "failing".
#' @param boot_samples Number of bootstrap samples to use.
#'
#' @return List with results for each parameter.
#' @export
time_steps_to_min_quantile_dist_bootstrap <- function(dist, D, rho, max_T = 10, boot_samples = 100){

  nz_index <- abs(dist) > 1e-04

  # not enough non-zeros to estimate
  if( sum(nz_index) < 3 )  return(list(prob = 0, iter = Inf, sufficient_iter = F))

  acceptance_rate <- mean( nz_index )

  boot_prob_gr_D <- vector(mode = "numeric", length = max_T)

  for(tT in 1:max_T){

    boot_prob_gr_D[tT] <- mean(rowSums(matrix(sample(dist[nz_index], size = boot_samples * tT, replace = T), ncol = tT)) > D)
    pr_D <- sum( ( 1 - boot_prob_gr_D[1:tT] ) * stats::dbinom(1:tT, size = tT, prob = acceptance_rate) )
    if(pr_D > rho) break;

  }

  res <- list(prob = pr_D, iter = tT, sufficient_iter = T)

  if(pr_D < rho) res$sufficient_iter <- F

  return(res)

}
