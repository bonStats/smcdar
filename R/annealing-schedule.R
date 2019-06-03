#' Functional for log geometric annealing schedule.
#'
#' \eqn{t \log \pi(x) + (1-t) \log \mu(x)}
#'
#' @param log_posterior A function for the log-posterior.
#' @param log_mu A function for the starting distribution.
#'
#' @return Function that outputs annealed log-likelihood.
#' @export
log_geometric_anneal_func <- function(log_posterior, log_mu){

  stopifnot(is.function(log_posterior),
            is.function(log_mu)
            )

  function(x, temp, ...){

    stopifnot(0 <= temp & temp <= 1)

    temp * log_posterior(x, ...) + (1 - temp) * log_mu(x, ...)

  }


}

#' Functional for log likelihood annealing schedule.
#'
#' \eqn{t \log \pi(x) + \log p(x)}
#'
#' @param log_likelihood A function for the log-likelihood.
#' @param log_prior A function for the starting distribution.
#'
#' @return Function that outputs annealed log-likelihood.
#' @export
log_likelihood_anneal_func <- function(log_likelihood, log_prior){

  stopifnot(is.function(log_likelihood),
            is.function(log_prior)
  )

  function(x, temp, ...){

    temp * log_likelihood(x, ...) +  log_prior(x, ...)

  }


}
