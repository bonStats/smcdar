#' Functional for log geometric annealing schedule.
#'
#' \eqn{t \log \pi(x) + (1-t) \log \mu(x)}
#'
#' @param log_posterior A function for the log-posterior.
#' @param log_mu A function for the starting distribution.
#' @param use_memoise Use memoise on posterior and mu functions?
#'
#' @return Function that outputs annealed log-likelihood.
#' @export
log_geometric_anneal_func <- function(log_posterior, log_mu, use_memoise = F){

  stopifnot(is.function(log_posterior),
            is.function(log_mu)
            )

  if(!use_memoise){

    f <- function(x, temp, ...){

      stopifnot(0 <= temp & temp <= 1)

      temp * log_posterior(x, ...) + (1 - temp) * log_mu(x, ...)

    }

  } else {

    mem_log_posterior <- memoise::memoise(log_posterior)
    mem_log_mu <- memoise::memoise(log_mu)

    f <- function(x, temp, ...){

      stopifnot(0 <= temp & temp <= 1)

      temp * mem_log_posterior(x, ...) + (1 - temp) * mem_log_mu(x, ...)

    }

  }

  return(f)

}

#' Functional for log likelihood annealing schedule.
#'
#' \eqn{t \log \pi(x) + \log p(x)}
#'
#' @param log_likelihood A function for the log-likelihood.
#' @param log_prior A function for the starting distribution.
#' @param use_memoise Use memoise on posterior and mu functions?
#'
#' @return Function that outputs annealed log-likelihood.
#' @export
log_likelihood_anneal_func <- function(log_likelihood, log_prior, use_memoise = F){

  stopifnot(is.function(log_likelihood),
            is.function(log_prior)
  )

  if(!use_memoise){

    f <- function(x, temp, ...){

      stopifnot(0 <= temp & temp <= 1)

      temp * log_likelihood(x, ...) +  log_prior(x, ...)

    }

  } else {

    mem_log_likelihood <- memoise::memoise(log_likelihood)
    mem_log_prior <- memoise::memoise(log_prior)

    f <- function(x, temp, ...){

      stopifnot(0 <= temp & temp <= 1)

      temp * mem_log_likelihood(x, ...) + mem_log_prior(x, ...)

    }

  }

  return(f)

}
