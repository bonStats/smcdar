#' Generic interface to annealed log-likeihood functions with memoisation: Likelihood annealing.
#'
#' Prior -> Posterior
#'
#' @param log_likelihood (Full) log-likelihood function.
#' @param log_like_approx Approximate log-likelihood function.
#' @param log_prior Log-prior function.
#'
#' @return Value of log-likelihood function.
#' @export
log_likelihood_anneal_func_da <- function(log_likelihood, log_like_approx, log_prior){

  mem_log_likelihood <- memoise::memoise(log_likelihood)
  mem_log_likelihood_approx <-  memoise::memoise(log_like_approx)
  mem_log_prior <- memoise::memoise(log_prior)

  unique_x <- NULL
  unique_llh_val <- NULL
  x_arg_name <- names(formals(mem_log_likelihood))

  evaluation_counts <- setNames( rep(list(0L),2), c("log_likelihood", "log_like_approx"))

  f <- function(x, temp = 1, lh_trans = identity, type = "full_posterior", ...){

    # make matrix of x values memoised, and count log-like evals
    if(type %in% c("full_posterior", "full_approx_lhood_ratio", "full_likelihood")){ # if will evaulate true likelihood

      hash_key <- get_hash_memoise(f = mem_log_likelihood, x = x, arg_name = x_arg_name)

      if( !has_hash_value_memoise(f = mem_log_likelihood, hash = hash_key) ){
        unique_x <<- rbind(unique_x, x) # matrix
        # count evals: log-like
        evaluation_counts$log_likelihood <<- evaluation_counts$log_likelihood + 1
      }

    }

    # count evals
    if(type %in% c("approx_posterior","approx_likelihood","full_approx_lhood_ratio")){

      hash_key <- get_hash_memoise(f = mem_log_likelihood_approx, x = x, arg_name = x_arg_name)

      if( !has_hash_value_memoise(f = mem_log_likelihood_approx, hash = hash_key) ){
        # count evals: approx log-like
        evaluation_counts$log_like_approx <<- evaluation_counts$log_like_approx + 1
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

#' Generic interface to annealed log-likeihood functions with memoisation: Approximate likelihood annealing.
#'
#' Approximate posterior -> Full posterior
#'
#' @param log_likelihood (Full) log-likelihood function.
#' @param log_like_approx Approximate log-likelihood function.
#' @param log_prior Log-prior function.
#'
#' @return Value of log-likelihood function.
#' @export
log_approx_likelihood_anneal_func_da <- function(log_likelihood, log_like_approx, log_prior){

  mem_log_likelihood <- memoise::memoise(log_likelihood)
  mem_log_likelihood_approx <-  memoise::memoise(log_like_approx)
  mem_log_prior <- memoise::memoise(log_prior)

  unique_x <- NULL
  unique_llh_val <- NULL
  x_arg_name <- names(formals(mem_log_likelihood))

  f <- function(x, temp = 1, lh_trans = identity, type = "full_posterior", ...){

    # make matrix of x values memoised, and count log-like evals
    if(type %in% c("full_posterior", "full_approx_lhood_ratio", "full_likelihood")){ # if will evaulate true likelihood

      hash_key <- get_hash_memoise(f = mem_log_likelihood, x = x, arg_name = x_arg_name)

      if( !has_hash_value_memoise(f = mem_log_likelihood, hash = hash_key) ){
        unique_x <<- rbind(unique_x, x) # matrix
        # count evals: log-like
        evaluation_counts$log_likelihood <<- evaluation_counts$log_likelihood + 1
      }

    }

    # count evals
    if(type %in% c("approx_posterior","approx_likelihood","full_approx_lhood_ratio")){

      hash_key <- get_hash_memoise(f = mem_log_likelihood_approx, x = x, arg_name = x_arg_name)

      if( !has_hash_value_memoise(f = mem_log_likelihood_approx, hash = hash_key) ){
        # count evals: approx log-like
        evaluation_counts$log_like_approx <<- evaluation_counts$log_like_approx + 1
      }
    }

    switch(type,
           approx_posterior = temp * mem_log_likelihood_approx(lh_trans(x), ...) + ( 1 - temp ) * mem_log_likelihood_approx(x, ...) + mem_log_prior(x, ...),
           full_approx_lhood_ratio = temp * (  mem_log_likelihood(x, ...) - mem_log_likelihood_approx(lh_trans(x), ...) ),
           approx_likelihood = temp * mem_log_likelihood_approx(lh_trans(x), ...),
           full_likelihood =  temp * mem_log_likelihood(x, ...),
           full_posterior = temp * mem_log_likelihood(x, ...) + ( 1 - temp ) * mem_log_likelihood_approx(x, ...) + mem_log_prior(x, ...),
           prior =  mem_log_prior(x, ...)
    )

  }

  return(f)

}

# extra helper functions

get_hash_memoise <- function(f, x, arg_name){

  called_args <- stats::setNames(list(x), arg_name)
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
