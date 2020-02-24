#' Generic interface to annealed log-likeihood functions with memoisation: Likelihood annealing.
#'
#' Prior -> Posterior
#'
#' @param log_likelihood (Full) log-likelihood function.
#' @param log_like_approx Approximate log-likelihood function.
#' @param log_prior Log-prior function.
#' @param cores Cores to use.
#'
#' @return Value of log-likelihood function.
#' @export
log_likelihood_anneal_func_da <- function(log_likelihood, log_like_approx, log_prior, cores = 1L){

  force(cores)
  force(log_likelihood)
  force(log_like_approx)
  force(log_prior)

  fun_cache <- list(
    log_likelihood = memoise::cache_memory(),
    log_like_approx = memoise::cache_memory()
  )

  get_key <- function(x, fun_name){
    papply(x, fun = function(x) fun_cache[[fun_name]]$digest(x))
  }

  check_hash <- function(x, fun_name){
    x_keys <- get_key(x, fun_name)
    sapply(x_keys, function(k) fun_cache[[fun_name]]$has_key(key = k))
  }

  add_hash <- function(x, values, fun_name){
    x_keys <- get_key(x, fun_name)
    invisible(
      mapply(function(k,v) fun_cache[[fun_name]]$set(key = k, value = v), k = x_keys, v = values)
    )
  }

  get_value <- function(x, fun_name, key_names = F){
    x_keys <- get_key(x, fun_name)
    vals <- sapply(x_keys, function(k) fun_cache[[fun_name]]$get(key = k))
    if(!key_names) names(vals) <- NULL
    vals
  }

  unique_x <- NULL
  evaluation_counts <- setNames( rep(list(0L),2), c("log_likelihood", "log_like_approx"))


  f <- function(particles, temp = 1, lh_trans = identity, type = "full_posterior", comp_time = F, ...){

    agg_comptime <- rep(0, times = num_particles(particles))
    agg_artiftime <- rep(0, times = num_particles(particles))

    if(type %in% c("full_posterior", "full_likelihood")){

      log_like_vals <- rep(NA_real_, times = num_particles(particles))

      has_hash <- check_hash(x = particles, fun_name = "log_likelihood")

      if(!all(has_hash)){
        new_vals <- papply(particles = particles[!has_hash, , drop = F], fun = log_likelihood, comp_time = comp_time, cores = cores, ...)
        add_hash(x = particles[!has_hash, , drop = F], values = new_vals, fun_name = "log_likelihood")
        unique_x <<- rbind(unique_x, particles[!has_hash, , drop = F])
        evaluation_counts$log_likelihood <<-  evaluation_counts$log_likelihood +
          sum(!has_hash)
        log_like_vals[!has_hash] <- new_vals

        if(comp_time){
          agg_comptime[!has_hash] <- agg_comptime[!has_hash] + attr(new_vals,"comptime")
          agg_artiftime[!has_hash] <- agg_artiftime[!has_hash] + attr(new_vals,"artiftime")
        }
      }

      if(any(has_hash)){
        cache_vals <- get_value(x = particles[has_hash, , drop = F], fun_name = "log_likelihood")
        log_like_vals[has_hash] <- cache_vals
      }

    }

    if(type %in% c("approx_posterior","approx_likelihood")){

      log_like_approx_vals <- rep(NA_real_, times = num_particles(particles))
      trans_x <- t(apply(X = particles, MARGIN = 1, FUN = lh_trans))
      has_hash <- check_hash(x = trans_x, fun_name = "log_like_approx")

      if(!all(has_hash)){

        new_vals <- papply(particles = trans_x[!has_hash, , drop = F], fun = log_like_approx, comp_time = comp_time, cores = cores, ...)

        add_hash(x = trans_x[!has_hash, , drop = F], values = new_vals, fun_name = "log_like_approx")
        evaluation_counts$log_like_approx <<- evaluation_counts$log_like_approx +
          sum(!has_hash)
        log_like_approx_vals[!has_hash] <- new_vals

        if(comp_time){
          agg_comptime[!has_hash] <- agg_comptime[!has_hash] + attr(new_vals,"comptime")
          agg_artiftime[!has_hash] <- agg_artiftime[!has_hash] + attr(new_vals,"artiftime")
        }
      }

      if(any(has_hash)){
        cache_vals <- get_value(x = trans_x[has_hash, , drop = F], fun_name = "log_like_approx")
        log_like_approx_vals[has_hash] <- cache_vals
      }

    }

    if(type %in% c("full_posterior","approx_posterior")){

      log_prior_vals <- papply(particles, fun = log_prior, cores = 1L, ...)

    }

    vals <- switch(type,
           approx_posterior = temp * log_like_approx_vals + log_prior_vals,
           approx_likelihood = temp * log_like_approx_vals,
           full_likelihood =  temp * log_like_vals,
           full_posterior = temp * log_like_vals + log_prior_vals,
           prior = log_prior_vals
      )

    if(comp_time){
      attr(vals, "comptime") <- agg_comptime
      attr(vals, "artiftime") <- agg_artiftime
    }

    return(vals)

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
#' @param max_approx_anneal Maximum value of temperature during approximate anneal.
#'
#' @return Value of log-likelihood function.
#' @export
log_approx_likelihood_anneal_func_da <- function(log_likelihood, log_like_approx, log_prior, max_approx_anneal){

  mem_log_likelihood <- memoise::memoise(log_likelihood)
  mem_log_likelihood_approx <-  memoise::memoise(log_like_approx)
  mem_log_prior <- memoise::memoise(log_prior)

  force(max_approx_anneal)

  unique_x <- NULL
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
           approx_posterior = temp * mem_log_likelihood_approx(lh_trans(x), ...) + ( 1 - temp ) * max_approx_anneal * mem_log_likelihood_approx(x, ...) + mem_log_prior(x, ...),
           approx_likelihood = temp * mem_log_likelihood_approx(lh_trans(x), ...),
           full_likelihood =  temp * mem_log_likelihood(x, ...),
           full_posterior = temp * mem_log_likelihood(x, ...) + ( 1 - temp ) * max_approx_anneal * mem_log_likelihood_approx(x, ...) + mem_log_prior(x, ...),
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
