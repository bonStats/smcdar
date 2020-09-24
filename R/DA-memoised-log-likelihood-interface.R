#' Generic interface to annealed log-likeihood functions with memoisation: Likelihood annealing.
#'
#' Prior -> Posterior or
#' Approx posterior -> Posterior
#'
#' @param log_likelihood (Full) log-likelihood function.
#' @param log_like_approx Approximate log-likelihood function.
#' @param cwise_log_like_approx Approximate log-likelihood function.
#' @param log_prior Log-prior function.
#' @param cores Cores to use.
#' @param start_from_approx Is the intial distribution an approximate model?
#' @param max_approx_anneal Maximum value of temperature reached during approximate anneal.
#' @param count_start Likelihood evaluation count. List with integer elements "log_likelihood" and "log_like_approx".
#'
#' @return Value of log-likelihood function.
#' @export
log_likelihood_anneal_func_da <- function(log_likelihood, log_like_approx, cwise_log_like_approx, log_prior, cores = 1L, start_from_approx = F, max_approx_anneal = 1, count_start = NULL){

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
  if(is.null(count_start)){
    evaluation_counts <- setNames( rep(list(0L),2), c("log_likelihood", "log_like_approx"))
  } else {
    evaluation_counts <- count_start
  }

  f <- function(particles, temp = 1, lh_trans = identity, type = "full_posterior", comp_time = F, weights = NULL, ...){

    agg_comptime <- rep(0, times = num_particles(particles))
    agg_artiftime <- rep(0, times = num_particles(particles))

    if(type %in% c("full_posterior", "full_likelihood")){

      log_like_vals <- rep(NA_real_, times = num_particles(particles))

      has_hash <- check_hash(x = particles, fun_name = "log_likelihood")

      if(!all(has_hash)){
        new_vals_ll <- papply(particles = particles[!has_hash, , drop = F], fun = log_likelihood, comp_time = comp_time, cores = cores, ...)
        add_hash(x = particles[!has_hash, , drop = F], values = new_vals_ll, fun_name = "log_likelihood")
        unique_x <<- rbind(unique_x, particles[!has_hash, , drop = F])
        evaluation_counts$log_likelihood <<-  evaluation_counts$log_likelihood +
          sum(!has_hash)
        log_like_vals[!has_hash] <- new_vals_ll

        if(comp_time){
          agg_comptime[!has_hash] <- agg_comptime[!has_hash] + attr(new_vals_ll,"comptime")
          agg_artiftime[!has_hash] <- agg_artiftime[!has_hash] + attr(new_vals_ll,"artiftime")
        }
      }

      if(any(has_hash)){
        cache_vals_ll <- get_value(x = particles[has_hash, , drop = F], fun_name = "log_likelihood")
        log_like_vals[has_hash] <- cache_vals_ll
      }

    }

    if(type %in% c("approx_posterior","approx_likelihood") & is.null(weights)){

      log_like_approx_vals <- rep(NA_real_, times = num_particles(particles))
      trans_x <- t(apply(X = particles, MARGIN = 1, FUN = lh_trans))
      has_hash <- check_hash(x = trans_x, fun_name = "log_like_approx")

      if(!all(has_hash)){

        new_vals_app <- papply(particles = trans_x[!has_hash, , drop = F], fun = log_like_approx, comp_time = comp_time, cores = cores, weights = weights, save_comps = F, ...)

        add_hash(x = trans_x[!has_hash, , drop = F], values = new_vals_app, fun_name = "log_like_approx")
        evaluation_counts$log_like_approx <<- evaluation_counts$log_like_approx +
          sum(!has_hash)
        log_like_approx_vals[!has_hash] <- new_vals_app

        if(comp_time){
          agg_comptime[!has_hash] <- agg_comptime[!has_hash] + attr(new_vals_app,"comptime")
          agg_artiftime[!has_hash] <- agg_artiftime[!has_hash] + attr(new_vals_app,"artiftime")
        }

      }

      if(any(has_hash)){
        cache_vals_app <- get_value(x = trans_x[has_hash, , drop = F], fun_name = "log_like_approx")
        log_like_approx_vals[has_hash] <- cache_vals_app
      }

    }

    if(type %in% c("approx_posterior","approx_likelihood") & !is.null(weights)){

      stopifnot(is.function(cwise_log_like_approx))

      log_like_approx_vals <- rep(NA_real_, times = num_particles(particles))
      trans_x <- t(apply(X = particles, MARGIN = 1, FUN = lh_trans))

      cheap_log_like_approx <- function(x) cwise_log_like_approx(x, calc_comps = weights > 0)

      log_like_approx_vals <- papply(particles = trans_x, fun = cheap_log_like_approx, comp_time = comp_time, cores = cores, weights = weights, save_comps = F, ...)

      evaluation_counts$log_like_approx <<- evaluation_counts$log_like_approx +
        length(log_like_approx_vals)

      if(comp_time){
        agg_comptime <- agg_comptime + attr(log_like_approx_vals,"comptime")
        agg_artiftime <- agg_artiftime + attr(log_like_approx_vals,"artiftime")
      }

    }

    if(type %in% c("full_posterior","approx_posterior") ){

      log_prior_vals <- papply(particles, fun = log_prior, cores = 1L, ...)

    }

    if(!start_from_approx){

      vals <- switch(type,
             approx_posterior = temp * log_like_approx_vals + log_prior_vals,
             approx_likelihood = temp * log_like_approx_vals,
             full_likelihood =  temp * log_like_vals,
             full_posterior = temp * log_like_vals + log_prior_vals,
             prior = log_prior_vals
        )

    } else {

      log_like_approx_vals_no_trans <- rep(NA_real_, times = num_particles(particles))
      has_hash <- check_hash(x = particles, fun_name = "log_like_approx")

      if(!all(has_hash)){

        new_vals_app_nt <- papply(particles = particles[!has_hash, , drop = F], fun = log_like_approx, comp_time = comp_time, cores = cores, ...)

        add_hash(x = particles[!has_hash, , drop = F], values = new_vals_app_nt, fun_name = "log_like_approx")
        evaluation_counts$log_like_approx <<- evaluation_counts$log_like_approx +
          sum(!has_hash)
        log_like_approx_vals_no_trans[!has_hash] <- new_vals_app_nt

        if(comp_time){
          agg_comptime[!has_hash] <- agg_comptime[!has_hash] + attr(new_vals_app_nt,"comptime")
          agg_artiftime[!has_hash] <- agg_artiftime[!has_hash] + attr(new_vals_app_nt,"artiftime")
        }
      }

      if(any(has_hash)){
        cache_vals_app_nt <- get_value(x = particles[has_hash, , drop = F], fun_name = "log_like_approx")
        log_like_approx_vals_no_trans[has_hash] <- cache_vals_app_nt
      }

      vals <- switch(type,
                    approx_posterior = temp * log_like_approx_vals + ( 1 - temp ) * max_approx_anneal * log_like_approx_vals_no_trans + log_prior_vals,
                    approx_likelihood = temp * log_like_approx_vals,
                    full_likelihood =  temp * log_like_vals,
                    full_posterior = temp * log_like_vals + ( 1 - temp ) * max_approx_anneal * log_like_approx_vals_no_trans + log_prior_vals,
                    prior =  log_prior_vals
      )
    }

    if(comp_time){
      attr(vals, "comptime") <- agg_comptime
      attr(vals, "artiftime") <- agg_artiftime
    }


    return(vals)

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
