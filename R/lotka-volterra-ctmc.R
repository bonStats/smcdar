#' Generator matrix for Lotka-Volterra continuous-time Markov chain
#'
#' @param theta Numeric vector of length 3.
#' @param y1_min Min for y1, integer, default is 0.
#' @param y1_max Max for y1, integer.
#' @param y2_min Min for y2, integer, default is 0.
#' @param y2_max Max for y2, integer.
#'
#' @return Sparse generator matrix
#'
#' @export
#'
generator_matrix_lotka_volterra_ctmc <- function(theta, y1_min = 0L, y1_max, y2_min = 0L, y2_max){

  stopifnot(length(theta) == 3)

  sparse_mat_list <- lotka_volterra_generator(theta = theta, y1_min = y1_min, y1_max = y1_max, y2_min = y2_min, y2_max = y2_max)

  Matrix::sparseMatrix(i = sparse_mat_list$rows,
                       j =  sparse_mat_list$cols,
                       x = sparse_mat_list$vals,
                       giveCsparse = FALSE
  )

}

#' Log-likelihood for Lotka-Volterra continuous-time Markov chain
#'
#' @param theta Parameters, numeric vector of length 3.
#' @param y1 Observed y1.
#' @param y2 Observed y2.
#' @param times Observed times.
#' @param y1_max Max for y1, integer.
#' @param y2_max Max for y2, integer.
#' @param ... To be passed to prob function.
#'
#' @return log-likelihood (numeric).
#'
#' @export
#'
log_lhood_lotka_volterra_ctmc <- function(theta, y1, y2, times, y1_max, y2_max, ...){

  Gmat <- generator_matrix_lotka_volterra_ctmc(theta = theta,
                                               y1_min = 0L, y1_max = y1_max,
                                               y2_min = 0L, y2_max = y2_max)

  flat_state_vec <- flattened_state_lotka_volterra(y1 = y1, y2 = y2,
                                                          y1_min = 0L, y2_min = 0L,
                                                          y1_max = y1_max, y2_max = y2_max)

  # parallelise if neccesary, if possible
  probs <- purrr::pmap_dbl(
    .l = list(
      flat_init_state = flat_state_vec[-length(flat_state_vec)],
      flat_final_state = flat_state_vec[-1],
      elapsed_time = diff(times)
        ),
    .f = prob_lotka_volterra_ctmc,
    generator_matrix = Gmat,
    ...
    )

  sum( log(probs) )

}

# obtain state of lotka volterra model in terms of flattened (2D -> 1D) variable
flattened_state_lotka_volterra <- function(y1, y2, y1_min, y1_max, y2_min, y2_max, mat){

  if(missing(y1)){
    y1 <- mat[,1]
    y2 <- mat[,2]
  }

  stopifnot(all( y1 >= y1_min),
            all( y2 >= y2_min),
            all( y1 <= y1_max),
            all( y2 <= y2_max)
            )

  return(
    (y1 - y1_min)*(y2_max - y2_min + 1) + (y2 - y2_min + 1)  # plus 1 for indices start at 1 in R
  )

}

unflattened_state_lotka_volterra <- function(z, y1_min, y1_max, y2_min, y2_max){

  y1 <- floor((z-1)/(y2_max-y2_min+1)) + y1_min
  y2 <- z - (y1 - y1_min)*(y2_max - y2_min + 1) + y2_min - 1

  return( as.matrix( data.frame(y1 = y1, y2 = y2) ))

}

prob_lotka_volterra_ctmc <- function(flat_init_state, flat_final_state, elapsed_time, generator_matrix, package = "expoRkit"){

  ffinal <- rep(0, times = nrow(generator_matrix))
  ffinal[flat_final_state] <- 1

  if(package == "expoRkit"){

    res <- expoRkit::expv(x = generator_matrix, v = ffinal, t = elapsed_time, Markov = T, transpose = F)[flat_init_state,]

  } else {

    # expm::expAtv could be sped up by using a sparse "v" since it is all 0s and one 1
    out <- expm::expAtv(A = generator_matrix, t = elapsed_time, v = ffinal)
    res <- out$eAtv[flat_init_state]
    if(abs(out$error) > 1e-03) warning("Numerical error above 0.001")

  }

  res

}

#' Simulate continuous-time Markov chain of Lotka-Volterra model
#'
#' @param theta Numeric vector of length 3.
#' @param y1_min Min for y1, integer, default is 0.
#' @param y1_max Max for y1, integer.
#' @param y2_min Min for y2, integer, default is 0.
#' @param y2_max Max for y2, integer.
#' @param initial_state_prob Probability matrix (y1, y2) for intial state of process.
#' @param initial_y1 Initial state of process y1.
#' @param initial_y2 Initial state of process y2.
#' @param total_time Total time that CTMC can run for.
#'
#' @return Trajectory (integer vector)
#'
#' @export
#'
simulate_lotka_volterra_ctmc <- function(theta, y1_min = 0, y2_min = 0, y1_max, y2_max, initial_y1, initial_y2, initial_state_prob, total_time){

  stopifnot(!missing(y1_max),
            !missing(y2_max),
            y1_min < y1_max,
            y2_min < y2_max
            )

  y1_len <- y1_max - y1_min + 1
  y2_len <- y2_max - y2_min + 1

  if(missing(initial_y1)){

    stopifnot(all( dim(initial_state_prob) == c(y1_max - y1_min, y2_max - y2_min) ))

    sample_y1 <- sample.int(n = y1_len, size = 1, prob = rowSums(initial_state_prob) ) - 1 + y1_min
    sample_y2 <- sample.int(n = y2_len, size = 1, prob = initial_state_prob[sample_y1,] )  - 1 + y2_min

    flat_initial_state <-
      flattened_state_lotka_volterra(y1 = sample_y1,
                                     y2 = sample_y2,
                                     y1_min = y1_min,
                                     y2_min = y2_min,
                                     y1_max = y1_max,
                                     y2_max = y2_max)

  } else {

    flat_initial_state <-
      flattened_state_lotka_volterra(y1 = initial_y1,
                                     y2 = initial_y2,
                                     y1_min = y1_min,
                                     y2_min = y2_min,
                                     y1_max = y1_max,
                                     y2_max = y2_max)

  }


  generator_matrix <- generator_matrix_lotka_volterra_ctmc(theta = theta,
                                                           y1_min = y1_min,
                                                           y2_min = y2_min,
                                                           y1_max = y1_max,
                                                           y2_max = y2_max)
  sim <- simulate_ctmc(generator_matrix = generator_matrix,
                       initial_state = flat_initial_state,
                       total_time = total_time)

  return(
    cbind(
      unflattened_state_lotka_volterra(
    sim$state, y1_min = y1_min, y1_max = y1_max,
    y2_min = y2_min, y2_max = y2_max),
    time = sim$time
    )
  )


}

