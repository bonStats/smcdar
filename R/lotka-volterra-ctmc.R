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

  Matrix::sparseMatrix(i = sparse_mat_list$rows + 1,
                       j =  sparse_mat_list$cols + 1,
                       x = sparse_mat_list$vals,
                       giveCsparse = FALSE
  )

}

#' Log-likelihood for Lotka-Volterra continuous-time Markov chain
#'
#' @param theta Numeric vector of length 3.
#' @param y1 Observed y1.
#' @param y2 Observed y2.
#' @param times Observed times.
#' @param y1_max Max for y1, integer.
#' @param y2_max Max for y2, integer.
#'
#' @return log-likelihood (numeric)
#'
#' @export
#'
log_lhood_lotka_volterra_ctmc <- function(theta, y1, y2, times, y1_max, y2_max){

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
    generator_matrix = Gmat
    )

  sum( log(probs) )

}

# obtain state of lotka volterra model in terms of flattened (2D -> 1D) variable
flattened_state_lotka_volterra <- function(y1, y2, y1_min, y1_max, y2_min, y2_max){

  return(
    (y1 - y1_min)*(y2_max - y2_min + 1) + (y2 - y2_min + 1) + 1 # plus 1 for indices start at 1 in R
  )

}

prob_lotka_volterra_ctmc <- function(flat_init_state, flat_final_state, elapsed_time, generator_matrix){

  ffinal <- rep(0, times = nrow(generator_matrix))
  ffinal[flat_final_state] <- 1

  # expm::expAtv could be sped up by using a sparse "v" since it is all 0s and one 1
  res <- expm::expAtv(A = generator_matrix, t = elapsed_time, v = ffinal)

  if(abs(res$error) > 1e-03) warning("Numerical error above 0.001")

  res$eAtv[flat_init_state]


}

#' Simulate continuous-time Markov chain
#'
#' @param generator_matrix Square generator matrix (possibly sparse).
#' @param initial_state_prob Probability vector for intial state of process.
#' @param initial_state Initial state of process.
#' @param total_time Total time that CTMC can run for.
#'
#' @return Trajectory (integer vector)
#'
#' @export
#'
simulate_ctmc <- function(generator_matrix, initial_state_prob, initial_state, total_time){

  # set initial state
  if(missing(initial_state)){
    if(missing(initial_state_prob)) stop("At least one of initial_state_prob and initial_state must be specified.")
    states <- sample(x = 1:length(initial_state_prob), prob = initial_state_prob, size = 1)
  } else {
    states <- as.integer(initial_state)
  }

  times <- 0 # start of time vec

  curr_state <- states
  curr_time <- times
  curr_time_rate <- -generator_matrix[curr_state, curr_state]

  while(curr_time < total_time){

    if(curr_time_rate > 0){

      times[length(times) + 1] <- curr_time + stats::rexp(n = 1, rate = -generator_matrix[curr_state, curr_state])

      states[length(states) + 1] <- sample_row_of_sparse_generator(generator_matrix = generator_matrix, i = curr_state)

      curr_state <- states[length(states)]
      curr_time <- times[length(times)]
      curr_time_rate <- -generator_matrix[curr_state, curr_state]

    } else {

      # reached holding...
      times[length(times) + 1] <- curr_time <- total_time
      states[length(states) + 1] <- curr_state

    }

  }

  return(
    list(
      time = times,
      state = states
    )
  )

}

sample_row_of_sparse_generator <-  function(generator_matrix, i){

  stopifnot(methods::is(generator_matrix, 'TsparseMatrix'))

  sum_rates <- - generator_matrix[i,i]

  if( abs(sum_rates) < 1e-10 ) return(i)

  in_row_i <- which( generator_matrix@i == ( as.integer(i) - 1L ) )  # 0-index

  # remove diagonal (and any others non-distiguishable from zero)
  in_row_i <- in_row_i[generator_matrix@x[in_row_i] > 0]

  # if only one option
  if(length(in_row_i) == 1) return(generator_matrix@j[in_row_i])

  probs_in_row_i <- generator_matrix@x[in_row_i] / sum_rates

  rand_row_j <- sample(x = generator_matrix@j[in_row_i], size = 1, prob = probs_in_row_i)

  return(rand_row_j)

}