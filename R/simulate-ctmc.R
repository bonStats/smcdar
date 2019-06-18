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

  in_row_i <- which( generator_matrix@i == as.integer(i - 1) ) # -1 = zero index

  # remove diagonal (and any others non-distiguishable from zero)
  in_row_i <- in_row_i[generator_matrix@x[in_row_i] > 0]

  # if only one option
  if(length(in_row_i) == 1) return(generator_matrix@j[in_row_i] + 1)  # +1 = zero index

  probs_in_row_i <- generator_matrix@x[in_row_i] / sum_rates

  rand_row_j <- sample(x = generator_matrix@j[in_row_i] + 1, size = 1, prob = probs_in_row_i)  # +1 = zero index

  return(rand_row_j)

}
