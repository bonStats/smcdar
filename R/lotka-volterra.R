#' Generator matrix for Lotka-Volterra continuous-time Markov chain
#'
#' @param theta Numeric vector of length 3.
#' @param y1_min Min for y1, integer, default is 0.
#' @param y1_max Max for y1, integer.
#' @param y2_min Min for y2, integer, default is 0.
#' @param y2_max Max for y2, integer.
#'
#' @return Sparse matrix
#'
#' @export
#'
generator_matrix_lotka_volterra_ctmc <- function(theta, y1_min = 0L, y1_max, y2_min = 0L, y2_max){

  stopifnot(length(theta) == 3)

  sparse_mat_list <- lotka_volterra_generator(theta = theta, y1_min = y1_min, y1_max = y1_max, y2_min = y2_min, y2_max = y2_max)

  Matrix::sparseMatrix(i = sparse_mat_list$rows + 1,
                       j =  sparse_mat_list$cols + 1,
                       x = sparse_mat_list$vals
  )

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
      flat_init_state = head(flat_state_vec, -1),
      flat_final_state = tail(flat_state_vec, -1),
      elapsed_time = diff(times)
        ),
    .f = prob_lotka_volterra_ctmc,
    generator_matrix = Gmat
    )

  sum( log(probs) )

}
