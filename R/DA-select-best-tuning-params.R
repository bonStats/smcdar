#' Find tuning parameter with minimum computation cost under DAMH,
#'
#' @param min_T Minimum time taken to perform MH iteration for each parameter.
#' @param surrogate_acceptance Average surrogate acceptance rate for each tuning parameter.
#' @param surrogate_cost Average surrogate acceptance rate for each tuning parameter.
#' @param full_cost Average cost of using surrogate.
#'
#' @return Index for which parameter to use.
#' @export
min_da_mh_cost <- function(min_T, surrogate_acceptance, surrogate_cost, full_cost){

  stopifnot(length(min_T) == length(surrogate_acceptance),
            length(surrogate_cost) == 1,
            length(full_cost) == 1
  )

  total_costs <- min_T * (surrogate_acceptance * full_cost + surrogate_cost)

  which.min(total_costs)

}

#' Find tuning parameter with minimum computation cost under MH
#'
#' @param min_T Minimum time taken to perform MH iteration for each parameter.
#'
#' @return Index for which parameter to use.
#' @export
min_mh_cost <- function(min_T){

  which.min(min_T)

}
