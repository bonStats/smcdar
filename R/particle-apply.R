#' Apply function to all particles
#'
#' Function must be formatted to take rows from particle matrix.
#'
#' @param particles Particle object.
#' @param fun Function to apply.
#' @param comp_time Record time taken for comutation? Logical.
#' @param ... Other arguments to pass to \code{fun}.
#'
#' @return Function evaluated at each row.
#' @export
#'
papply <- function(particles, fun, comp_time = F, ...){

  if(!comp_time){

    apply(particles, MARGIN = 1, FUN = fun, ...)

  } else {

    fun_time <- function(x, ...){
      tic <- Sys.time()
      res <- fun(x, ...)
      toc <- Sys.time()
      art <- ifelse(is.null(attr(res, "comptime")), 0, attr(res, "comptime"))
      tim <- toc - tic + art
      list(res = res, time = tim, artifical_time = art)
    }

    res_list <- apply(particles, MARGIN = 1, FUN = fun_time, ...)

    res <- simplify2array(lapply(res_list, getElement, name = "res"))
    timing <- sapply(res_list, getElement, name = "time")
    artiftiming <- sapply(res_list, getElement, name = "artifical_time")

    attr(res, "comptime") <- timing
    attr(res, "artiftime") <- artiftiming

    return(res)
  }

}

#' papply2 <- function(particles, fun, comp_time = F, groups, ...){
#'
#'   if(!comp_time){
#'
#'     apply(particles, MARGIN = 1, FUN = fun, ...)
#'
#'   } else {
#'
#'     if(missing(groups)) groups <- rep(1, num_particles(particles))
#'
#'     res <- vector(mode = "list", length = num_particles(particles))
#'     ugroups <- unique(groups)
#'     timing <- setNames(rep(NA_real_, length(ugroups)), ugroups)
#'     len <- rep(NA_integer_, length(ugroups))
#'
#'     for(tgp in ugroups){
#'       tic <- Sys.time()
#'       res[tgp == groups] <- apply(particles[tgp == groups, , drop = F], MARGIN = 1, FUN = fun, ...)
#'       len <- sum(tgp == groups)
#'       toc <- Sys.time()
#'       timing[tgp] <- toc - tic
#'     }
#'
#'     res <- simplify2array(res)
#'     attr(res, "comptime") <- timing / len
#'     return(res)
#'   }
#'
#' }
