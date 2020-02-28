#' Apply function to all particles
#'
#' Function must be formatted to take rows from particle matrix.
#'
#' @param particles Particle object.
#' @param fun Function to apply.
#' @param comp_time Record time taken for comutation? Logical.
#' @param cores  Use multicore to evaluate? Uses parallel package.
#' @param ... Other arguments to pass to \code{fun}.
#'
#' @return Function evaluated at each row.
#' @export
#'
papply <- function(particles, fun, comp_time = F, cores = 1L, ...){

  stopifnot(is.matrix(particles))

  if(cores == 1L){

      res <- papply_1core(particles, fun = fun, comp_time = comp_time, ...)

  } else {

      grps <- allocate_to_cores(p_num = num_particles(particles), core_num = cores)
      split_particles <- lapply(1:cores, FUN = function(i) subset(particles, subset = grps == i))
      res_list <- parallel::mclapply(split_particles, FUN = papply_1core,
                         fun = fun,
                         ...,
                         comp_time = comp_time,
                         mc.cores = cores)

      are_errors <- sapply(res_list, class) == "try-error"

      if(any(are_errors)) stop("\n",sapply(res_list[are_errors], print))

      res <- c(res_list, recursive = T)

      if(comp_time){
        timing <- c(lapply(res_list, FUN = attr, which = "comptime"), recursive = T)
        artiftiming <- c(lapply(res_list, FUN = attr, which = "artiftime"), recursive = T)
        attr(res, "comptime") <- timing
        attr(res, "artiftime") <- artiftiming
      }

    }

  if(!is.null(dim(res))) warning("dim(res) != NULL: papply designed for functions that map each particle to a scalar.")

  return(res)

}

papply_1core <- function(particles, fun, comp_time = F, ...){

  if(!comp_time){

      res <- apply(particles, MARGIN = 1, FUN = fun, ...)

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

  }

  return(res)

}

allocate_to_cores <- function(p_num, core_num){

  if( p_num <= core_num){
    1:p_num
  } else {
    sort(rep(1:core_num, times = p_num %/% core_num + 1)[1:p_num])
  }

}
