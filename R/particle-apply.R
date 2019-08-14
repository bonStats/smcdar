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

    fun_t <- function(...){
      tout <- system.time({
        out <- fun(...)
      })

      attr(out, "t") <- tout["user.self"]
      return(out)
    }

    out <- lapply(split(particles[], row(particles[])), fun_t, ...)

    vals <- simplify2array(out)
    attr(vals, "comptime") <- unname(sapply(out, attr, which = "t"))

    vals

  }

}
