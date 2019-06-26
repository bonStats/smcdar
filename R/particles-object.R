#' Create particles object.
#'
#' \code{particles} returns a particle with specified names.
#'
#' @param ... Names and types of object that is being used for each component in a particle.
#' @param weights (Un)normalised weights of particles being created.
#'
#' @export
#' @examples
#' num_particles <- 10
#' len <- 2
#' rv <- matrix(rnorm(num_particles * len), nrow = num_particles, ncol = len)
#' prts <- particles(beta = rv)
particles <- function(..., weights = 1){

  prt_list <- list(...)

  stopifnot(!is.null(names(prt_list)))

  stopifnot( components_have_same_particle_number(prt_list) )

  dim_list <- p_dim(prt_list)

  # attributes for components of particles
  comp_attr <- list(
    dim = dim_list,
    len = component_length(dim_list),
    pos = absolute_position(dim_list),
    nam = names(prt_list)
  )

  prts_mat <- do.call("cbind",
                lapply(prt_list, FUN = function(s) p_shape(s, ... = NULL))
              )

  stopifnot(length(weights) == 1 | length(weights) == nrow(prts_mat))

  unn_log_weights <- rep(0, times = nrow(prts_mat))
  unn_log_weights[] <- log(weights)

  class(prts_mat) <- c("particles","matrix")
  attr(prts_mat, "components") <- comp_attr
  attr(prts_mat, "log_weights") <- log_weights_normalise(unn_log_weights)
  return(prts_mat)

}

#' Print particles object.
#'
#' \code{particles} returns a particle with specified names.
#'
#' @param x particles to be printed.
#' @param ... Arguments to be passed to \code{str}.
#'
#' @export
#' @examples
#' x <- particles(beta = matrix(rnorm(20), nrow = 10, ncol = 2))
#' print(x)

print.particles <- function(x, ...){

  cat(num_particles(x), "particles with", length(x), "components:\n")
  compd <- p_components(x)

  print_len <- min(3, num_particles(x))

  for(comp_name in compd$nam){
    cat(comp_name,"[ length =", compd$len[comp_name],"]", x[1:print_len, compd$pos[comp_name]], "...\n")
  }

}

#' Particle components.
#'
#' Extract information about the components of a \code{particles} object.
#'
#' @param x particles object.
#' @param attr the particular attribute of the component to access.
#'
#' @export
p_components <- function(x, attr = NULL){

  if(is.null(attr)){

    attr(x, which = "components")

  } else {

    attr(x, which = "components")[[attr]]

  }

}

#' @export
length.particles <- function(x){

  length(p_components(x, attr = "nam"))

}


#' Number of particles per component.
#'
#' Get the number of particles for a particular object.
#'
#' @param x particles object.
#'
#' @export
num_particles <- function(x){

  nrow(x)

}

get_comp_details <- function(pcomp, i){

  if(is.numeric(i)){
    comp_i <- i
  } else {
    comp_i <- which( pcomp$nam ==  i)
  }

  stopifnot(length(comp_i) == 1L, comp_i > 0, comp_i <= length(pcomp$nam))

  list(name = pcomp$nam[comp_i],
       start = pcomp$pos[comp_i],
       stop = pcomp$pos[comp_i] + pcomp$len[comp_i] - 1,
       dim = pcomp$dim[[comp_i]]
       )

}

#' @export
`[[.particles` <- function(x, i, p_i, ...){

  stopifnot(length(i) == 1)

  pcomp <- p_components(x)

  compd <- get_comp_details(pcomp, i)

  if(missing(p_i)){

    p_i <- 1:num_particles(x)

  }

  if( length(compd$dim) == 1){

    x[p_i,compd$start]

  } else if(length(compd$dim) == 2){

    matrix(x[p_i,compd$start:compd$stop], nrow = length(p_i), ncol = compd$dim[2])

  } else if(length(compd$dim) >= 3){

    array(x[,compd$start:compd$stop], dim = c(length(p_i), compd$dim[-1]) )

  }

}

#' @export
`[[<-.particles` <- function(x, i, p_i, ..., value){

  stopifnot(length(i) == 1)

  pcomp <- p_components(x)

  compd <- get_comp_details(pcomp, i)

  if(missing(p_i)){

    p_i <- 1:num_particles(x)

  }

  x[p_i,compd$start:compd$stop] <- p_shape(value,...)

  return(x)

}

#' Get particle weights
#'
#' @param object Particles object.
#' @param log Return log weights?
#'
#' @return Weights
#' @export
#'
weights.particles <- function(object, log = FALSE){

  log_w <- attr(object, "log_weights")

  if(log){
    log_w
  } else {
    exp(log_w)
  }

}

#' Set  weights
#'
#' @param object Object.
#' @param ... Other parameters.
#' @param value New weights.
#'
#' @export
`weights<-` <- function(object, ..., value){

  UseMethod("weights<-")

}


#' Set particle weights
#'
#' @param object Particles object.
#' @param value New particle weights.
#' @param log Input is log weights?
#'
#' @export
#'
`weights<-.particles` <- function(object, log = FALSE, value){

  if(length(value)!= num_particles(object)) value <- rep(value, length.out = num_particles(object))

  if(log){
    log_w <- value
  } else {
    stopifnot( all(value >= 0) )
    log_w <- log(value)
  }

  attr(object, "log_weights") <- log_weights_normalise(log_w)

  return(object)

}

#' Effective sample size
#'
#' @param object
#'
#' @return ESS
#' @export
#'
ess <- function(object){

  UseMethod("ess")

}

#' @export
ess.particles <- function(object){

  w <- weights(object)

  # sum(w)^2 / sum(w^2)
  # sum(w) == 1 always
  1 / sum(w^2)

}

#' Select new set of particles from old.
#'
#' @param object Particle object
#' @param index Indices to use.
#' @param reweight Should new particles be assinged a uniform weight?
#'
#' @return Particles object.
#' @export
#'
select_particles <- function(object, index, reweight = T){

  object[] <- object[index,]

  if(reweight) weights(object) <- 1

  return(object)

}
