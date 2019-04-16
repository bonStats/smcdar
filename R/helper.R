components_have_same_particle_number <- function(prt_list){

  all( sapply(prt_list, check_p_num, N = p_num(prt_list[[1]])))

}

p_num <- function(x){

  p_dim(x)[1]

}

check_p_num <- function(x, N){

  p_num(x) == N

}

component_length <- function(dim_list){

   sapply(dim_list, FUN =
        function(x) prod(x) / x[1]
            )

}

absolute_position <- function(dim_list){

  c(1, cumsum(component_length(dim_list))[-1])

}
