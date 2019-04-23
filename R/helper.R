
is_particles_obj <- function(particles){

  "particles" %in% class(particles)

}



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

  comp_len <- component_length(dim_list)
  pos <- cumsum(comp_len) - comp_len + 1
  names(pos) <- names(dim_list)
  pos

}

is_upper_triangular_matrix <- function(mat, eps = 1e-10){

  sum( abs( mat[lower.tri(mat)]) ) < eps

}
