#' Time-derivative for LNA of Lotka-Volterra CTMC
#'
#' @param t time parameter (not used, but needed for ode solver).
#' @param y Population vector.
#' @param theta Parameter vector.
#'
#' @return Derivative vector for LNA approximation with components: (y1, y2, v11, v12, v22).
#' @export
d_dt_lotka_volterra_lna <- function(t, y, theta){

  stopifnot(
    length(y) == 5,
    length(theta) == 3
  )

  ydot <- rep(NA_real_, times = 5)

  h <- c(theta[1] * y[1],
         theta[2] * y[1] * y[2],
         theta[3] * y[2])

  S <- matrix(c(1, -1, 0, 0, 1, -1), byrow = T, nrow = 2)

  ydot[1:2] <- S %*% h

  Fmat <- matrix(
    c( theta[1] - theta[2]*y[2], -1*theta[2]*y[1],
       theta[2]*y[2], theta[2]*y[1] - theta[3] ),
    byrow = T,
    nrow = 2
  )

  V <- matrix(
    c( y[3], y[4],
       y[4], y[5]),
    byrow = T,
    nrow = 2
  )

  V_F_trans <- tcrossprod(V, Fmat)

  #F_V = Fmat %*% V

  S_H_S <- S %*% diag(h) %*% t(S) # can be optimised.

  ydot[3] <- V_F_trans[1,1] + S_H_S[1,1] + V_F_trans[1,1] # F_V[1,1] == V_F_trans[1,1]
  ydot[4] <- V_F_trans[1,2] + S_H_S[1,2] + V_F_trans[2,1] # F_V[1,2] == V_F_trans[2,1]
  ydot[5] <- V_F_trans[2,2] + S_H_S[2,2] + V_F_trans[2,2] # F_V[2,2] == V_F_trans[2,2]
   # are these equations correct?

  return(list(ydot))

}


#' Log-likelihood for Lotka-Volterra Linear Noise Approximation of CTMC
#'
#' @param theta Parameters, numeric vector of length 3.
#' @param y1 Observed y1.
#' @param y2 Observed y2.
#' @param times Observed times.
#' @param ... Arguments to pass to ode solver.
#'
#' @return log-likelihood (numeric)
#' @export
#'
log_lhood_lotka_volterra_lna <- function(theta, y1, y2, times, ...){

  # ode outside of loop?

  d_times <- diff(times)

  loglike <- 0

  for(t in 2:length(times) ){

    if(y1[t] == 0 & y2[t] == 0) return(loglike) # all extinct

    osol <- deSolve::ode(y = c(y1[t-1],y2[t-1],0,0,0),
                times = c(0, d_times[t-1]),
                func = d_dt_lotka_volterra_lna,
                parms = theta,
                method = "ode45", ...)

    Sigma_t <- matrix(c(osol[2,4],osol[2,5],
                       osol[2,5],osol[2,6]),
                     ncol = 2, byrow = T)

    if( det(Sigma_t) <= 0 ){
      #return(-Inf) # log-likelihood becomes -Inf
      next # ignore observation's contribution to log-likelihood
    }

    Mu_t <- osol[2,2:3]

    if( y1[t] == 0 & y2[t] != 0){ # prey extinct

      loglike <- loglike - 0.5*((y2[t] - Mu_t[2])^2) / Sigma_t[2,2]

    } else if( y1[t] != 0 & y2[t] == 0){ # predator extinct

      loglike <- loglike - 0.5*((y1[t] - Mu_t[1])^2) / Sigma_t[1,1]

    } else {

      chol_Sigma <- chol(Sigma_t)

     #t(c(y1[t], y2[t]) - Mu_t)%*% solve(Sigma_t) %*% (c(y1[t], y2[t]) - Mu_t)

      loglike <- loglike - 0.5 * crossprod(backsolve(chol_Sigma, c(y1[t], y2[t]) - Mu_t, transpose = T))[1]

    }

  }

  return(loglike)

}


