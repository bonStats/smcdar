#' Time-derivative for LNA of Lotka-Volterra CTMC
#'
#' @param y Population vector.
#' @param theta Parameter vector.
#'
#' @return
#' @export
dy_dt_lotka_volterra_lna <- function(y, theta){

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

  return(ydot)

}


# log_lhood_lotka_volterra_lna <- function(theta, y, max_time, dt)
#
#   y1 = y(1,1); y2 = y(1,2);
#
#   [~,sol] = ode45(@(t,yi)lotka_volterra_lna(t,yi,theta),0:dt:dt,[y1;y2;0;0;0]);
#   loglike = 0.0;
#   for t = 2:(max_time+1)
#   Sigma_t = [sol(end,3) sol(end,4); sol(end,4) sol(end,5)];
#   if (any(diag(Sigma_t)<0))
#     loglike = -Inf;
#   return;
#   end
#   mu_t = [sol(end,1); sol(end,2)];
#   y_t = [y(t,1); y(t,2)];
#
#   if (y(t-1,1) == 0 && y(t-1,2) ~= 0) % prey extinct
#   loglike = loglike - 0.5/Sigma_t(2,2)*(y_t(2) - sol(end,2))^2;
#   elseif (y(t-1,1) ~= 0 && y(t-1,2) == 0) % predator extinct
#   loglike = loglike - 0.5/Sigma_t(1,1)*(y_t(1) - sol(end,1))^2;
#   elseif (y(t-1,1) == 0 && y(t-1,2) == 0) % both populations extinct.  Point mass transition density
#   return;
#   else
#     loglike = loglike -0.5*(y_t-mu_t)'*inv(Sigma_t)*(y_t-mu_t);
#   end
#
#   [~,sol] = ode45(@(t,yi)lotka_volterra_lna(t,yi,theta),0:dt:dt,[y_t;0;0;0]);
# end
#
