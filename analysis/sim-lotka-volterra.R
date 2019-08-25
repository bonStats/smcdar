.libPaths(Sys.getenv("R_LIB_USER"))
outdir <- Sys.getenv("OUT_DIRECTORY")

devtools::install_github("bonStats/smcdar")

#### load packages ####

library(smcdar)
library(dplyr)
#library(tidyr)
#library(ggplot2)
library(future)
library(future.apply)

####

#### Settings ####

  ## Priors, likelihoods

  # prior: theta ~ logNormal(-1,0.5) or log(theta) ~ Normal(-1,0.5)
  log_prior <- function(log_theta){

    sum( dnorm(log_theta, mean = c(-2,-5,-2), sd = rep(0.5, 3), log = T) )

  }

  log_like <- function(log_theta, sim_data){
    log_lhood_lotka_volterra_ctmc(
      theta = exp(log_theta),
      y1 = sim_data[,"y1"], y2 = sim_data[,"y2"],
      y1_max = 100, y2_max = 100,
      times = sim_data[,"time"]
    )
  }

  log_like_approx <- function(log_theta, sim_data){
    log_lhood_lotka_volterra_lna(
      theta = exp(log_theta),
      y1 = sim_data[,"y1"], y2 = sim_data[,"y2"],
      times = sim_data[,"time"]
    )
  }

  # true_theta = (birth rate prey, death rate prey/ birth rate pred, death rate pred)
  g_pars <- list(N = 100, N_approx = 20, true_theta = c(0.17,0.01,0.2),
                 log_prior = log_prior, log_like = log_like, log_like_approx = log_like_approx
                 )

  sim_settings <- rep(list(list(f_pars = NULL, g_pars = g_pars)), 2)

  sim_settings[[1]]$f_pars <- list(
    num_p = 100,
    step_scale_set = c(0.01, 0.05, 0.1, 0.2),
    b_s_start = c(0,0,0,0,0,0),
    refresh_ejd_threshold = 0.01
  )

  sim_settings[[2]]$f_pars <- list(
    num_p = 200,
    step_scale_set = c(0.01, 0.05, 0.1, 0.2),
    b_s_start = c(0,0,0,0,0,0),
    refresh_ejd_threshold = 0.01
  )

####



#### Code for sim ####

  run_sim <- function(ss, verbose = F){

    ## simulate data

    sim <- simulate_lotka_volterra_ctmc(
      ss$g_pars$true_theta,
      y1_min = 0,
      y2_min = 0,
      y1_max = 100,
      y2_max = 100,
      initial_y1 = 40,
      initial_y2 = 25,
      total_time = 50
    )

    sim_data <- sim[1:min(ss$g_pars$N, nrow(sim)),]

    log_prior <- ss$g_pars$log_prior
    log_like <- function(log_theta) ss$g_pars$log_like(log_theta, sim_data = sim_data)
    log_like_approx <- function(log_theta) ss$g_pars$log_like_approx(log_theta, sim_data = sim_data[1:ss$g_pars$N_approx, ])

    ## more code
    source("analysis/lotka-volterra-ctmc-smc-da.R")

    ## run
    smc_da <- with(ss$f_pars, smc_lotka_volterra_da(
      use_da = T, use_approx = F,
      num_p = num_p,
      step_scale_set = step_scale_set,
      b_s_start = b_s_start,
      refresh_ejd_threshold = refresh_ejd_threshold,
      log_prior = log_prior,
      log_like = log_like,
      log_like_approx = log_like_approx,
      sim_data = sim_data,
      verbose = verbose
    )
    )
    cat(1)
    smc_full <- with(ss$f_pars, smc_lotka_volterra_da(
      use_da = F, use_approx = F,
      num_p = num_p,
      step_scale_set = step_scale_set,
      b_s_start = b_s_start,
      refresh_ejd_threshold = refresh_ejd_threshold,
      log_prior = log_prior,
      log_like = log_like,
      log_like_approx = log_like_approx,
      sim_data = sim_data
    )
    )
    cat(2)
    smc_approx <- with(ss$f_pars, smc_lotka_volterra_da(
      use_da = F, use_approx = T,
      num_p = num_p,
      step_scale_set = step_scale_set,
      b_s_start = b_s_start,
      refresh_ejd_threshold = refresh_ejd_threshold,
      log_prior = log_prior,
      log_like = log_like,
      log_like_approx = log_like_approx,
      sim_data = sim_data
    )
    )


    return(
      list(sim_settings = ss,
           sim_data = sim_data,
           smc_da = smc_da,
           smc_full = smc_full,
           smc_approx = smc_approx
      )
    )

  }

####

#### Run sim ####

  plan(multicore)

  for(i in seq_along(sim_settings)){

    res <- future_replicate(n = 1,
                                 run_sim(ss = sim_settings[[i]], verbose = T),
                                 simplify = F)

    #saveRDS(res, file.path(outdir, paste0("sim_res_",i,"_20190816.rds")))

  }


# optimisation likely too long with 100 particles? But already sampling...
# try simple mean model
# add function evaluations counts
