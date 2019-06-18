#### load packages ####
  library(smcdar)
  library(dplyr)
  library(tidyr)
  library(ggplot2)

#### simulate data ####
  # (birth prey, death prey/ birth pred, death pred)
  true_theta <- c(0.17,0.01,0.2)
  sim <- simulate_lotka_volterra_ctmc(
    true_theta,
    y1_min = 0,
    y2_min = 0,
    y1_max =100,
    y2_max=100,
    initial_y1 = 40,
    initial_y2 =25,
    total_time=50
    )

  tidy_sim <- sim %>% as.data.frame() %>%
    gather(pop,level,-time) %>%
    mutate(pop_name = ifelse(pop == "y1", "prey", "pred"))

  tidy_sim %>% ggplot() +
    geom_line(aes(x = time, y = level, colour = pop_name)) +
    theme_bw()

#### SMC ####

  # setup
  num_p <- 100



