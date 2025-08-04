library(ggplot2)
library(dplyr)
library(tidyr)
library(parallel)
library(DEoptim)
library(purrr)
library(ggpubr)

# --- Simulation function ---
simulate_organoid <- function(mean_kill, sd_kill, p_evade) {
  duration <- 585  # minutes
  time_points <- seq(15, duration, by = 15)
  contact_rate <- 0.025
  total_apoptosis <- 0
  trace <- numeric(length(time_points))
  time_index <- 1
  for (t in 15:duration) {
    n_contacts <- rpois(1, contact_rate)
    kill_effects <- rnorm(n_contacts, mean = mean_kill, sd = sd_kill) * (1-p_evade)
    kill_effects <- pmax(0, kill_effects)  # no negative apoptosis
    total_apoptosis <- total_apoptosis + sum(kill_effects) 
    if (t %% 15 == 0) {
      trace[time_index] <- total_apoptosis 
      time_index <- time_index + 1
    }
  }
  return(trace)
}

simulate_population <- function(n, mean_kill, sd_kill, p_evade, label) {
  time_points <- seq(15, 585, by = 15)  # needed inside function

  cl <- makeCluster(8)
  clusterEvalQ(cl, library(MASS))

  clusterExport(cl, varlist = c("simulate_organoid"), envir = environment())

  # Pass required parameters as part of the function closure
  traces <- t(parSapply(cl, 1:n, function(i) {
    simulate_organoid(mean_kill = mean_kill,
                      sd_kill = sd_kill,
                      p_evade = p_evade)
  }))

  stopCluster(cl)

  df <- as.data.frame(traces)
  colnames(df) <- paste0("t", time_points)
  df$organoid_id <- paste0(1:n, "_", label)
  df$group <- label
  return(df)
}

calculate_mse <- function(obs,pred,return_df = FALSE) {
  obs_df <- obs %>% filter(treatment!='control' & !is.na(value)) %>% dplyr::select(time,value)
  pred_df <- pred %>% dplyr::select(time_min,apoptotic_area) %>%
    filter(time_min %in% unique(obs_df$time))
  colnames(pred_df) <- c('time','value')
  obs_counts <- obs_df %>%
    group_by(time) %>%
    summarise(n_obs = n(), .groups = 'drop')
  
  pred_downsampled <- pred_df %>%
    inner_join(obs_counts, by = "time") %>%
    group_by(time,n_obs) %>%
    nest() %>%
    mutate(sampled = map2(data, n_obs, ~ slice_sample(.x, n = .y))) %>%
    unnest(sampled) %>%
    ungroup() %>%
    dplyr::select(time, value)
  
  # Sort and bind obs and pred
  pred_downsampled_sort <- pred_downsampled %>% group_by(time) %>%
    arrange(time,value)
  colnames(pred_downsampled_sort) <-  c('pred_time','pred_apoptotic_area')
  obs_sort <- obs_df %>% group_by(time) %>%
    arrange(time,value)
  colnames(obs_sort) <-  c('time','observed_apoptotic_area')
  # Downsample to minimum number available at all time points to avoid timepoint biases
  obs_pred <- cbind(obs_sort,pred_downsampled_sort) %>% group_by(time) %>% slice_sample(n=3000)
  #colnames(obs_pred) <- c('time','observed_apoptotic_area','pred_time','pred_apoptotic_area')
  
  mean_squared_error = sum((obs_pred$observed_apoptotic_area - obs_pred$pred_apoptotic_area)^2)/nrow(obs_pred)
  #out_mse <- list('mse'=mean_squared_error,'mse_df'=obs_pred)
  return(mean_squared_error)
  if (return_df) {
    return(obs_pred)
  } else {
    return(mean_squared_error)
  }
}

sim_wrapper <- function(n_sensitive,
                        n_resistant,
                        mean_kill,
                        sd_kill,
                        p_evade_sensitive,
                        p_evade_resistant,
                        underestimation_factor = 0) {

  df_sensitive <- simulate_population(n_sensitive, mean_kill, sd_kill, p_evade = 0, "Sensitive")
  df_resistant <- simulate_population(n_resistant, mean_kill, sd_kill, p_evade_resistant, "Resistant")

  df_all <- bind_rows(df_sensitive, df_resistant)

  df_long <- df_all %>%
    pivot_longer(cols = starts_with("t"), names_to = "time", values_to = "apoptotic_area") %>%
    mutate(
      time_min = as.numeric(gsub("t", "", time)),
      apoptotic_area = pmax(0, apoptotic_area - underestimation_factor)  # Clipped at 0
    )
  return(df_long)
}

optim_wrapper <- function(par,
                          param_names,
                          apoptosis_obs,
                          n_sensitive,
                          n_resistant,
                          p_evade_sensitive) {

  # Convert parameter vector to named list
  params <- setNames(par, param_names)
  # Fill in defaults if a parameter isn't included
  mean_kill <- params[["mean_kill"]]
  sd_kill <- params[["sd_kill"]]

  p_evade_resistant <- if ("p_evade_resistant" %in% names(params)) {
    params[["p_evade_resistant"]]
  } else {
    0
  }

  underestimation_factor <- if ("underestimation_factor" %in% names(params)) {
    params[["underestimation_factor"]]
  } else {
    0
  }

  # Repeat simulation and average MSE
  mse_values <- replicate(3, {
    sim_df <- sim_wrapper(
      n_sensitive = n_sensitive,
      n_resistant = n_resistant,
      mean_kill = mean_kill,
      sd_kill = sd_kill,
      p_evade_sensitive = p_evade_sensitive,
      p_evade_resistant = p_evade_resistant,
      underestimation_factor = underestimation_factor
    )
    calculate_mse(apoptosis_obs, sim_df)
  })

  avg_mse <- mean(mse_values)
  #print(c(par, avg_mse))
  return(avg_mse)
}
  
######################
# Load observed data #
######################

apoptosis_obs <- read.delim('apoptosis_all_data.txt',sep='\t',header=T)

############################
# Define model parameters  #
############################
#set.seed(1234)
n_organoids <- 6000
resistant_fraction <- 0.073 
n_resistant <- round(n_organoids * resistant_fraction)
n_sensitive <- n_organoids - n_resistant
#duration <- 585  # minutes
#time_points <- seq(15, duration, by = 15)
#contact_rate <- 0.025 # mean NK contacts per min

mean_kill_lower <- 0
sd_kill_lower <- 0
underestimation_factor_lower <- 0

mean_kill_upper <- 4
sd_kill_upper <- 8
underestimation_factor_upper <- 10

p_evade_resistant_lower <- 0.5 
p_evade_resistant_upper <- 1
p_evade_sensitive <- 0



############################
# Optimize model 1 (Null)  #
############################


param_names_null <- c("mean_kill", "sd_kill")

sim_deopt_null <- DEoptim(
  fn = function(par) {
    optim_wrapper(
      par = par,
      param_names = param_names_null,
      apoptosis_obs = apoptosis_obs,
      n_sensitive = n_sensitive,
      n_resistant = n_resistant,
      p_evade_sensitive = p_evade_sensitive
    )
  },
  lower = c(mean_kill_lower, sd_kill_lower),
  upper = c(mean_kill_upper, sd_kill_upper),
  control = DEoptim.control(NP = 30, itermax = 200, trace = TRUE,parallelType = 0)
)

saveRDS(sim_deopt_null, file = "null_deoptim_full_output.rds")
cat("Optimization complete. Results saved to 'null_deoptim_full_output.rds'\n")

#######################################
# Optimize model 2 (Null with error)  #
#######################################


param_names_error <- c("mean_kill", "sd_kill", "underestimation_factor")

sim_deopt_error <- DEoptim(
  fn = function(par) {
    optim_wrapper(
      par = par,
      param_names = param_names_error,
      apoptosis_obs = apoptosis_obs,
      n_sensitive = n_sensitive,
      n_resistant = n_resistant,
      p_evade_sensitive = p_evade_sensitive
    )
  },
  lower = c(mean_kill_lower, sd_kill_lower, underestimation_factor_lower),
  upper = c(mean_kill_upper, sd_kill_upper, underestimation_factor_upper),
  control = DEoptim.control(NP = 30, itermax = 200, trace = TRUE,parallelType = 0)
)

saveRDS(sim_deopt_error, file = "null_error_deoptim_full_output.rds")
cat("Optimization complete. Results saved to 'null_error_deoptim_full_output.rds'\n")

###########################################
# Optimize model 3 (Null with resistant)  #
###########################################

param_names_resistant <- c("mean_kill", "sd_kill", "p_evade_resistant")

sim_deopt_resistant <- DEoptim(
  fn = function(par) optim_wrapper(
    par = par,
    param_names = param_names_resistant,
    apoptosis_obs = apoptosis_obs,
    n_sensitive = n_sensitive,
    n_resistant = n_resistant,
    p_evade_sensitive = p_evade_sensitive
  ),
  lower = c(mean_kill_lower, sd_kill_lower, p_evade_resistant_lower),
  upper = c(mean_kill_upper, sd_kill_upper, p_evade_resistant_upper),
  control = DEoptim.control(NP = 30, itermax = 200, trace = TRUE, parallelType = 0)
)

saveRDS(sim_deopt_resistant, file = "null_resistant_deoptim_full_output.rds")
cat("Optimization complete. Results saved to 'null_resistant_deoptim_full_output.rds'\n")

###########################################
# Resample from optimal parameters        #
###########################################

best_null <- sim_wrapper(n_sensitive,n_resistant,0.022038,7.779459, p_evade_sensitive,p_evade_resistant)
obs_pred_long_null <- calculate_mse(apoptosis_obs,best_null, return_df = TRUE) %>% select(time,observed_apoptotic_area,pred_apoptotic_area) %>%
  melt(id.vars = 'time')

best_null_mse_values <- replicate(100, {
  sim_df <- sim_wrapper(
    n_sensitive = n_sensitive,
    n_resistant = n_resistant,
    mean_kill = 0.022038,
    sd_kill = 7.779459,
    p_evade_sensitive = p_evade_sensitive,
    p_evade_resistant = p_evade_resistant
  )
  # Calculate MSE between simulated and observed data
  calculate_mse(apoptosis_obs, sim_df)
})

# Base with error, optimized params
best_error <- sim_wrapper_error(n_sensitive,n_resistant,0.058709,7.923023, p_evade_sensitive,p_evade_resistant,1.191737)
obs_pred_long_error <- calculate_mse(apoptosis_obs,best_error, return_df = TRUE) %>% select(time,observed_apoptotic_area,pred_apoptotic_area) %>%
  melt(id.vars = 'time')

best_error_mse_values <- replicate(100, {
  sim_df <- sim_wrapper_error(
    n_sensitive = n_sensitive,
    n_resistant = n_resistant,
    mean_kill = 0.058709,
    sd_kill = 7.923023,
    p_evade_sensitive = p_evade_sensitive,
    p_evade_resistant = p_evade_resistant,
    underestimation_factor = 1.191737
  )
  # Calculate MSE between simulated and observed data
  calculate_mse(apoptosis_obs, sim_df)
})

best_resistant <- sim_wrapper(n_sensitive,n_resistant, 0.353396,7.727807, p_evade_sensitive,p_evade_resistant = 0.946613)
obs_pred_long_resistant <- calculate_mse(apoptosis_obs,best_resistant, return_df = TRUE) %>% select(time,observed_apoptotic_area,pred_apoptotic_area) %>%
  melt(id.vars = 'time')

best_kill_res_mse_values <- replicate(100, {
  sim_df <- sim_wrapper(
    n_sensitive = n_sensitive,
    n_resistant = n_resistant,
    mean_kill = 0.353396,
    sd_kill = 7.727807,
    p_evade_sensitive = p_evade_sensitive,
    p_evade_resistant = 0.946613
  )
  # Calculate MSE between simulated and observed data
  calculate_mse(apoptosis_obs, sim_df)
})

opt_mse_samples <- data.frame(Null = best_kill_base_mse_values,'Null_with_error' = best_kill_base_error_mse_values, 'Null_with_resistant' = best_kill_res_mse_values)
opt_mse_samples_long <- melt(opt_mse_samples)
opt_mse_samples_long$variable <-gsub("_"," ",opt_mse_samples_long$variable)

my_comparisons <- list( c("Null", "Null with error"), c("Null with error", "Null with resistant"), c("Null", "Null with resistant") )
ggplot(opt_mse_samples_long,aes(variable,value)) + 
  geom_boxplot()+
  ylab('Mean squared error')+
  xlab('')+
  theme_bw()+
  theme(text=element_text(size=24))+
  stat_compare_means(comparisons = my_comparisons, label.y = c(9.25, 9.5, 9.75),method = "t.test")
#+
  #stat_compare_means(method = "t.test")
ggsave('MSE_Comparison_Null_versus_Resistant_Model.pdf')
ggsave('MSE_Comparison_Null_versus_Resistant_Model.png')


obs_pred_long_null$group <- 'null'
obs_pred_long_resistant$group <- 'null_resistant'
obs_pred_long_error$group <- 'null_error' 
obs_pred_long_all <- rbind(obs_pred_long_resistant,obs_pred_long_null,obs_pred_long_error)

ggplot(filter(obs_pred_long_all,time %in% c(60,300,585)),aes(x=value,fill=variable)) + geom_histogram(alpha=0.5,position = 'identity')+
  facet_wrap(group ~ time,scales = "free")+
  xlab('Apoptotic surface area (%)')+
  scale_fill_manual(values = c(observed_apoptotic_area = "red", pred_apoptotic_area = "blue"))+
  theme_bw()+
  theme(text=element_text(size=18),legend.position = 'bottom')
ggsave('Histogram_Comparison_Obs_Pred_Null_versus_Resistant_Model.pdf')
ggsave('Histogram_Comparison_Obs_Pred_Null_versus_Resistant_Model.png')

