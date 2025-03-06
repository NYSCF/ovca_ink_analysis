library(brms)
library(ggplot2)
library(ggpubr)
library(tidybayes)
library(marginaleffects)
library(reshape2)
library(dplyr)
library(tidyr)
library(stringr)
library(broom.mixed)
library(forcats)
library(scales)
library(investr)
library(report)
library(cowplot)

convert_time_column_to_minutes <- function(df, time_column) {
  convert_to_minutes <- function(time_string) {
    hours <- str_extract(time_string, "\\d+(?=h)")  # Extract hours
    minutes <- str_extract(time_string, "\\d+(?=min)")  # Extract minutes
    # Convert to numeric and handle NA
    hours <- ifelse(is.na(hours), 0, as.numeric(hours))
    minutes <- ifelse(is.na(minutes), 0, as.numeric(minutes))
    # Calculate total minutes
    total_minutes <- hours * 60 + minutes
    return(total_minutes)
  }
  # Apply the conversion function to the specified column
  df <- df %>%
    mutate(!!time_column := sapply(!!sym(time_column), convert_to_minutes))
  return(df)
}

# Define the logistic function
logistic_function <- function(time, lmax, k, x0) {
  lmax / (1 + exp(-k * (time - x0)))
}
# lmax: asymptote (e.g. highest value of apoptotic organoids)
# k: slope
# x0: inflection point

# Define the Hill model function
hill_model <- function(x, a, K, n) {
  a * (x^n) / (K^n + x^n)
}
# Another definition of the the Hill model function with other param names
# and the offset c instead of Vmax 
hill_function <- function(x, c, Vmax, Km, alpha) {
  (1-c) * (Vmax * x^alpha / (x^alpha + Km^alpha))
}
#x: Value of the effector.
#Vmax	: Maximal value of the function.
#Km:	 Value of the effector at which the value of the function is half-maximal.
# alpha: Hill coefficient defining the steepness of the function. Default is 1.

generalized_logistic_function <- function(time, y0,mumax,K,beta) {
  K*(1-exp(-beta * mumax * time)*(1-(y0/K)^-beta))^(-1/beta)
}
#https://rdrr.io/cran/growthrates/man/grow_richards.html
# Richards, F. J. (1959) A Flexible Growth Function for Empirical Use. Journal of Experimental Botany 10 (2): 290–300.
# Tsoularis, A. (2001) Analysis of Logistic Growth Models. Res. Lett. Inf. Math. Sci, (2001) 2, 23–46.
#y0 initial value of abundance,
#mumax maximum growth rate (note different interpretation compared to exponential growth),
#K carrying capacity (max. total concentration of cells),
#beta shape parameter determining the curvature.

# Define the Gompertz function
gompertz_function <- function(time, a, b, c) {
  a * exp(-b * exp(-c * time))
}

# a: carrying capacity (asymptote)
# b: displacement along x-axis
# c growth rate

gompertz_upper_lower_function <- function(time, L, U, b, c) {
  L + (U - L) * exp(-b * exp(-c * time))
}
#b is a scaling parameter,
#c controls the growth rate,
#t is time.
#L is the lower asymptote (left asymptote)
#U is the upper asymptote.


# Modified Gompertz function
modified_gompertz_function <- function(time, y0,mumax,K,lambda) {
  y0*(K/y0)^(exp(-exp((exp(1)*mumax*(lambda - time))/log(K/y0)+1)))
}
#y0 initial value of abundance,
#mumax maximum growth rate (1/time),
#K maximum abundance (carrying capacity),
#lambda time of lag phase of the 3 parameter Gompertz model .

simulate_for_pred <- function(t0,tn,sim_n,patient_var,fit_model) {
  new_data <- data.frame(time = seq(t0, tn, length.out = sim_n))
  new_data$patient <- patient_var
  # Sample from the posterior distribution to generate predictions
  posterior_samples <- posterior_predict(fit_model, newdata = new_data, draws = sim_n)
  # Calculate the mean and credible intervals
  predicted_means <- apply(posterior_samples, 2, mean)
  lower_ci <- apply(posterior_samples, 2, quantile, probs = 0.05)
  upper_ci <- apply(posterior_samples, 2, quantile, probs = 0.95)
  predictions <- data.frame(time = new_data$time, 
                            value = predicted_means, 
                            lower_ci = lower_ci, 
                            upper_ci = upper_ci)
  predictions <- predictions %>% mutate(lower_ci = ifelse(lower_ci<0,0,lower_ci),
                                       upper_ci = ifelse(upper_ci>1,1,upper_ci))
  return(predictions)
}

simulate_for_epred <- function(t0,tn,sim_n,patient_var,fit_model) {
  new_data <- data.frame(time = seq(t0, tn, length.out = sim_n))
  new_data$patient <- patient_var
  # Sample from the posterior distribution to generate predictions
  posterior_samples <- posterior_epred(fit_model, newdata = new_data, draws = sim_n)
  # Calculate the mean and credible intervals
  predicted_means <- apply(posterior_samples, 2, mean)
  lower_ci <- apply(posterior_samples, 2, quantile, probs = 0.05)
  upper_ci <- apply(posterior_samples, 2, quantile, probs = 0.95)
  predictions <- data.frame(time = new_data$time, 
                            value = predicted_means, 
                            lower_ci = lower_ci, 
                            upper_ci = upper_ci)
  predictions <- predictions %>% mutate(lower_ci = ifelse(lower_ci<0,0,lower_ci),
                                        upper_ci = ifelse(upper_ci>1,1,upper_ci))
  return(predictions)
}

simulate_for_nls_pred <- function(t0,tn,fit_model) {
  new.data <- data.frame(time=seq(t0, tn, by = 1))
  interval <- as_tibble(predFit(fit_model, newdata = new.data, interval = "confidence", level= 0.99)) %>% 
    mutate(time = new.data$time)
}

#################
### Load data ###
#################

## Load data ##
pd1 <- read.delim(gzfile("../data/csv/organoid_apoptosis_pd01.csv.gz"),header = TRUE,sep=',')
pd2 <- read.delim(gzfile("../data/csv/organoid_apoptosis_pd02.csv.gz"),header = TRUE,sep=',')
pd1 <- melt(pd1)
pd2 <- melt(pd2)
pd1$group <- "p1"
pd2$group <- "p2"
pd <- rbind(pd1,pd2)

pd$variable <- gsub("Control","control",pd$variable)
pd$variable <- gsub("control.","control_",pd$variable)
pd$variable <- gsub("NK.","NK_",pd$variable)
pd$variable <- gsub(".exp","_exp",pd$variable)

# split on underscore
pd <- separate_wider_delim(pd, cols = variable, delim = "_", names = c("treatment","time", "expt"))
pd$time <- gsub("\\.","",pd$time)
pd$expt <- gsub("\\.","",pd$expt)

pd <- convert_time_column_to_minutes(pd, "time")
pd$time <- as.character(pd$time)
s=0.5
colnames(pd) <- c("treatment","time","expt","value","patient","cluster","k1","k2")

#######################################################
### Prepare apoptotic proportion data for modelling ###
#######################################################

# Hardcoded threshold from previous analysis
resistant_threshold <- 8.5812

pd_class <- pd[pd$treatment=="NK",] %>% mutate(component=ifelse(value<=resistant_threshold,"resistant","susceptible"))
pd_all <- pd %>% mutate(component=ifelse(value<=resistant_threshold,"resistant","susceptible"))%>%
  mutate(component=ifelse(treatment!="NK","Control",component))

pd_class$time <- as.numeric(pd_class$time)
pd_all$time <- as.numeric(pd_all$time)

# Compare proportions of apoptotic and non-apoptotic organoids over time
pd_theta_stat <- pd_all[complete.cases(pd_all),] %>%
  mutate(component=ifelse(value<=resistant_threshold,"resistant","susceptible")) %>%
  group_by(treatment,time,expt,patient) %>%
  summarize(n_apoptotic = sum(component == "susceptible"),
            n_nonapoptotic = sum(component == "resistant")) %>%
  mutate(n_apoptotic_perc = (n_apoptotic/(n_apoptotic+n_nonapoptotic))*100,
         n_nonapoptotic_perc = (n_nonapoptotic/(n_apoptotic+n_nonapoptotic))*100)

pd_theta_stat_err <- pd_theta_stat  %>%
  group_by(treatment,time,patient) %>%
  summarize(mean_apoptotic_perc = mean(n_apoptotic_perc),
            mean_nonapoptotic_perc = mean(n_nonapoptotic_perc),
            max_apoptotic_perc = max(n_apoptotic_perc),
            min_apoptotic_perc = min(n_apoptotic_perc),
            max_nonapoptotic_perc = max(n_nonapoptotic_perc),
            min_nonapoptotic_perc = min(n_nonapoptotic_perc))

pd_theta_stat$patient <- gsub("p","PDO",pd_theta_stat$patient)
pd_theta_stat <- melt(pd_theta_stat,id.vars = c("treatment",
                                                "time",
                                                "expt",
                                                "patient",
                                                "n_apoptotic",        
                                                "n_nonapoptotic"))
pd_theta_stat$patient <- gsub("p","PDO",pd_theta_stat$patient)
pd_theta_stat_err$patient <- gsub("p","PDO",pd_theta_stat_err$patient)

############################################
### Fit models by Ordinary Least Squares ###
############################################

# Fit logistic function with offset to proportion data

# Prepare data for fitting
prop_data <- pd_theta_stat[pd_theta_stat$treatment=="NK" & pd_theta_stat$variable=="n_apoptotic_perc",]
prop_data$value <- prop_data$value / 100

# Fit logistic model by ordinary least squares to get information on useful priors
logistic_fit_ols <- nls(value ~ SSlogis(time, Asym, xmid, scal), data = prop_data)

# Note that while Asym = lmax and xmid=x0, scal!=k 
# scal in SSlogis is related to k but is its inverse, so large scal corresponds to small k

fit_hill <- nls(value ~ hill_model(time, a, K, n), data = prop_data,
                start = list(a = 1, K = 100, n = 1))

# Fit Gompertz
fit_gompertz <- nls(value ~ gompertz_function(time, a, b, c), data = prop_data,
                start = list(a = 0.92, b = 10, c = 0.0099))

# Create predictions for plotting
prop_data$predicted_hill <- predict(fit_hill)
prop_data$predicted_logistic <- predict(logistic_fit_ols)
prop_data$predicted_gompertz <- predict(fit_gompertz)

# Create confidence intervals for plotting
gompertz_conf <- simulate_for_nls_pred(0,600,fit_gompertz) 
logistic_conf <- simulate_for_nls_pred(0,600,logistic_fit_ols) 
hill_conf <- simulate_for_nls_pred(0,600,fit_hill) 

# Plot the original data and the fitted Hill curve
colors <- c("Hill" = "blue", "Logistic" = "orange", "Gompertz" = "darkgreen")
ggplot(prop_data, aes(x = time, y = value*100)) +
  geom_point(alpha = 0.5) +
  geom_line(aes(y = predicted_hill*100,color = "Hill"),  linewidth = 1.2) +
  geom_line(aes(y = predicted_logistic*100,color = "Logistic"),  linewidth = 1.2) +
  geom_line(aes(y = predicted_gompertz*100,color = "Gompertz"),  linewidth = 1.2) +
  geom_ribbon(data=logistic_conf, aes(x=time, ymin=lwr*100, ymax=upr*100), alpha=0.5, inherit.aes=F, fill="orange")+
  geom_ribbon(data=hill_conf, aes(x=time, ymin=lwr*100, ymax=upr*100), alpha=0.5, inherit.aes=F, fill="blue")+
  geom_ribbon(data=gompertz_conf, aes(x=time, ymin=lwr*100, ymax=upr*100), alpha=0.5, inherit.aes=F, fill="darkgreen")+
  labs(title = "Comparison of two sigmoid fits",
       x = "Time",
       y = "Percentage apoptotic organoids", 
       color = "Sigmoid type") +
  scale_color_manual(values = colors)+
  theme_minimal()+
  theme(text = element_text(size=18))
ggsave("../figures/NLS_model_comparison_ordinary_least_squares.pdf")


#################################################
### Fit models by Bayesian regression in brms ###
#################################################

## Now use a beta error model
#fit_sigmoid_c_beta <- brm(
#  bf(value ~ (1-c) * (1 / (1 + exp(-k * (time - x0)))),c ~1 ,k ~1, x0 ~1,nl=TRUE),
#  data = prop_data,
#  prior = c(prior(gamma(0.01, 0.01), class = "phi"),
#            prior(normal(0.1, 0.3), class = "b", nlpar = "c",ub=1,lb=0),   # Prior for the asymptote (lmax)
#            prior(normal(0, 0.5), class = "b", nlpar = "k"),   # Prior for the slope (k)
#            prior(normal(162, 10), class = "b", nlpar = "x0") # Prior for the midpoint (x0)
#  ), 
#  family = brmsfamily(family="Beta", link="identity", link_phi="identity"),
#  #init = 0, 
#  chains = 4,
#  cores = 4,
#  threads = 4,
#  iter = 4000,
#  warmup = 2000,
#  control = list(adapt_delta = 0.95,max_treedepth = 12)
#)
#
##summary(fit_sigmoid_c_beta )
##pp_check(fit_sigmoid_c_beta )
#
## Next we can infer the confidence intervals of the resistant proportion 
## and the difference between patients
#
#fit_sigmoid_c_patient_beta <- brm(
#  bf(value ~ (1-c) * (1 / (1 + exp(-k * (time - x0)))),c ~ patient ,k ~1, x0 ~1,nl=TRUE),
#  data = prop_data,
#  prior = c(prior(gamma(0.01, 0.01), class = "phi"),
#            prior(normal(0.1, 0.05), class = "b", nlpar = "c",ub=1,lb=0),   # Prior for the asymptote (lmax)
#            prior(normal(0, 0.5), class = "b", nlpar = "k"),   # Prior for the slope (k)
#            prior(normal(162, 10), class = "b", nlpar = "x0") # Prior for the midpoint (x0)
#  ), 
#  family = brmsfamily(family="Beta", link="identity", link_phi="identity"),
#  #init = 0, 
#  chains = 4,
#  cores = 2,
#  threads = 4,
#  iter = 8000,
#  warmup = 4000,
#  control = list(adapt_delta = 0.95,max_treedepth = 12)
#)
#
##summary(fit_sigmoid_c_patient_beta)
##pp_check(fit_sigmoid_c_patient_beta)
#
#fit_hill_c_beta <- brm(
#  bf(value ~ (1-c) * (1*time^alpha / (time^alpha + Km^alpha)),c ~ 1 ,Km ~1, alpha ~1,nl=TRUE),
#  data = prop_data,
#  prior = c(prior(gamma(0.01, 0.01), class = "phi"),
#            prior(normal(0.1, 0.05), class = "b", nlpar = "c",ub=1,lb=0),   # Prior for the asymptote (lmax)
#            prior(normal(200, 100), class = "b", nlpar = "Km"),   # Prior for the slope (k)
#            prior(normal(2, 1), class = "b", nlpar = "alpha") # Prior for the midpoint (x0)
#  ), 
#  family = brmsfamily(family="Beta", link="identity", link_phi="identity"),
#  #init = 0, 
#  chains = 4,
#  cores = 4,
#  threads = 4,
#  iter = 8000,
#  warmup = 4000,
#  control = list(adapt_delta = 0.95,max_treedepth = 12)
#)
#
##summary(fit_hill_c_beta)
##pp_check(fit_hill_c_beta)
#
#fit_hill_c_patient_beta <- brm(
#  bf(value ~ (1-c) * (1*time^alpha / (time^alpha + Km^alpha)),c ~ patient ,Km ~1, alpha ~1,nl=TRUE),
#  data = prop_data,
#  prior = c(prior(gamma(0.01, 0.01), class = "phi"),
#            prior(normal(0.1, 0.05), class = "b", nlpar = "c",ub=1,lb=0),   # Prior for the asymptote (lmax)
#            prior(normal(200, 100), class = "b", nlpar = "Km"),   # Prior for the slope (k)
#            prior(normal(2, 1), class = "b", nlpar = "alpha") # Prior for the midpoint (x0)
#  ), 
#  family = brmsfamily(family="Beta", link="identity", link_phi="identity"),
#  #init = 0, 
#  chains = 4,
#  cores = 4,
#  threads = 4,
#  iter = 8000,
#  warmup = 4000,
#  control = list(adapt_delta = 0.95,max_treedepth = 12)
#)
#
##summary(fit_hill_c_patient_beta)
##pp_check(fit_hill_c_patient_beta)
#
#fit_gompertz_c_beta <- brm(
#  bf(value ~ (1-c) * 1 * exp(-b * exp(-d * time)),c ~ 1 ,d ~1, b ~1,nl=TRUE),
#  data = prop_data,
#  prior = c(prior(gamma(0.01, 0.01), class = "phi"),
#            prior(normal(0.1, 0.05), class = "b", nlpar = "c",ub=1,lb=0),   # Prior for the asymptote
#            prior(normal(3, 1), class = "b", nlpar = "b"),   # Prior for the displacement along x-axis
#            prior(normal(0.01, 0.1), class = "b", nlpar = "d") # Prior for the growth rate
#  ), 
#  family = brmsfamily(family="Beta", link="identity", link_phi="identity"),
#  chains = 4,
#  cores = 4,
#  threads = 2,
#  iter = 4000,
#  warmup = 2000,
#  control = list(adapt_delta = 0.95,max_treedepth = 12)
#)
#
##summary(fit_gompertz_c_beta)
##pp_check(fit_gompertz_c_beta)
##predictions <- simulate_for_pred(t0,tn,sim_n,patient_var,fit_gompertz_c_beta)
#
##summary(fit_gompertz_c_patient)
##pp_check(fit_gompertz_c_patient)
#
#fit_gompertz_c_patient_beta <- brm(
#  bf(value ~ (1-c) * 1 * exp(-b * exp(-d * time)),c ~ patient ,d ~1, b ~1,nl=TRUE),
#  data = prop_data,
#  prior = c(prior(gamma(0.01, 0.01), class = "phi"),
#            #prior(normal(0, 0.1), class = "b", coef = "patientPDO2",nlpar = "c"),
#            prior(normal(0.1, 0.05), class = "b", nlpar = "c",ub=1,lb=0),   # Prior for the asymptote
#            prior(normal(3, 1), class = "b", nlpar = "b"),   # Prior for the displacement along x-axis
#            prior(normal(0.01, 0.1), class = "b", nlpar = "d") # Prior for the growth rate
#  ), 
#  family = brmsfamily(family="Beta", link="identity", link_phi="identity"),
#  #init = 0, 
#  chains = 4,
#  cores = 4,
#  threads = 4,
#  iter = 8000,
#  warmup = 4000,
#  control = list(adapt_delta = 0.95,max_treedepth = 12)
#)
#
##summary(fit_gompertz_c_patient_beta)
##pp_check(fit_gompertz_c_patient_beta)

#saveRDS(fit_hill_c_beta, "brms_fit_hill_c_beta.rds")
#saveRDS(fit_sigmoid_c_beta, "brms_fit_sigmoid_c_beta.rds")
#saveRDS(fit_gompertz_c_beta, "brms_fit_gompertz_c_beta.rds")

#saveRDS(fit_hill_c_patient_beta, "brms_fit_hill_c_patient_beta.rds")
#saveRDS(fit_sigmoid_c_patient_beta, "brms_fit_sigmoid_c_patient_beta.rds")
#saveRDS(fit_gompertz_c_patient_beta, "brms_fit_gompertz_c_patient_beta.rds")

####################################
### Compare model fits using loo ###
####################################

fit_hill_c_beta <- readRDS("../data/models/brms_fit_hill_c_beta.rds")
fit_sigmoid_c_beta <- readRDS("../data/models/brms_fit_sigmoid_c_beta.rds")
fit_gompertz_c_beta <- readRDS("../data/models/brms_fit_gompertz_c_beta.rds")

fit_hill_c_beta <- add_criterion(fit_hill_c_beta, "loo")
fit_sigmoid_c_beta <- add_criterion(fit_sigmoid_c_beta, "loo")
fit_gompertz_c_beta <- add_criterion(fit_gompertz_c_beta, "loo")
#loo_compare(fit_hill_c_beta,fit_sigmoid_c_beta,fit_gompertz_c_beta)

fit_hill_c_patient_beta <- readRDS("../data/models/brms_fit_hill_c_patient_beta.rds")
fit_sigmoid_c_patient_beta <- readRDS("../data/models/brms_fit_sigmoid_c_patient_beta.rds")
fit_gompertz_c_patient_beta <- readRDS("../data/models/brms_fit_gompertz_c_patient_beta.rds")

fit_hill_c_patient_beta <- add_criterion(fit_hill_c_patient_beta, "loo")
fit_sigmoid_c_patient_beta <- add_criterion(fit_sigmoid_c_patient_beta, "loo")
fit_gompertz_c_patient_beta <- add_criterion(fit_gompertz_c_patient_beta, "loo")
#loo_compare(fit_hill_c_patient_beta,fit_sigmoid_c_patient_beta,fit_gompertz_c_patient_beta)

#loo_compare(fit_hill_c_beta,
#            fit_sigmoid_c_beta,
#            fit_gompertz_c_beta,
#            fit_hill_c_patient_beta,
#            fit_sigmoid_c_patient_beta,
#            fit_gompertz_c_patient_beta
#)

# Note the model with lowest elpd_diff is best

#################################################################
### Visualize the posterior of the basic logistic with offset ###
#################################################################

# We can sample from the posterior to get values for an error ribbon
#summary(fit_sigmoid_c_beta)
lmax <- 0.90
k <- 0.02 
x0 <- 165.45
results <- data.frame()
time_values <- seq(0, 600, length.out = 100)
# Calculate logistic values for each k
values <- logistic_function(time_values, lmax, k, x0)
sigmoid_results <- rbind(results, data.frame(time = time_values, value = values))

t0 = 0
tn = 600
sim_n = 100
patient_var = "PDO1"
sigmoid_predictions <- simulate_for_pred(t0,tn,sim_n,patient_var,fit_sigmoid_c_beta)

##########################################################################
### Visualize the posterior of the patient effect logistic with offset ###
##########################################################################

#summary(fit_sigmoid_c_patient_beta)
# We can sample from the posterior to get values for an error ribbon
lmax <- 0.91
lmax_p2 <- 0.90 
k <- 0.02 
x0 <- 165.47
results <- data.frame()
time_values <- seq(0, 600, length.out = 100)
# Calculate logistic values for each k
values_p1 <- logistic_function(time_values, lmax, k, x0)
sigmoid_p_results_p1 <- rbind(results, data.frame(time = time_values, value = values_p1))
values_p2 <- logistic_function(time_values, lmax_p2, k, x0)
sigmoid_p_results_p2 <- rbind(results, data.frame(time = time_values, value = values_p2))

patient_var = "PDO1"
sigmoid_p_predictions_p1 <- simulate_for_pred(t0,tn,sim_n,patient_var,fit_sigmoid_c_patient_beta)
patient_var = "PDO2"
sigmoid_p_predictions_p2 <- simulate_for_pred(t0,tn,sim_n,patient_var,fit_sigmoid_c_patient_beta)

###############################################################################
### Visualize the posterior of the patient effect Hill function with offset ###
###############################################################################

#summary(fit_hill_c_beta)
# We can sample from the posterior to get values for an error ribbon
Km = 158.64
Vmax <- 1
alpha <- 1.99
c <- 0.0

hill_results <- data.frame()

# Calculate logistic values for each k
time_values <- seq(0, 600, length.out = 100)
values <- hill_function(time_values, c, Vmax, Km, alpha)
hill_results<- rbind(results, data.frame(time = time_values, value = values))

patient_var <- "PDO1"
hill_predictions <- simulate_for_pred(t0,tn,sim_n,patient_var,fit_hill_c_beta)

##################################################################################
### Visualize the posterior of the patient effect Hill function without offset ###
##################################################################################

#summary(fit_hill_c_patient_beta)
# We can sample from the posterior to get values for an error ribbon
Km = 157.80
Vmax <- 1
Vmax_p2 <- 0.99
alpha <- 2.01
c <- 0

results <- data.frame()

# Calculate logistic values for each k
time_values <- seq(0, 600, length.out = 100)
values_p1 <- hill_function(time_values, c, Vmax, Km, alpha)
hill_p_results_p1 <- rbind(results, data.frame(time = time_values, value = values_p1))
values_p2 <- hill_function(time_values, c, Vmax_p2, Km, alpha)
hill_p_results_p2 <- rbind(results, data.frame(time = time_values, value = values_p2))

patient_var = "PDO1"
hill_p_predictions_p1 <- simulate_for_pred(t0,tn,sim_n,patient_var,fit_hill_c_patient_beta)
patient_var = "PDO2"
hill_p_predictions_p2 <- simulate_for_pred(t0,tn,sim_n,patient_var,fit_hill_c_patient_beta)

#####################################################################
### Visualize the posterior of the  Gompertz function with offset ###
#####################################################################

#summary(fit_gompertz_c_beta)
# We can sample from the posterior to get values for an error ribbon
c <- 0.07
d <- 0.01
b <- 3.43

results <- data.frame()

# Calculate logistic values for each k
time_values <- seq(0, 600, length.out = 100)
values <- gompertz_function(time_values, 1-c, b, d)
gompertz_results <- rbind(results, data.frame(time = time_values, value = values))

patient_var <- "PDO1"
gompertz_predictions <- simulate_for_pred(t0,tn,sim_n,patient_var,fit_gompertz_c_beta)

###################################################################################
### Visualize the posterior of the patient effect Gompertz function with offset ###
###################################################################################

#summary(fit_gompertz_c_patient_beta)
# We can sample from the posterior to get values for an error ribbon
c <- 0.06
c_p2 <- 0.07
d <- 0.01
b <- 3.44

results <- data.frame()

# Calculate logistic values for each k
time_values <- seq(0, 600, length.out = 100)
values_p1 <- gompertz_function(time_values, 1-c, b, d)
gompertz_p_results_p1 <- rbind(results, data.frame(time = time_values, value = values_p1))
values_p2 <- gompertz_function(time_values, 1-c_p2, b, d)
gompertz_p_results_p2 <- rbind(results, data.frame(time = time_values, value = values_p2))

patient_var = "PDO1"
gompertz_p_predictions_p1 <- simulate_for_pred(t0,tn,sim_n,patient_var,fit_gompertz_c_patient_beta)
patient_var = "PDO2"
gompertz_p_predictions_p2 <- simulate_for_pred(t0,tn,sim_n,patient_var,fit_gompertz_c_patient_beta)

################################################
### Plot posterior of resistant proportion c ###
################################################

posterior_all <- fit_gompertz_c_patient_beta %>% 
  gather_draws(`b_.*`, regex = TRUE)

posterior_all_cast <- posterior_all %>% pivot_wider(names_from = .variable,
                                                    values_from = .value) %>%
                      mutate(b_c_Intercept_withPDO2 = b_c_Intercept +b_c_patientPDO2) %>%
                     select(-c(b_c_patientPDO2,b_d_Intercept,b_b_Intercept)) %>%
                     pivot_longer(cols = starts_with("b"), names_to = ".variable", values_to = ".value")

posterior_all_cast$.variable <- gsub("b_c_Intercept_withPDO2","PDO2",posterior_all_cast$.variable)
posterior_all_cast$.variable <- gsub("b_c_Intercept","PDO1",posterior_all_cast$.variable)

# Now add the results from the model without the patient coefficient
posterior_all_nopatient <- fit_gompertz_c_beta %>% 
  gather_draws(`b_.*`, regex = TRUE)
posterior_c_nopatient <- posterior_all_nopatient[posterior_all_nopatient$.variable=="b_c_Intercept",]
posterior_c_nopatient$.variable <- gsub("b_c_Intercept","PDO1+PDO2",posterior_c_nopatient$.variable)

posterior_c_comp <- rbind(posterior_c_nopatient,posterior_all_cast)

###############################
###  Make all final figures ###
###############################

model_comparison_loo <- loo_compare(fit_hill_c_beta,
                                    fit_sigmoid_c_beta,
                                    fit_gompertz_c_beta,
                                    fit_hill_c_patient_beta,
                                    fit_sigmoid_c_patient_beta,
                                    fit_gompertz_c_patient_beta
)

#report(model_comparison_loo)

# The rule of thumb is that the models are "very similar" if |elpd_diff| (the absolute value of elpd_diff) is less than 4 (Sivula, Magnusson and Vehtari, 2020)
model_comparison_loo_df <- as.data.frame(model_comparison_loo)
model_comparison_loo_df  <- tibble::rownames_to_column(model_comparison_loo_df , "model")

# Fig 1 Logistic model
sigmoid_elpd <- round(model_comparison_loo_df[model_comparison_loo_df$model=="fit_sigmoid_c_beta",]$elpd_loo,2)
f1 <- ggplot(prop_data, aes(x = time, y = value)) +
  geom_point(alpha=0.6) +
  geom_line(data = sigmoid_results, aes(x = time, y = value), color = "black",linewidth=0.7) +
  geom_ribbon(data = sigmoid_predictions, aes(x = time, ymin = lower_ci, ymax = upper_ci), 
              fill = "lightgrey", alpha = 0.5,color="lightgrey",linewidth=0.5,linetype = 2) +
  xlab("Time (minutes)")+
  ylab("Apoptotic organoids")+
  ggtitle(paste("Logistic ELPD=",sigmoid_elpd ,sep=""))+
  scale_y_continuous(labels = label_percent()) +
 #scale_color_manual(values = c("PDO1"="orange","PDO2"="blue"))+
  theme_minimal()+
  theme(text = element_text(size=12),plot.title = element_text(size=12),legend.position="none")

# Fig 2 Logist model with patient coef

sigmoid_p_elpd <- round(model_comparison_loo_df[model_comparison_loo_df$model=="fit_sigmoid_c_patient_beta",]$elpd_loo,2)
f2 <- ggplot(prop_data, aes(x = time, y = value,color = patient)) +
  geom_point(alpha=0.6) +
  geom_line(data = sigmoid_p_results_p1, aes(x = time, y = value), color = "orange", linewidth=0.7) +
  geom_line(data = sigmoid_p_results_p2, aes(x = time, y = value), color = "blue",linewidth=0.7) +
  geom_ribbon(data = sigmoid_p_predictions_p1, aes(x = time, ymin = lower_ci, ymax = upper_ci), 
              fill = "lightgrey", alpha = 0.2,color="orange",linewidth=0.5,linetype = 2) +
  geom_ribbon(data = sigmoid_p_predictions_p2, aes(x = time, ymin = lower_ci, ymax = upper_ci), 
              fill = "lightgrey", alpha = 0.2,color="blue",linewidth=0.5,linetype = 2) +
  guides(fill = "none") +
  xlab("Time (minutes)")+
  ylab("Apoptotic organoids")+
  scale_y_continuous(labels = label_percent()) +
  ggtitle(paste("Logistic ELPD=",sigmoid_elpd ,sep=""))+
  scale_color_manual(values = c("PDO1"="orange","PDO2"="blue"))+
  theme_minimal()+
  theme(text = element_text(size=12),plot.title = element_text(size=12),legend.position="none")

hill_elpd <- round(model_comparison_loo_df[model_comparison_loo_df$model=="fit_hill_c_beta",]$elpd_loo,2)
f3 <- ggplot(prop_data, aes(x = time, y = value,)) +
  geom_point(alpha=0.6) +
  geom_line(data = hill_results, aes(x = time, y = value), color = "black",linewidth=0.7) +
  geom_ribbon(data = hill_predictions, aes(x = time, ymin = lower_ci, ymax = upper_ci), 
              fill = "lightgrey", alpha = 0.5,color="lightgrey",linewidth=0.5,linetype = 2) +
  guides(fill = "none") +
  xlab("Time (minutes)")+
  ylab("Apoptotic organoids")+
  scale_y_continuous(labels = label_percent()) +
  ggtitle(paste("Hill ELPD=",hill_elpd,sep=""))+
  #scale_color_manual(values = c("PDO1"="orange","PDO2"="blue"))+
  theme_minimal()+
  theme(text = element_text(size=12),plot.title = element_text(size=12),legend.position="none")

hill_p_elpd <- round(model_comparison_loo_df[model_comparison_loo_df$model=="fit_hill_c_patient_beta",]$elpd_loo,2)
f4 <- ggplot(prop_data, aes(x = time, y = value,color = patient)) +
  geom_point(alpha=0.6) +
  geom_line(data = hill_p_results_p1, aes(x = time, y = value), color = "orange",linewidth=0.7) +
  geom_line(data = hill_p_results_p2, aes(x = time, y = value), color = "blue",linewidth=0.7) +
  geom_ribbon(data = hill_p_predictions_p1, aes(x = time, ymin = lower_ci, ymax = upper_ci), 
              fill = "lightgrey", alpha = 0.2,color="orange",linewidth=0.5,linetype = 2) +
  geom_ribbon(data = hill_p_predictions_p2, aes(x = time, ymin = lower_ci, ymax = upper_ci), 
              fill = "lightgrey", alpha = 0.2,color="blue",linewidth=0.5,linetype = 2) +
  xlab("Time (minutes)")+
  ylab("Apoptotic organoids")+
  scale_y_continuous(labels = label_percent()) +
  ggtitle(paste("Hill ELPD=",hill_p_elpd,sep=""))+
  scale_color_manual(values = c("PDO1"="orange","PDO2"="blue"))+
  theme_minimal()+
  theme(text = element_text(size=12),plot.title = element_text(size=12),legend.position="none")

gompertz_elpd <- round(model_comparison_loo_df[model_comparison_loo_df$model=="fit_gompertz_c_beta",]$elpd_loo,2)
f5 <- ggplot(prop_data, aes(x = time, y = value)) +
  geom_point(alpha=0.6) +
  geom_line(data = gompertz_results, aes(x = time, y = value), color = "black",linewidth=0.7) +
  geom_ribbon(data = gompertz_predictions, aes(x = time, ymin = lower_ci, ymax = upper_ci), 
              fill = "lightgrey", alpha = 0.5,color="lightgrey",linewidth=0.5,linetype = 2) +
  xlab("Time (minutes)")+
  ylab("Apoptotic organoids")+
  scale_y_continuous(labels = label_percent()) +
  ggtitle(paste("Gompertz ELPD=",gompertz_elpd,sep=""))+
  #scale_color_manual(values = c("PDO1"="orange","PDO2"="blue"))+
  theme_minimal()+
  theme(text = element_text(size=12),plot.title = element_text(size=12),legend.position="none")

gompertz_p_elpd <- round(model_comparison_loo_df[model_comparison_loo_df$model=="fit_gompertz_c_patient_beta",]$elpd_loo,2)
f6 <- ggplot(prop_data, aes(x = time, y = value,color = patient)) +
  geom_point(alpha=0.6) +
  geom_line(data = gompertz_p_results_p1, aes(x = time, y = value), color = "orange",linewidth=0.7) +
  geom_line(data = gompertz_p_results_p2, aes(x = time, y = value), color = "blue",linewidth=0.7) +
  geom_ribbon(data = gompertz_p_predictions_p1, aes(x = time, ymin = lower_ci, ymax = upper_ci), 
              fill = "lightgrey", alpha = 0.2,color="orange",linewidth=0.5,linetype = 2) +
  geom_ribbon(data = gompertz_p_predictions_p2, aes(x = time, ymin = lower_ci, ymax = upper_ci), 
              fill = "lightgrey", alpha = 0.2,color="blue",linewidth=0.5,linetype = 2) +
  xlab("Time (minutes)")+
  ylab("Apoptotic organoids")+
  scale_y_continuous(labels = label_percent()) +
  ggtitle(paste("Gompertz ELPD=",gompertz_p_elpd,sep=""))+
  scale_color_manual(values = c("PDO1"="orange","PDO2"="blue"))+
  theme_minimal()+
  theme(text = element_text(size=12),plot.title = element_text(size=12),legend.position="none")

f0 <- ggplot(prop_data, aes(x = time, y = value,color = patient)) +
  geom_point() +
  scale_color_manual(values = c("PDO1"="orange","PDO2"="blue"))+
  guides(color = guide_legend(nrow = 1)) +
  theme_minimal()+
  theme(text = element_text(size=12),plot.title = element_text(size=12),legend.position = "bottom",legend.title=element_blank())
legend_f0 <- cowplot::get_plot_component(f0, 'guide-box-bottom', return_all = TRUE)

p1 <- plot_grid(f1,f2,f3,f4,f5,f6, ncol=2,labels = "AUTO")
plot_grid(p1,legend_f0,ncol=1,rel_heights = c(0.95, 0.05))
ggsave("../figures/Sigmoid_model_comparison.pdf", height = 12, width = 8)

f7 <- ggplot(prop_data, aes(x = time, y = value,color = patient)) +
  geom_point(alpha=0.6,size=2.5) +
  geom_line(data = gompertz_p_results_p1, aes(x = time, y = value), color = "orange",linewidth=1) +
  geom_line(data = gompertz_p_results_p2, aes(x = time, y = value), color = "blue",linewidth=1) +
  geom_ribbon(data = gompertz_p_predictions_p1, aes(x = time, ymin = lower_ci, ymax = upper_ci), 
              fill = "lightgrey", alpha = 0.2,color="orange",linewidth=0.7,linetype = 2) +
  geom_ribbon(data = gompertz_p_predictions_p2, aes(x = time, ymin = lower_ci, ymax = upper_ci), 
              fill = "lightgrey", alpha = 0.2,color="blue",linewidth=0.7,linetype = 2) +
  xlab("Time (minutes)")+
  ylab("Apoptotic organoids")+
  scale_y_continuous(labels = label_percent()) +
  scale_color_manual(values = c("PDO1"="orange","PDO2"="blue"))+
  theme_minimal()+
  theme(text = element_text(size=22),legend.position="bottom",legend.title=element_blank())
ggsave("../figures/Gompertz_best_fit_model.pdf")

f8 <- ggplot(posterior_c_comp, aes(x =.value, y=.variable,fill=.variable)) +
  stat_halfeye(point_interval = "median_hdi",slab_alpha = 0.75) +
  ylab("Density")+
  xlab("Resistant organoids")+
  guides(fill = "none") +
  scale_x_continuous(labels = label_percent()) +
  #scale_fill_manual(values = c("PDO1"="#fc8d62","PDO2"="#66c2a5","PD01+PDO2"="#e7298a")) +
  scale_fill_manual(values = c("#fc8d62","#e7298a","#66c2a5")) +
  theme_bw()+
  theme(text = element_text(size = 18))
ggsave("../figures/Gompertz_best_fit_model_resistant_organoids_by_patient.pdf")

