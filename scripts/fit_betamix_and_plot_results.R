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
library(cowplot)
library(grid)

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

###############################
## Load and pre-process data ##
###############################

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
pd_non01 <- pd %>%
  filter(!is.na(value))  %>%
  filter(value<=99) %>%
  group_by(treatment,time,expt,patient) %>%
  mutate(value = (((value/100)*(n()-1))+s)/n()) %>%
  mutate(cluster = ifelse(value<=0.05,1,2)) %>%
  mutate(k1 = ifelse(value<=0.01,1,ifelse(value<=0.05,0.95,0)),
         k2 = 1 - k1 ) %>%
  filter(!is.na(value))
full_pd <- pd_non01[!is.na(pd_non01$value) & pd_non01$time == 585 & pd_non01$treatment == "NK",]

############################
## Bayesian model fitting ##
############################

mix <- mixture(Beta,Beta)

# The Intercept priors mu1 and mu2 are logit transformed while the other parameters are not
# plogis() will convert it back, qloqis() is the inverse function that is applied for the initial transformation
priors_notheta_notail_norm <- c(
  prior(normal(-2,0.1), class = "Intercept", dpar = mu1),
  prior(normal(0, 2), class = "Intercept", dpar = mu2),
  prior(gamma(0.01, 0.01), class="phi1",lb=1.5,ub=3),
  prior(gamma(0.15, 0.01), class="phi2",lb=5)
)
# IMPORTANT: When estimating mixing proportions in brms, 
# estimates are indeed on the logit scale and priors should be specified as such.
# the value theta2_Intercept can thus be transformed with plogis() to get the proportion

#fit_fixall_notailprior_norm <- brm(bf(value ~ 1,theta2 ~ patient),
#                              data= full_pd, family = mix,init = 0, 
#                              control = list(adapt_delta = 0.97,max_treedepth = 12),
#                              prior = priors_notheta_notail_norm, chains = 4, cores = 4,threads = 2,iter = 4000)

#summary(fit_fixall_notailprior_norm)
#saveRDS(fit_fixall_notailprior_norm, "../data/models/brms_fit_fixall_notailprior_norm_21102024.rds")

######################################
### Process and plot fitted models ###
######################################

# The model will be loaded from a pre-generated file to ensure consistency
fit_fixall_notailprior_norm <- readRDS("../data/models/brms_fit_fixall_notailprior_norm_21102024.rds")

pred_beta_theta2 <- fit_fixall_notailprior_norm %>% 
  epred_draws(newdata = expand_grid(patient = c("p1", "p2")),
              dpar = "theta2",re_formula=NA)
mfx_patient_theta2 <- pred_beta_theta2 %>% 
  compare_levels(variable = theta2, by = patient)
label_pp <- label_number(accuracy = 1, scale = 100, 
                         suffix = " pp.", style_negative = "minus")

# Annotate resistant and susceptible
fit_mix <- as.data.frame(pp_mixture(fit_fixall_notailprior_norm))
full_pd$k1_pp <- fit_mix$`Estimate.P(K = 1 | Y)`
full_pd <- full_pd %>%
  mutate(component = ifelse(k1_pp>=0.5,"resistant","susceptible"))

### Hard coded model params ###

# These params can also be extracted from the model

mu1 = 0.1202569
phi1 = 2.35
mu2 = 0.4378235
phi2 = 6.85

# Convert to alternative beta parametrization
a1 <- mu1 * phi1
b1 <- (1 - mu1) * phi1
a2 <- mu2 * phi2
b2 <- (1 - mu2) * phi2

### Show the changes in mu1 and mu2 over time ###

# Method 1
# Use endpoint fit to infer thresholds for resistant and susceptible groups

resistant_threshold <- max(full_pd[full_pd$component=="resistant",]$value) * 100

pd_class <- pd[pd$treatment=="NK",] %>% mutate(component=ifelse(value<=resistant_threshold,"resistant","susceptible"))
pd_all <- pd %>% mutate(component=ifelse(value<=resistant_threshold,"resistant","susceptible"))%>%
  mutate(component=ifelse(treatment!="NK","Control",component))

pd_class$time <- as.numeric(pd_class$time)
pd_all$time <- as.numeric(pd_all$time)

pd_class_min <- pd_class %>% filter(time %in% c(15,165,270,375,480,585))

pd_class_stat <- pd_class[complete.cases(pd_class),] %>%
  group_by(time,expt,component,patient) %>%
  summarize(mean_apoptosis = mean(value,na.rm = TRUE),
            sd_apoptosis = sd(value,na.rm = TRUE),
            N_N=n(),  
            se=sd_apoptosis/sqrt(N_N),  
            upper_limit=mean_apoptosis+se,  
            lower_limit=mean_apoptosis-se )

pd_all_stat <- pd_all[complete.cases(pd_all),] %>%
  group_by(treatment,time,expt,component,patient) %>%
  summarize(mean_apoptosis = mean(value,na.rm = TRUE),
            sd_apoptosis = sd(value,na.rm = TRUE),
            N_N=n(),  
            se=sd_apoptosis/sqrt(N_N),  
            upper_limit=mean_apoptosis+se,  
            lower_limit=mean_apoptosis-se )

pd_all_stat_pool <- pd_all[complete.cases(pd_all),] %>%
  group_by(treatment,time,component,patient) %>%
  summarize(mean_apoptosis = mean(value,na.rm = TRUE),
            sd_apoptosis = sd(value,na.rm = TRUE),
            N_N=n(),  
            se=sd_apoptosis/sqrt(N_N),  
            upper_limit=mean_apoptosis+(se*2),  
            lower_limit=mean_apoptosis-(se*2))

# Compare proportions of apoptotic and non-apoptotic organoids over time
pd_theta_stat <- pd_all[complete.cases(pd_all),] %>%
  mutate(component=ifelse(value<=resistant_threshold,"resistant","susceptible")) %>%
  group_by(treatment,time,expt,patient) %>%
  summarize(n_apoptotic = sum(component == "susceptible"),
            n_nonapoptotic = sum(component == "resistant")) %>%
  mutate(n_apoptotic_perc = (n_apoptotic/(n_apoptotic+n_nonapoptotic))*100,
         n_nonapoptotic_perc = (n_nonapoptotic/(n_apoptotic+n_nonapoptotic))*100)

# Summarize reps
pd_theta_stat_err <- pd_theta_stat  %>%
  group_by(treatment,time,patient) %>%
  summarize(mean_apoptotic_perc = mean(n_apoptotic_perc),
            mean_nonapoptotic_perc = mean(n_nonapoptotic_perc),
            max_apoptotic_perc = max(n_apoptotic_perc),
            min_apoptotic_perc = min(n_apoptotic_perc),
            max_nonapoptotic_perc = max(n_nonapoptotic_perc),
            min_nonapoptotic_perc = min(n_nonapoptotic_perc))

pd_class_stat_min <- pd_class_stat %>% filter(time %in% c(15,165,270,375,480,585))

# Plots for report 

# Plot 1 Apoptosis time series split by control/resistant/susceptible
pd_all_stat$patient <- gsub("p","PDO",pd_all_stat$patient)
pd_all_stat$component <- gsub("Control","control",pd_all_stat$component)
p1_cols <- c("control" = "grey", "resistant" = "blue", "susceptible" = "orange")

ggplot(pd_all_stat, aes(x=time, y=mean_apoptosis,ymax=upper_limit,ymin=lower_limit,fill=component,shape=expt)) +  
  geom_line() +  
  geom_ribbon(alpha=0.5)+
  facet_wrap(patient ~.)+
  ylab("Mean organoid apoptotic surface area (%)")+
  xlab("Time (minutes)")+
  scale_fill_manual(values = p1_cols)+
  theme_pubclean()+
  theme(text = element_text(size=18),legend.title=element_blank())

ggsave("../figures/Mean_Trajectory_Apoptotic_SA.pdf")

# Plot 2
full_pd$component_expt <- paste(full_pd$component,full_pd$expt)
p2_cols <- c("susceptible exp15"= "orange",
             "resistant exp15"= "blue",
             "susceptible exp17"= "orange",
             "resistant exp17"= "blue",
             "susceptible exp18"= "orange",
             "resistant exp18"= "blue",  
             "susceptible exp13"= "orange",
             "resistant exp13"= "blue",
             "susceptible exp14"= "orange",
             "resistant exp14"= "blue",
             "susceptible exp16"= "orange",
             "resistant exp16"= "blue")

# Add the fitted curve

bw=0.05
theta2 <- 0.9315024
theta1 <- 1-theta2
theta2_p2_effect <- plogis(2.61 + -0.82) - plogis(2.61)

# p1 and p2 are different wrt theta
hist1 <- ggplot(full_pd[full_pd$time== 585 & full_pd$patient=="p1",],aes(x=value,fill=component_expt))+ 
  geom_histogram(alpha=0.4,position="identity",binwidth=bw,color="grey")+
  ylab("Count")+
  xlab("Proportion organoid apoptotic area")+
  ylim(0,125)+
  xlim(NA,1)+
  scale_fill_manual(values = p2_cols)+
  stat_function(fun = function(x) 
    dbeta(x, shape1 = a1, shape2=b1)*  bw * length(full_pd[full_pd$time== 585 & full_pd$expt=="exp17",]$value)*theta1,color="black")+
  stat_function(fun = function(x) 
    dbeta(x, shape1 = a2, shape2=b2)* bw * length(full_pd[full_pd$time== 585 & full_pd$expt=="exp17",]$value)*theta2,color="black")+
  theme_pubclean()+
  theme(legend.position="none",text=element_text(size=18))

hist2 <- ggplot(full_pd[full_pd$time== 585 & full_pd$patient=="p2",],aes(x=value,fill=component_expt))+ 
  geom_histogram(alpha=0.4,position="identity",binwidth=bw,color="grey")+
  ylab("Count")+
  xlab("Proportion organoid apoptotic area")+
  ylim(0,125)+
  xlim(NA,1)+
  scale_fill_manual(values = p2_cols)+
  stat_function(fun = function(x) 
    dbeta(x, shape1 = a1, shape2=b1)*  bw * length(full_pd[full_pd$time== 585 & full_pd$expt=="exp16",]$value)*(theta1-theta2_p2_effect),color="black")+
  stat_function(fun = function(x) 
    dbeta(x, shape1 = a2, shape2=b2)* bw * length(full_pd[full_pd$time== 585 & full_pd$expt=="exp16",]$value)*(theta2+theta2_p2_effect),color="black")+
  theme_pubclean()+
  theme(legend.position="none",text=element_text(size=18))

plot_grid(hist1,hist2)
ggsave("../figures/Histograms_Apoptotic_SA.pdf",width=10)

#  Plot 3 Difference in theta between p1 and p2 (Marginal effect of patient2)
ggplot(mfx_patient_theta2, aes(x = -theta2)) +
  #scale_x_continuous(labels = label_pp) +
  stat_halfeye(.width = c(0.8, 0.95), point_interval = "median_hdi",
               fill = "#fb9e07")+
  scale_x_continuous(labels = label_percent()) +
  ylab("Density")+
  xlab("Increase in resistant organoids in PDO2")+
  theme_pubclean()+
  theme(text = element_text(size = 18))
ggsave("../figures/Marginal_effect_of_patient2_coeff.pdf")

# Plot 4 Compare percent apoptotic and nonapoptotic

pd_theta_stat$patient <- gsub("p","PDO",pd_theta_stat$patient)
pd_theta_stat <- melt(pd_theta_stat,id.vars = c("treatment",
                                                "time",
                                                "expt",
                                                "patient",
                                                "n_apoptotic",        
                                                "n_nonapoptotic"))

p1_cols <- c("control" = "grey", "NK" = "orange")
pd_theta_stat$patient <- gsub("p","PDO",pd_theta_stat$patient)
pd_theta_stat_err$patient <- gsub("p","PDO",pd_theta_stat_err$patient)

pd_theta_stat_err$treatment <- gsub("NK","co-culture",pd_theta_stat_err$treatment)
ggplot(pd_theta_stat_err, aes(x=time, y=mean_apoptotic_perc,ymax=max_apoptotic_perc,ymin=min_apoptotic_perc,fill=treatment)) +  
  geom_line() +  
  geom_ribbon(alpha=0.5)+
  facet_wrap(patient ~ .,scales="free",axis.labels = "all")+
  ylab("Apoptotic organoids (%)")+
  xlab("Time (minutes)")+
  scale_fill_manual(values = c("co-culture"="orange","control"="grey"))+
  theme_pubclean()+
  theme(text = element_text(size=18),legend.title=element_blank(),panel.spacing = unit(2, "lines"))
ggsave("../figures/Apoptotic_organoids_percent_over_time.pdf",width=8)

ggplot(pd_theta_stat_err, aes(x=time, y=mean_nonapoptotic_perc,ymax=max_nonapoptotic_perc,ymin=min_nonapoptotic_perc,fill=treatment)) +  
  geom_line() +  
  geom_ribbon(alpha=0.5)+
  facet_wrap(patient ~ .,scales="free",axis.labels = "all")+
  ylab("Non-apoptotic organoids (%)")+
  xlab("Time (minutes)")+
  scale_fill_manual(values = c("co-culture"="orange","control"="grey"))+
  theme_pubclean()+
  theme(text = element_text(size=18),legend.title=element_blank(),panel.spacing = unit(2, "lines"))
ggsave("../figures/Nonapoptotic_organoids_percent_over_time.pdf",width=8)

# Reformat table to also use apoptotic and nonapoptotic as group
pd_theta_stat_err_melt <- melt(pd_theta_stat_err, id.vars = c("treatment",
                                                              "time",
                                                              "patient"
                                                              ))
pd_theta_stat_err_melt <- separate_wider_delim(pd_theta_stat_err_melt, cols = variable, delim = "_", names = c("state", "variable","other")) %>%
  select(-other)
pd_theta_stat_err_dcast <- pivot_wider(pd_theta_stat_err_melt,names_from = state,values_from = value)
pd_theta_stat_err_dcast$treatment_variable <- paste(pd_theta_stat_err_dcast$treatment,pd_theta_stat_err_dcast$variable)
pd_theta_stat_err_dcast$patient <- gsub("p","PDO",pd_theta_stat_err_dcast$patient)

