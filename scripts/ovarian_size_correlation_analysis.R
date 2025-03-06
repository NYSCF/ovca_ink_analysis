library(dplyr)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(cowplot)
library(tidyr)
library(ggExtra)
library(ggpointdensity)
library(viridis)
library(brms)

################################
### Load and preprocess data ### 
################################

resistant_threshold <- 8.5812

ov <- read.delim(gzfile("../data/csv/organoid_sizes.csv.gz"),header=T,sep=',')
# remove rows with missing data
ov <- ov[complete.cases(ov), ]
ov <- ov %>% mutate(patient=ifelse(expt %in% c("exp13","exp14","exp16"),"p2","p1"))
ov$patient <- gsub("p","PDO",ov$patient)

ov_end_nk <- ov[ov$timepoint==39 & ov$treatment=="NK",]
ov_end_ctrl <- ov[ov$timepoint==39 & ov$treatment=="control",]
ov_end <- ov[ov$timepoint==39,]
ov_mid <- ov[ov$timepoint==20,]
ov_start <- ov[ov$timepoint==5,]

###############################
### Plot organoid size data ###
###############################

ggplot(ov_start[ov_start$treatment=="NK",],aes(apop_area,log(diameter_in_micrometre),group=patient)) + 
  geom_point(alpha=0.5) +
  geom_smooth(method='lm',color="black")+
  #stat_cor(method = "pearson", label.x = 3, label.y = 8.5) +
  stat_cor(method = "pearson", label.y = 4, label.x = 70) +
  geom_pointdensity() +
  scale_color_viridis()+
  ylab(expression(paste(
    "log diameter (",mu,")", sep="")))+
  xlab("Organoid apoptotic surface area (%)")+
  facet_wrap(patient ~ expt)+
  ggtitle("NK Timepoint 5")+
  coord_flip()+
  theme_pubclean()+
  theme(legend.position = "none",text = element_text(size=16),strip.background = element_rect(colour="white", fill="white"))
ggsave("../figures/Apoptosis_Size_Control_Time5.pdf")

ggplot(ov_mid[ov_mid$treatment=="NK",],aes(apop_area,log(diameter_in_micrometre),group=patient)) + 
  geom_point(alpha=0.5) +
  geom_smooth(method='lm',color="black")+
  geom_pointdensity() +
  scale_color_viridis()+
  #stat_cor(method = "pearson", label.x = 3, label.y = 8.5) +
  stat_cor(method = "pearson", label.y = 4, label.x = 80) +
  ylab(expression(paste(
    "log diameter (",mu,")", sep="")))+
  xlab("Organoid apoptotic surface area (%)")+
  facet_wrap(patient ~ expt)+
  ggtitle("NK Timepoint 20")+
  coord_flip()+
  theme_pubclean()+
  theme(legend.position = "none",text = element_text(size=16),strip.background = element_rect(colour="white", fill="white"))
ggsave("../figures/Apoptosis_Size_Control_Time20.pdf")

ggplot(ov_end_nk,aes(apop_area,log(diameter_in_micrometre),group=patient)) + 
geom_point(alpha=0.5) +
#geom_smooth(method='lm',color="darkgrey")+
#stat_cor(method = "pearson", label.x = 3, label.y = 8.5) +
stat_cor(method = "pearson", label.y = 4, label.x = 110,p.accuracy = 0.0001) +
geom_pointdensity() +
scale_color_viridis()+
scale_x_continuous(breaks = seq(0, 100, len = 5))+
scale_y_continuous(breaks = seq(3, 10, len = 8))+
facet_wrap(patient ~ expt,axis.labels="all",axes = "all")+
ggtitle("NK Treatment - Timepoint 39")+
ylab(expression(paste(
"log diameter (",mu,"m)", sep="")))+
xlab("Organoid apoptotic surface area (%)")+
coord_flip()+
theme_pubclean()+
theme(legend.position = "none",text = element_text(size=16),strip.background = element_rect(colour="white", fill="white"))
ggsave("../figures/Apoptosis_Size_NK_Time39.pdf")

ggplot(ov_end_nk,aes(apop_area,log(diameter_in_micrometre),group=patient)) + 
geom_point(alpha=0.5) +
#stat_cor(method = "pearson", label.x = 3, label.y = 8.5) +
stat_cor(method = "pearson", label.y = 4, label.x = 110,p.digits = 2, size = 3) +
geom_pointdensity() +
scale_color_viridis()+
geom_smooth(method='lm',color="tomato1",se=F,linewidth=0.7,alpha=0.1,linetype="dashed")+
scale_x_continuous(breaks = seq(0, 100, len = 5))+
scale_y_continuous(breaks = seq(3, 10, len = 8))+
facet_wrap(patient ~ expt,axis.labels="all",axes = "all")+
#facet_wrap(patient ~ expt, scales="free_y")+
ggtitle("NK Treatment - Timepoint 39")+
ylab(expression(paste(
"log diameter (",mu,"m)", sep="")))+
xlab("Organoid apoptotic surface area (%)")+
coord_flip()+
theme_pubclean()+
theme(text = element_text(size=16),legend.text=element_text(size=10),
	strip.background = element_rect(colour="white", fill="white"))+
labs(color='Organoids') 
ggsave("../figures/Apoptosis_Size_NK_Time39_withlegend.pdf")

ggplot(ov_end[ov_end$treatment=="control",],aes(apop_area,log(diameter_in_micrometre),group=patient)) + 
  geom_point(alpha=0.5) +
  #geom_smooth(method='lm',color="darkgrey")+
  #stat_cor(method = "pearson", label.x = 3, label.y = 8.5) +
  stat_cor(method = "pearson", label.y = 4, label.x = 110,p.digits = 4) +
  geom_pointdensity() +
  scale_color_viridis()+
  #facet_wrap(patient ~ expt,scales="free_y")+
  scale_x_continuous(breaks = seq(0, 100, len = 5))+
  scale_y_continuous(breaks = seq(3, 10, len = 8))+
  facet_wrap(patient ~ expt,axis.labels="all",axes = "all")+
  ggtitle("Control Treatment - Timepoint 39")+
  ylab(expression(paste(
    "log diameter (",mu,"m)", sep="")))+
  xlab("Organoid apoptotic surface area (%)")+
  coord_flip()+
  theme_pubclean()+
  theme(legend.position = "none",text = element_text(size=16),strip.background = element_rect(colour="white", fill="white"))
ggsave("../figures/Apoptosis_Size_Control_Time39.pdf")

ggplot(ov_end_nk,aes(log(diameter_in_micrometre),apop_area,group=patient)) + 
  geom_point(alpha=0.5) +
  #stat_cor(method = "pearson", label.x = 3, label.y = 8.5) +
  #stat_cor(method = "pearson", label.y = 4, label.x = 110,p.digits = 2, size = 3) +
  stat_cor(method = "pearson", label.x = 4, label.y = 110,p.digits = 2, size = 3) +
  geom_pointdensity() +
  scale_color_viridis()+
  #geom_smooth(method='lm',color="tomato1",se=F,linewidth=0.7,alpha=0.1)+
  geom_smooth(method='lm', color="tomato1",se=F,linewidth=0.7,alpha=0.1)+
  scale_y_continuous(breaks = seq(0, 100, len = 5))+
  #scale_x_continuous(breaks = seq(3, 11, len = 9))+
  xlim(3,10)+
  facet_wrap(patient ~ expt,axis.labels="all",axes = "all")+
  #facet_wrap(patient ~ expt, scales="free_y")+
  ggtitle("NK Treatment - Timepoint 39")+
  xlab(expression(paste(
    "log diameter (",mu,"m)", sep="")))+
  ylab("Organoid apoptotic surface area (%)")+
  #coord_flip()+
  theme_pubclean()+
  theme(text = element_text(size=16),legend.text=element_text(size=10),
        strip.background = element_rect(colour="white", fill="white"))+
  labs(color='Organoids') 
ggsave("../figures/Apoptosis_Size_NK_Time39_withlegend.pdf",width = 8,height=8)

ov_end$log_diameter_in_micrometre <- log(ov_end$diameter_in_micrometre)
ggplot(ov_end[ov_end$treatment=="control",],aes(log_diameter_in_micrometre,apop_area,group=patient)) + 
  geom_point(alpha=0.5) +
  #stat_cor(method = "pearson", label.x = 3, label.y = 8.5) +
  #stat_cor(data=ov_end[ov_end$treatment=="control",],method = "pearson", label.x = 4, label.y = 110,p.digits = 2,size = 3) +
  stat_cor(method = "pearson", label.x = 4, label.y = 110,p.digits = 2,size = 3) +
  geom_pointdensity() +
  scale_color_viridis()+
  geom_smooth(data=ov_end[ov_end$treatment=="control",],method='lm',color="tomato1",se=F,linewidth=0.7,alpha=0.1)+
  scale_y_continuous(breaks = seq(0, 100, len = 5))+
  #scale_x_continuous(breaks = seq(3, 10, len = 8))+
  xlim(3,10)+
  facet_wrap(patient ~ expt,axis.labels="all",axes = "all")+
  ggtitle("Control Treatment - Timepoint 39")+
  xlab(expression(paste(
    "log diameter (",mu,"m)", sep="")))+
  ylab("Organoid apoptotic surface area (%)")+
  #coord_flip()+
  theme_pubclean()+
  theme(text = element_text(size=16),legend.text=element_text(size=12),
	strip.background = element_rect(colour="white", fill="white"))+ 
  labs(color='Organoids') 

ggsave("../figures/Apoptosis_Size_Control_Time39_withlegend.pdf",width = 8,height=8)

##############################
### Run brms model fitting ###
##############################

mix <- mixture(gaussian,gaussian())
priors <- c(
  prior(normal(3,1), class = "Intercept", dpar = mu1),
  prior(normal(5,1), class = "Intercept", dpar = mu2)
)

#fit_size <- brm(bf(log(diameter_in_micrometre) ~ 1),
#                              data=ov_end_nk, family = mix,init = 0, 
#                              control = list(adapt_delta = 0.97,max_treedepth = 12),
#                              prior = priors, chains = 4, cores = 4,threads = 2,iter = 4000)

#summary(fit_size)
#saveRDS(fit_size, "../data/models/brms_fit_organoid_size_21012025.rds")

##############################
### Plot modelling results ###
##############################

# The model will be loaded from a pre-generated file to ensure consistency
fit_size <- readRDS("../data/models/brms_fit_organoid_size_21012025.rds")

# Annotate resistant and susceptible
fit_mix <- as.data.frame(pp_mixture(fit_size))
ov_end_nk$k1_pp <- fit_mix$`Estimate.P(K = 1 | Y)`
ov_end_nk <- ov_end_nk %>%
  mutate(component = ifelse(k1_pp>=0.5,"small organoids","large organoids"))

## Histogram of sizes

ggplot(ov_end_nk,aes(log(diameter_in_micrometre),group=component,fill = component)) + 
  geom_histogram(alpha=0.5,position="identity")+
  scale_fill_manual(values = c("darkblue","darkred"))+
  geom_vline(xintercept = 3.33,linetype="dashed",color = "darkred")+
  geom_vline(xintercept = 5.51,linetype="dashed", color = "darkblue")+
  facet_wrap(patient ~ expt,axis.labels="all",axes = "all")+
  ggtitle("NK Treatment - Timepoint 39")+
  xlab(expression(paste(
    "log diameter (",mu,"m)", sep="")))+
  theme_pubclean()+
  theme(text = element_text(size=16),legend.text=element_text(size=10),legend.title=element_blank(),
        strip.background = element_rect(colour="white", fill="white"))+
  labs(color='Organoids') 

small_diam_threshold <- log(max(ov_end_nk[ov_end_nk$component=="small organoids",]$diameter_in_micrometre))
ggsave("../figures/Organoid_size_histogram_Coculture_Time39.pdf")

###################################
### Conduct Fisher's exact test ###
###################################

# Are small organoids more resistant to therapeutics than large organoids?

#resistant_small <- dim(ov_end_nk[ov_end_nk$apop_area<resistant_threshold & ov_end_nk$component=="small organoids",])[1]
#resistant_large <- dim(ov_end_nk[ov_end_nk$apop_area<resistant_threshold & ov_end_nk$component=="large organoids",])[1]
#susceptible_small <- dim(ov_end_nk[ov_end_nk$apop_area>=resistant_threshold & ov_end_nk$component=="small organoids",])[1]
#susceptible_large <- dim(ov_end_nk[ov_end_nk$apop_area>=resistant_threshold & ov_end_nk$component=="large organoids",])[1]
#
#dat <- data.frame(
#  "resistant" = c(resistant_small, resistant_large ),
#  "susceptible" = c(susceptible_small, susceptible_large),
#  row.names = c("Small", "Large"),
#  stringsAsFactors = FALSE
#)
#dat <- data.frame(
#  "small" = c(resistant_small, susceptible_small ),
#  "large" = c(resistant_large, susceptible_large),
#  row.names = c("Resistant", "Susceptible"),
#  stringsAsFactors = FALSE
#)
#mosaicplot(dat,
#           main = "Mosaic plot",
#           color = TRUE
#)
#
#fe_test <- fisher.test(dat)
