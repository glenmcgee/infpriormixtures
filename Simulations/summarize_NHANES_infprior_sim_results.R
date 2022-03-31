# ## code to run simulations
# sbatch --array=1-500 submitJob_NHANES_infprior_sim.sh 100 5 1 1 1
# sbatch --array=1-500 submitJob_NHANES_infprior_sim.sh 100 5 1 1 2
# sbatch --array=1-500 submitJob_NHANES_infprior_sim.sh 100 5 1 1 3
# sbatch --array=1-500 submitJob_NHANES_infprior_sim.sh 100 5 1 1 4
# sbatch --array=1-500 submitJob_NHANES_infprior_sim.sh 100 5 1 1 5
# sbatch --array=1-500 submitJob_NHANES_infprior_sim.sh 100 5 1 1 6
# sbatch --array=1-500 submitJob_NHANES_infprior_sim.sh 100 5 1 1 7
# sbatch --array=1-500 submitJob_NHANES_infprior_sim.sh 100 5 1 1 8
# 
# sbatch --array=1-500 submitJob_NHANES_infprior_sim.sh 100 5 1 2 1
# sbatch --array=1-500 submitJob_NHANES_infprior_sim.sh 100 5 1 2 2
# sbatch --array=1-500 submitJob_NHANES_infprior_sim.sh 100 5 1 2 3
# sbatch --array=1-500 submitJob_NHANES_infprior_sim.sh 100 5 1 2 4
# sbatch --array=1-500 submitJob_NHANES_infprior_sim.sh 100 5 1 2 5
# sbatch --array=1-500 submitJob_NHANES_infprior_sim.sh 100 5 1 2 6
# sbatch --array=1-500 submitJob_NHANES_infprior_sim.sh 100 5 1 2 7
# sbatch --array=1-500 submitJob_NHANES_infprior_sim.sh 100 5 1 2 8
# 
# sbatch --array=1-500 submitJob_NHANES_infprior_sim.sh 100 5 1 3 1
# sbatch --array=1-500 submitJob_NHANES_infprior_sim.sh 100 5 1 3 2
# sbatch --array=1-500 submitJob_NHANES_infprior_sim.sh 100 5 1 3 3
# sbatch --array=1-500 submitJob_NHANES_infprior_sim.sh 100 5 1 3 4
# sbatch --array=1-500 submitJob_NHANES_infprior_sim.sh 100 5 1 3 5
# sbatch --array=1-500 submitJob_NHANES_infprior_sim.sh 100 5 1 3 6
# sbatch --array=1-500 submitJob_NHANES_infprior_sim.sh 100 5 1 3 7
# sbatch --array=1-500 submitJob_NHANES_infprior_sim.sh 100 5 1 3 8

## code to concatenate results
# source("combineResults.R")
# combine_dat("bkmr_infprior_n100_POSlowcorr_sd005_iter*")
# combine_dat("SIM_infprior_n100_POSlowcorr_sd005_iter*")
# combine_dat("SIMpos_infprior_n100_POSlowcorr_sd005_iter*")
# combine_dat("SIMdirequal_infprior_n100_POSlowcorr_sd005_iter*")
# combine_dat("SIMdir_infprior_n100_POSlowcorr_sd005_iter*")
# combine_dat("SIMdirvarsel_infprior_n100_POSlowcorr_sd005_iter*")
# combine_dat("SIMrank_infprior_n100_POSlowcorr_sd005_iter*")
# combine_dat("TEQ_infprior_n100_POSlowcorr_sd005_iter*")
# 
# 
# combine_dat("bkmr_infprior_n100_POSwronglowcorr_sd005_iter*")
# combine_dat("SIM_infprior_n100_POSwronglowcorr_sd005_iter*")
# combine_dat("SIMpos_infprior_n100_POSwronglowcorr_sd005_iter*")
# combine_dat("SIMdirequal_infprior_n100_POSwronglowcorr_sd005_iter*")
# combine_dat("SIMdir_infprior_n100_POSwronglowcorr_sd005_iter*")
# combine_dat("SIMdirvarsel_infprior_n100_POSwronglowcorr_sd005_iter*")
# combine_dat("SIMrank_infprior_n100_POSwronglowcorr_sd005_iter*")
# combine_dat("TEQ_infprior_n100_POSwronglowcorr_sd005_iter*")
# 
# 
# combine_dat("bkmr_infprior_n100_NEGlowcorr_sd005_iter*")
# combine_dat("SIM_infprior_n100_NEGlowcorr_sd005_iter*")
# combine_dat("SIMpos_infprior_n100_NEGlowcorr_sd005_iter*")
# combine_dat("SIMdirequal_infprior_n100_NEGlowcorr_sd005_iter*")
# combine_dat("SIMdir_infprior_n100_NEGlowcorr_sd005_iter*")
# combine_dat("SIMdirvarsel_infprior_n100_NEGlowcorr_sd005_iter*")
# combine_dat("SIMrank_infprior_n100_NEGlowcorr_sd005_iter*")
# combine_dat("TEQ_infprior_n100_NEGlowcorr_sd005_iter*")
# 


## Code to summarize results of NHANES simulations

library(tidyverse)
loadpath=""

## Simulation A

bkmr_infprior_n100_POSlowcorr_sd005 <- read_csv(paste0(loadpath,"Results/bkmr_infprior_n100_POSlowcorr_sd005.csv"))
SIM_infprior_n100_POSlowcorr_sd005 <- read_csv(paste0(loadpath,"Results/SIM_infprior_n100_POSlowcorr_sd005.csv"))
SIMpos_infprior_n100_POSlowcorr_sd005 <- read_csv(paste0(loadpath,"Results/SIMpos_infprior_n100_POSlowcorr_sd005.csv"))
SIMdirequal5_infprior_n100_POSlowcorr_sd005 <- read_csv(paste0(loadpath,"Results/SIMdirequal5_infprior_n100_POSlowcorr_sd005.csv"))
SIMdir5_infprior_n100_POSlowcorr_sd005 <- read_csv(paste0(loadpath,"Results/SIMdir5_infprior_n100_POSlowcorr_sd005.csv"))
SIMdirvarsel5_infprior_n100_POSlowcorr_sd005 <- read_csv(paste0(loadpath,"Results/SIMdirvarsel5_infprior_n100_POSlowcorr_sd005.csv"))
SIMdirequal10_infprior_n100_POSlowcorr_sd005 <- read_csv(paste0(loadpath,"Results/SIMdirequal10_infprior_n100_POSlowcorr_sd005.csv"))
SIMdir10_infprior_n100_POSlowcorr_sd005 <- read_csv(paste0(loadpath,"Results/SIMdir10_infprior_n100_POSlowcorr_sd005.csv"))
SIMdirvarsel10_infprior_n100_POSlowcorr_sd005 <- read_csv(paste0(loadpath,"Results/SIMdirvarsel10_infprior_n100_POSlowcorr_sd005.csv"))
SIMrank_infprior_n100_POSlowcorr_sd005 <- read_csv(paste0(loadpath,"Results/SIMrank_infprior_n100_POSlowcorr_sd005.csv"))
SIMrank24_infprior_n100_POSlowcorr_sd005 <- read_csv(paste0(loadpath,"Results/SIMrank24_infprior_n100_POSlowcorr_sd005.csv"))
TEQ_infprior_n100_POSlowcorr_sd005 <- read_csv(paste0(loadpath,"Results/TEQ_infprior_n100_POSlowcorr_sd005.csv"))
## setting c=100 (very high certainty with variable selection)
SIMdirvarsel100_infprior_n100_POSlowcorr_sd005 <- read_csv(paste0(loadpath,"Results/SIMdirvarsel100_infprior_n100_POSlowcorr_sd005.csv"))


rawA <- rbind(
  apply(bkmr_infprior_n100_POSlowcorr_sd005,2,mean)[2:13],
  apply(SIM_infprior_n100_POSlowcorr_sd005,2,mean)[2:13],
  apply(SIMpos_infprior_n100_POSlowcorr_sd005,2,mean)[2:13],
  apply(SIMdirequal5_infprior_n100_POSlowcorr_sd005,2,mean)[2:13],
  apply(SIMdir5_infprior_n100_POSlowcorr_sd005,2,mean)[2:13],
  apply(SIMdirvarsel5_infprior_n100_POSlowcorr_sd005,2,mean)[2:13],
  # apply(SIMdirequal10_infprior_n100_POSlowcorr_sd005,2,mean)[2:13],
  # apply(SIMdir10_infprior_n100_POSlowcorr_sd005,2,mean)[2:13],
  # apply(SIMdirvarsel10_infprior_n100_POSlowcorr_sd005,2,mean)[2:13], ## less certainty with variable selection
  apply(SIMrank_infprior_n100_POSlowcorr_sd005,2,mean)[2:13],
  # apply(SIMrank24_infprior_n100_POSlowcorr_sd005,2,mean)[2:13],
  apply(TEQ_infprior_n100_POSlowcorr_sd005,2,mean)[2:13]  )

# round(rawA,2)
resA <- rawA/matrix(rawA[2,],ncol=12,nrow=8,byrow=T)
resA[,c(3,6,9,12)] <- rawA[,c(3,6,9,12)]
round(resA,2)



## look at distribution of posterior means
datalist <- list(
  SIMpos_infprior_n100_POSlowcorr_sd005,
  SIMdirvarsel5_infprior_n100_POSlowcorr_sd005,
  SIMdir5_infprior_n100_POSlowcorr_sd005,
  SIMrank_infprior_n100_POSlowcorr_sd005,
  TEQ_infprior_n100_POSlowcorr_sd005
)
namelist <- c(#"(i) Unconstrained",
              "(ii) Constrained",
              "(iii) Dirichlet",
              "(iv) Dir (No selection)",
              "(v) Ranked",
              "(vi) TEQ")
longdf <- c()
for(ll in 1:length(datalist)){
  df_wide <- datalist[[ll]]
  df_wide <- data.frame(df_wide[,grepl("weight_avg", colnames(df_wide), fixed = TRUE)])
  colnames(df_wide) <- as.character(1:8)
  df <- gather(df_wide,component,weight)
  df$Model <- namelist[ll]#ll
  longdf <- rbind(longdf,df)
}
longdf$Model <- as.factor(longdf$Model)
#
weights_box_A <- ggplot(data=longdf,aes(y=weight,x=component,fill=Model))+
  geom_boxplot( outlier.alpha = 0.25,outlier.size = 0.25,lwd=0.25)+
  # labs(title=expression("Informative Dirichlet ("*L[m]*"=3)"),x="Component", y = "Weight")+
  labs(title=paste("Simulation A"),x=expression("Component "*(italic(l))), y = expression("Weight "*(italic(w[ml]))))+
  ylim(c(0,1))+
  theme_classic()
wA <- c(0.50, 0.25, 0.10, 0.05, 0.05, 0.02, 0.02, 0.01) 
for(j in 1:8){
  weights_box_A <- weights_box_A+geom_segment(size=0.25,aes_(x=j-0.5,xend=j+0.5,y=wA[j],yend=wA[j]))
}
weights_box_A
ggsave(paste0(loadpath,"Plots/weights_box_POS.pdf"),weights_box_A,width=6,height=4)


###### Supplementary plot investigating high certainty with variable selection
## look at distribution of posterior means
datalist <- list(
  SIMdirvarsel5_infprior_n100_POSlowcorr_sd005,
  SIMdirvarsel10_infprior_n100_POSlowcorr_sd005,
  SIMdirvarsel100_infprior_n100_POSlowcorr_sd005,
  TEQ_infprior_n100_POSlowcorr_sd005
)
namelist <- c(#"(i) Unconstrained",
  "(iv-a) Dirichlet-5",
  "(iv-b) Dirichlet-10",
  "(iv-c) Dirichlet-100",
  "(vi) TEQ")
longdf <- c()
for(ll in 1:length(datalist)){
  df_wide <- datalist[[ll]]
  df_wide <- data.frame(df_wide[,grepl("weight_avg", colnames(df_wide), fixed = TRUE)])
  colnames(df_wide) <- as.character(1:8)
  df <- gather(df_wide,component,weight)
  df$Model <- namelist[ll]#ll
  longdf <- rbind(longdf,df)
}
longdf$Model <- as.factor(longdf$Model)
#
weights_box_A100 <- ggplot(data=longdf,aes(y=weight,x=component,fill=Model))+
  geom_boxplot( outlier.alpha = 0.25,outlier.size = 0.25,lwd=0.25)+
  # labs(title=expression("Informative Dirichlet ("*L[m]*"=3)"),x="Component", y = "Weight")+
  labs(title=paste("Simulation A"),x=expression("Component "*(italic(l))), y = expression("Weight "*(italic(w[ml]))))+
  ylim(c(0,1))+
  theme_classic()
wA <- c(0.50, 0.25, 0.10, 0.05, 0.05, 0.02, 0.02, 0.01) 
for(j in 1:8){
  weights_box_A100 <- weights_box_A100+geom_segment(size=0.25,aes_(x=j-0.5,xend=j+0.5,y=wA[j],yend=wA[j]))
}
weights_box_A100
ggsave(paste0(loadpath,"Plots/weights100_box_POS.pdf"),weights_box_A100,width=6,height=4)






##############################################################
## Simulation B

bkmr_infprior_n100_POSwronglowcorr_sd005 <- read_csv(paste0(loadpath,"Results/bkmr_infprior_n100_POSwronglowcorr_sd005.csv"))
SIM_infprior_n100_POSwronglowcorr_sd005 <- read_csv(paste0(loadpath,"Results/SIM_infprior_n100_POSwronglowcorr_sd005.csv"))
SIMpos_infprior_n100_POSwronglowcorr_sd005 <- read_csv(paste0(loadpath,"Results/SIMpos_infprior_n100_POSwronglowcorr_sd005.csv"))
SIMdirequal5_infprior_n100_POSwronglowcorr_sd005 <- read_csv(paste0(loadpath,"Results/SIMdirequal5_infprior_n100_POSwronglowcorr_sd005.csv"))
SIMdir5_infprior_n100_POSwronglowcorr_sd005 <- read_csv(paste0(loadpath,"Results/SIMdir5_infprior_n100_POSwronglowcorr_sd005.csv"))
SIMdirvarsel5_infprior_n100_POSwronglowcorr_sd005 <- read_csv(paste0(loadpath,"Results/SIMdirvarsel5_infprior_n100_POSwronglowcorr_sd005.csv"))
SIMdirequal10_infprior_n100_POSwronglowcorr_sd005 <- read_csv(paste0(loadpath,"Results/SIMdirequal10_infprior_n100_POSwronglowcorr_sd005.csv"))
SIMdir10_infprior_n100_POSwronglowcorr_sd005 <- read_csv(paste0(loadpath,"Results/SIMdir10_infprior_n100_POSwronglowcorr_sd005.csv"))
SIMdirvarsel10_infprior_n100_POSwronglowcorr_sd005 <- read_csv(paste0(loadpath,"Results/SIMdirvarsel10_infprior_n100_POSwronglowcorr_sd005.csv"))
SIMrank_infprior_n100_POSwronglowcorr_sd005 <- read_csv(paste0(loadpath,"Results/SIMrank_infprior_n100_POSwronglowcorr_sd005.csv"))
TEQ_infprior_n100_POSwronglowcorr_sd005 <- read_csv(paste0(loadpath,"Results/TEQ_infprior_n100_POSwronglowcorr_sd005.csv"))
## setting c=100 (very high certainty with variable selection)
SIMdirvarsel100_infprior_n100_POSwronglowcorr_sd005 <- read_csv(paste0(loadpath,"Results/SIMdirvarsel100_infprior_n100_POSwronglowcorr_sd005.csv"))

rawB <- rbind(
  apply(bkmr_infprior_n100_POSwronglowcorr_sd005,2,mean)[2:13],
  apply(SIM_infprior_n100_POSwronglowcorr_sd005,2,mean)[2:13],
  apply(SIMpos_infprior_n100_POSwronglowcorr_sd005,2,mean)[2:13],
  apply(SIMdirequal5_infprior_n100_POSwronglowcorr_sd005,2,mean)[2:13],
  apply(SIMdir5_infprior_n100_POSwronglowcorr_sd005,2,mean)[2:13],
  apply(SIMdirvarsel5_infprior_n100_POSwronglowcorr_sd005,2,mean)[2:13],
  # apply(SIMdirequal10_infprior_n100_POSwronglowcorr_sd005,2,mean)[2:13],
  # apply(SIMdir10_infprior_n100_POSwronglowcorr_sd005,2,mean)[2:13],
  # apply(SIMdirvarsel10_infprior_n100_POSwronglowcorr_sd005,2,mean)[2:13],
  apply(SIMrank_infprior_n100_POSwronglowcorr_sd005,2,mean)[2:13],
  apply(TEQ_infprior_n100_POSwronglowcorr_sd005,2,mean)[2:13]  )

# round(rawB,2)
resB <- rawB/matrix(rawB[2,],ncol=12,nrow=8,byrow=T)
resB[,c(3,6,9,12)] <- rawB[,c(3,6,9,12)]
round(resB,2)

## look at distribution of posterior means
datalist <- list(
  SIMpos_infprior_n100_POSwronglowcorr_sd005,
  SIMdirvarsel5_infprior_n100_POSwronglowcorr_sd005,
  SIMdir5_infprior_n100_POSwronglowcorr_sd005,
  SIMrank_infprior_n100_POSwronglowcorr_sd005,
  TEQ_infprior_n100_POSwronglowcorr_sd005
)
namelist <- c(#"(i) Unconstrained",
  "(ii) Constrained",
  "(iii) Dirichlet",
  "(iv) Dir (No selection)",
  "(v) Ranked",
  "(vi) TEQ")
longdf <- c()
for(ll in 1:length(datalist)){
  df_wide <- datalist[[ll]]
  df_wide <- data.frame(df_wide[,grepl("weight_avg", colnames(df_wide), fixed = TRUE)])
  colnames(df_wide) <- as.character(1:8)
  df <- gather(df_wide,component,weight)
  df$Model <- namelist[ll]#ll
  longdf <- rbind(longdf,df)
}
longdf$Model <- as.factor(longdf$Model)
#
weights_box_B <- ggplot(data=longdf,aes(y=weight,x=component,fill=Model))+
  geom_boxplot( outlier.alpha = 0.25,outlier.size = 0.25,lwd=0.25)+
  # labs(title=expression("Informative Dirichlet ("*L[m]*"=3)"),x="Component", y = "Weight")+
  labs(title=paste("Simulation B"),x=expression("Component "*(italic(l))), y = expression("Weight "*(italic(w[ml]))))+
  ylim(c(0,1))+
  theme_classic()
wB <- c(0.10, 0.25, 0.50, 0.05, 0.05, 0.02, 0.02, 0.01)
for(j in 1:8){
  weights_box_B <- weights_box_B+geom_segment(size=0.25,aes_(x=j-0.5,xend=j+0.5,y=wB[j],yend=wB[j]))
}
weights_box_B
ggsave(paste0(loadpath,"Plots/weights_box_POSwrong.pdf"),weights_box_B,width=6,height=4)










########################################################
## Simulation C

bkmr_infprior_n100_NEGlowcorr_sd005 <- read_csv(paste0(loadpath,"Results/bkmr_infprior_n100_NEGlowcorr_sd005.csv"))
SIM_infprior_n100_NEGlowcorr_sd005 <- read_csv(paste0(loadpath,"Results/SIM_infprior_n100_NEGlowcorr_sd005.csv"))
SIMpos_infprior_n100_NEGlowcorr_sd005 <- read_csv(paste0(loadpath,"Results/SIMpos_infprior_n100_NEGlowcorr_sd005.csv"))
SIMdirequal5_infprior_n100_NEGlowcorr_sd005 <- read_csv(paste0(loadpath,"Results/SIMdirequal5_infprior_n100_NEGlowcorr_sd005.csv"))
SIMdir5_infprior_n100_NEGlowcorr_sd005 <- read_csv(paste0(loadpath,"Results/SIMdir5_infprior_n100_NEGlowcorr_sd005.csv"))
SIMdirvarsel5_infprior_n100_NEGlowcorr_sd005 <- read_csv(paste0(loadpath,"Results/SIMdirvarsel5_infprior_n100_NEGlowcorr_sd005.csv"))
SIMdirequal10_infprior_n100_NEGlowcorr_sd005 <- read_csv(paste0(loadpath,"Results/SIMdirequal10_infprior_n100_NEGlowcorr_sd005.csv"))
SIMdir10_infprior_n100_NEGlowcorr_sd005 <- read_csv(paste0(loadpath,"Results/SIMdir10_infprior_n100_NEGlowcorr_sd005.csv"))
SIMdirvarsel10_infprior_n100_NEGlowcorr_sd005 <- read_csv(paste0(loadpath,"Results/SIMdirvarsel10_infprior_n100_NEGlowcorr_sd005.csv"))
SIMrank_infprior_n100_NEGlowcorr_sd005 <- read_csv(paste0(loadpath,"Results/SIMrank_infprior_n100_NEGlowcorr_sd005.csv"))
TEQ_infprior_n100_NEGlowcorr_sd005 <- read_csv(paste0(loadpath,"Results/TEQ_infprior_n100_NEGlowcorr_sd005.csv"))
## setting c=100 (very high certainty with variable selection)
SIMdirvarsel100_infprior_n100_NEGlowcorr_sd005 <- read_csv(paste0(loadpath,"Results/SIMdirvarsel100_infprior_n100_NEGlowcorr_sd005.csv"))

rawC <- rbind(
  apply(bkmr_infprior_n100_NEGlowcorr_sd005,2,mean)[2:13],
  apply(SIM_infprior_n100_NEGlowcorr_sd005,2,mean)[2:13],
  apply(SIMpos_infprior_n100_NEGlowcorr_sd005,2,mean)[2:13],
  apply(SIMdirequal5_infprior_n100_NEGlowcorr_sd005,2,mean)[2:13],
  apply(SIMdir5_infprior_n100_NEGlowcorr_sd005,2,mean)[2:13],
  apply(SIMdirvarsel5_infprior_n100_NEGlowcorr_sd005,2,mean)[2:13],
  # apply(SIMdirequal10_infprior_n100_NEGlowcorr_sd005,2,mean)[2:13],
  # apply(SIMdir10_infprior_n100_NEGlowcorr_sd005,2,mean)[2:13],
  # apply(SIMdirvarsel10_infprior_n100_NEGlowcorr_sd005,2,mean)[2:13],
  apply(SIMrank_infprior_n100_NEGlowcorr_sd005,2,mean)[2:13],
  apply(TEQ_infprior_n100_NEGlowcorr_sd005,2,mean)[2:13]  )

# round(rawC,2)
resC <- rawC/matrix(rawC[2,],ncol=12,nrow=8,byrow=T)
resC[,c(3,6,9,12)] <- rawC[,c(3,6,9,12)]
round(resC,2)


## look at distribution of posterior means
datalist <- list(
  SIMpos_infprior_n100_NEGlowcorr_sd005,
  SIMdirvarsel5_infprior_n100_NEGlowcorr_sd005,
  SIMdir5_infprior_n100_NEGlowcorr_sd005,
  SIMrank_infprior_n100_NEGlowcorr_sd005,
  TEQ_infprior_n100_NEGlowcorr_sd005
)
namelist <- c(#"(i) Unconstrained",
  "(ii) Constrained",
  "(iii) Dirichlet",
  "(iv) Dir (No selection)",
  "(v) Ranked",
  "(vi) TEQ")
longdf <- c()
for(ll in 1:length(datalist)){
  df_wide <- datalist[[ll]]
  df_wide <- data.frame(df_wide[,grepl("weight_avg", colnames(df_wide), fixed = TRUE)])
  colnames(df_wide) <- as.character(1:8)
  df <- gather(df_wide,component,weight)
  df$Model <- namelist[ll]#ll
  longdf <- rbind(longdf,df)
}
longdf$Model <- as.factor(longdf$Model)
#
weights_box_C <- ggplot(data=longdf,aes(y=weight,x=component,fill=Model))+
  geom_boxplot( outlier.alpha = 0.25,outlier.size = 0.25,lwd=0.25)+
  # labs(title=expression("Informative Dirichlet ("*L[m]*"=3)"),x="Component", y = "Weight")+
  labs(title=paste("Simulation C"),x=expression("Component "*(italic(l))), y = expression("Weight "*(italic(w[ml]))))+
  ylim(c(0,1))+
  theme_classic()
# wC <- c(0.50, -0.25, 0.10, 0.05, 0.05, 0.02, 0.02, 0.01)
# for(j in 1:8){
#   weights_box_C <- weights_box_C+geom_segment(size=0.5,aes_(x=j-0.5,xend=j+0.5,y=wC[j],yend=wC[j]))
# }
weights_box_C
ggsave(paste0(loadpath,"Plots/weights_box_NEG.pdf"),weights_box_C,width=6,height=4)



## Priors

## look at prior
priorSIMpos<- investigate_priors(8,4000,constraint=1,prior_slabpos=c(1.6,8))
priorSIMdir5<- investigate_priors(8,4000,constraint=2,prior_slabpos=c(1.6,8),prior_alphas=5*wA,prior_slabrho=c(1,1))
priorSIMdirequal5<- investigate_priors(8,4000,constraint=2,prior_slabpos=c(1.6,8),prior_alphas=5*rep(1/8,8),prior_slabrho=c(1,1))
priorSIMdirvarsel5<- investigate_priors(8,4000,constraint=1,prior_slabpos=c(1.6,8),prior_slabpos_shape_inf=(5*wA))
priorSIMdir10<- investigate_priors(8,4000,constraint=2,prior_slabpos=c(1.6,8),prior_alphas=10*wA,prior_slabrho=c(1,1))
priorSIMdirequal10<- investigate_priors(8,4000,constraint=2,prior_slabpos=c(1.6,8),prior_alphas=10*rep(1/8,8),prior_slabrho=c(1,1))
priorSIMdirvarsel10<- investigate_priors(8,4000,constraint=1,prior_slabpos=c(1.6,8),prior_slabpos_shape_inf=(10*wA))
priorSIMrank <- investigate_priors(8,4000,constraint=1,prior_slabpos=c(3.2,24)) 
priorSIMrank$thetaPOS <- priorSIMpos$thetaPOS%*%t(getbasis(diag(8),basis.opts=list(type="RANKED"))$psi)
priorSIMrank$thetaPOS <- priorSIMrank$thetaPOS/apply(priorSIMrank$thetaPOS,1,sum)
priorTEQ <- list(thetaPOS=matrix(wA,byrow=TRUE,nrow=10,ncol=8))

prior_list <- list(priorSIMpos,
                   # priorSIMdirequal5,
                   priorSIMdirvarsel5,
                   priorSIMdir5,
                   # priorSIMdirequal10,
                   # priorSIMdirvarsel10,
                   # priorSIMdir10,
                   priorSIMrank)
                   #priorTEQ)
namelist <- c(#"(i) Unconstrained",
  "(ii) Constrained",
  "(iii) Dirichlet (5)",
  # "(iii) Dirichlet (10)",
  "(iv) Dir (No selection; 5)",
  # "(iv) Dir (No selection; 10)",
  "(v) Ranked")
  #"(vi) TEQ")
prior_plots <- list()
for(ll in 1:length(prior_list)){
  priordf <- data.frame(prior_list[[ll]]$thetaPOS); colnames(priordf) <- as.character(1:8)
  #priordf <- priordf[apply(priordf,1,function(x) sum(x==0))==0,]
  priordf <- gather(priordf,component,weight)
  prior_plots[[ll]] <- ggplot(data=priordf,aes(y=weight,x=factor(component)))+
    geom_violin(trim=TRUE,fill="gray", scale = "width")+
    # labs(title=expression("Informative Dirichlet ("*L[m]*"=3)"),x="Component", y = "Weight")+
    labs(title=namelist[ll],x=expression("Component "*(italic(l))), y = expression("Weight "*(italic(w[ml]))))+
    ylim(c(0,1))+
    theme_classic()
}

## thetastar
prior_plotsstar <- list()
for(ll in 1:length(prior_list)){
  priordf <- prior_list[[ll]]$thetastar
  if(namelist[ll]=="(v) Ranked"){ ## apply transformation for ranked method
    priordf <- data.frame(priordf%*%lower.tri(diag(8),diag=TRUE)+0)
  }
  priordf <- data.frame(priordf)
  colnames(priordf) <- as.character(1:8)
  #priordf <- priordf[apply(priordf,1,function(x) sum(x==0))==0,]

  priordf <- gather(priordf,component,weight)
  prior_plotsstar[[ll]] <- ggplot(data=priordf,aes(y=weight,x=factor(component)))+
    geom_violin(trim=TRUE,fill="gray", scale = "width")+
    # labs(title=expression("Informative Dirichlet ("*L[m]*"=3)"),x="Component", y = "Weight")+
    labs(title=namelist[ll],x=expression("Component "*(italic(l))), y = expression((italic(theta[ml]))) )+
    ylim(c(0,2))+
    theme_classic()
}
prior_plotsstar[[4]]

