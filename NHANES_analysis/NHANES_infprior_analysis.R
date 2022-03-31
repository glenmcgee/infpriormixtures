################################################
###     NHANES Analysis with constraints    ###
################################################
## Oct 25, 2021
## Run NHANES analysis using different constraint specifications
## do just for 3rd index
## also do the same on all 3 indices
## Edited Oct 25, 2021: trying variable selection methods with lower selection (i.e 80% prior inclusion probability)

## set to TRUE to run locally ## FALSE is on cluster
runLOCAL=FALSE 

## params
R <- 120000            ## no. of iterations
burn <- 0.40            ## percent burn-in
thin <- 20             ## thinning number
doLog <- FALSE         ## dont log transform the exposures
# swapX <- TRUE          ## permute exposures (within groups) to avoid artifacts as a result of unique correlation structure
# dir_alpha <- 10        ## alpha parameter for dirichlet (i.e. (alpha*(p1,p2,...p_k)) s.t. sum p_j =1 )
# folds <- 4             ## no. of folds for CV
jump <- 0.35           ## sd of random walk for theta* bsmim
sel <- seq(burn*R+1,R,by=thin) 
dat_names <- c("index3","all") ## different index settings
mod_names <- c("unconstrained","constrained","ordered","dirichlet","dirichlet_varsel","TEQ") ## different models to fit
prior_b0 <- 0.5 ## controls amount of variable selection: Beta(1,prior_b0)
if(mod_names=="dirichlet"){
  R <- 90000  ## dont need as many
  burn <- 0.4            ## percent burn-in
  thin <- 15             ## thinning number
  sel <- seq(burn*R+1,R,by=thin) 
}

set.seed(1000)

### load libraries
library(devtools)
devtools::install_github("glenmcgee/bsmim2") ## install package to fit BMIMs (only need to do this once)
library(bsmim2)
library(bkmr)
library(GIGrvg)
library(refund)
library(ggplot2)
library(gridExtra)
library(tidyverse)
library(splines)
library(MASS)
library(qgcomp)


#### if runLOCAL=FALSE -->set up to run on compute cluster (via slurm)
## submit job using "submitJob_NHANES_infprior_sim.sh"
## may need to change paths depending on where files are saved
#### otherwise run single analysis locally
## "studypop.csv" should be saved in the NHANES_analysis folder
if(runLOCAL==TRUE){
  path <- "Results/" ## path for results
  # n <- 300          ##  sample size
  # sd <- 0.5         ##  of errors 
  # low_corr <- 0     ## high correlation
  # dat_version <- "POS" ## positive indices
  mod_version <- 1  ##  model version
  # iter_no <- 0      ##  iteration number
  suffix <- paste0("_infprior","_n",n,"_sd0",sub("\\.","",sd),"_iter") #paste0("_analysis","_n",n) ## output file doesnt have a suffix
  nhanes <- na.omit(read.csv("studypop.csv"))
  # nhanes <- nhanes_raw
  # set.seed(1000+iter_no)
} else{
  path <- "Results/" ## path for results
  args <- commandArgs(trailingOnly = TRUE) ## collect arguments from shell script
  # n <- as.integer(args[1])            ## get sample size via slurm
  # sd <- as.integer(args[2])/10        ## sd of errors via slurm
  # low_corr <- as.integer(args[3])==1  ## 1=low correlation (index 3), 0=high correlation (index 1)
  dat_version <- dat_names[as.integer(args[1])]
  mod_version <- mod_names[as.integer(args[2])]  ## get model version via slurm
  # iter_no <- as.integer(args[4])      ## get iteration number
  # if(low_corr==TRUE){
  #   low_corr_suff <- "lowcorr"
  # }else{
  #   low_corr_suff <- ""
  # }
  # low_corr_suff <- ""
  suffix <- paste0(mod_version,"_",dat_version) ## append output file names with the iteration number, to be combined later
  print(paste0("model: ",suffix))
  nhanes <- na.omit(read.csv("studypop.csv")) ## read in complete data only 
  # set.seed(1000+iter_no) 
  # nhanes <- nhanes[sample(nrow(nhanes)),] ## randomly reorder data
}


# ##########################
# ###  Pre-process data  ###  
# ##########################
# ## drop some large outliers for simulations
# # nhanes <- nhanes[-(1:10),]
# ## draw set of covariates for use in simulations
# set.seed(0) 
# # nhanes <- nhanes[1:500,]
# nhanes <- nhanes[sample(nrow(nhanes))[1:500],]
# set.seed(2000+iter_no) ## set,seed to permute exposures
# print(2000+iter_no)
# source("NHANES_cleandat.R")
# set.seed(1000+iter_no) ## different outcomes/errors for each dataset
# resample_ids <- 1:n  ## same exposures each iteration
# dat <- prep_data_split(resample_ids)
# # resample_ids <- sample(nrow(nhanes))[1:n]
# # dat <- prep_data_split(resample_ids)
# # y <- dat$y
# # y_TEST <- dat$y_TEST

##########################
###  Pre-process data  ###  
##########################
source("NHANES_cleandat_infprior.R")
dat <- prep_data_full()
y <- dat$y



##########################
###  HANDLE TEQS       ###  
##########################
# ## Dictionary
# ## PCBs (excluding 126 and 169) were converted to 1000s in SAS
#   ### i.e. theyre all measured in picograms
# LBX105LA # PCB_105
# LBX118LA # PCB_118
# LBXPCBLA # PCB_126
# LBX156LA # PCB_156
# LBX157LA # PCB_157
# LBX167LA # PCB_167
# LBXHXCLA # PCB_169
# LBX189LA # PCB_189
# LBXTCDLA # Dioxin_2378TCDD
# LBXD01LA # Dioxin_12378PeCDD
# LBXD02LA # Dioxin_123478HxCDD
# LBXD03LA # Dioxin_123678HxCDD
# LBXD04LA # Dioxin_123789HxCDD
# LBXD05LA # Dioxin_1234678HpCDD
# LBXD07LA # Dioxin_12346789OCDD
# LBXF01LA # Furan_2378TCDF
# LBXF02LA # Furan_12378PeCDF
# LBXF03LA # Furan_23478PeCDF
# LBXF04LA # Furan_123478HxCDF
# LBXF05LA # Furan_123678HxCDF
# LBXF06LA # Furan_123789HxCDF
# LBXF07LA # Furan_234678HxCDF
# LBXF08LA # Furan_1234678HpCDF
# LBXF09LA # Furan_1234789HpCDF

# ## Dictionary
# (LBX105LA*0.00003) +  # PCB_105
# (LBX118LA*0.00003) +  # PCB_118
# (LBXPCBLA*0.1) +      # PCB_126
# (LBX156LA*0.00003) +  # PCB_156
# (LBX157LA*0.00003) +  # PCB_157
# (LBX167LA*0.00003) +  # PCB_167
# (LBXHXCLA*0.03) +     # PCB_169
# (LBX189LA*0.00003) +  # PCB_189
# (LBXTCDLA*1) +        # Dioxin_2378TCDD
# (LBXD01LA*1) +        # Dioxin_12378PeCDD
# (LBXD02LA*0.1) +      # Dioxin_123478HxCDD
# (LBXD03LA*0.1) +      # Dioxin_123678HxCDD
# (LBXD04LA*0.1) +      # Dioxin_123789HxCDD
# (LBXD05LA*0.01) +     # Dioxin_1234678HpCDD
# (LBXD07LA*0.0003) +   # Dioxin_12346789OCDD
# (LBXF01LA*0.1) +      # Furan_2378TCDF
# (LBXF02LA*0.03) +     # Furan_12378PeCDF
# (LBXF03LA*0.3) +      # Furan_23478PeCDF
# (LBXF04LA*0.1) +      # Furan_123478HxCDF
# (LBXF05LA*0.1) +      # Furan_123678HxCDF
# (LBXF06LA*0.1) +      # Furan_123789HxCDF
# (LBXF07LA*0.1) +      # Furan_234678HxCDF
# (LBXF08LA*0.01) +     # Furan_1234678HpCDF
# (LBXF09LA*0.01)       # Furan_1234789HpCDF


# ## Exposure Mixture
# rawnames <- c(
#   "LBX105LA", # PCB_105
#   "LBX118LA", # PCB_118
#   "LBXPCBLA", # PCB_126
#   "LBX156LA", # PCB_156
#   "LBX157LA", # PCB_157
#   "LBX167LA", # PCB_167
#   "LBXHXCLA", # PCB_169
#   "LBX189LA", # PCB_189
#   "LBXTCDLA", # Dioxin_2378TCDD
#   "LBXD01LA", # Dioxin_12378PeCDD
#   "LBXD02LA", # Dioxin_123478HxCDD
#   "LBXD03LA", # Dioxin_123678HxCDD
#   "LBXD04LA", # Dioxin_123789HxCDD
#   "LBXD05LA", # Dioxin_1234678HpCDD
#   "LBXD07LA", # Dioxin_12346789OCDD
#   "LBXF01LA", # Furan_2378TCDF
#   "LBXF02LA", # Furan_12378PeCDF
#   "LBXF03LA", # Furan_23478PeCDF
#   "LBXF04LA", # Furan_123478HxCDF
#   "LBXF05LA", # Furan_123678HxCDF
#   "LBXF06LA", # Furan_123789HxCDF
#   "LBXF07LA", # Furan_234678HxCDF
#   "LBXF08LA", # Furan_1234678HpCDF
#   "LBXF09LA") # Furan_1234789HpCDF

TEQ_coefs <- c(
  # 0.00003,    # PCB_105
  0.00003,    # PCB_118
  # 0.1,        # PCB_126
  # 0.00003,    # PCB_156
  # 0.00003,    # PCB_157
  # 0.00003,    # PCB_167
  # 0.03,       # PCB_169
  # 0.00003,    # PCB_189
  # 1,          # Dioxin_2378TCDD
  # 1,          # Dioxin_12378PeCDD
  # 0.1,        # Dioxin_123478HxCDD
  0.1,        # Dioxin_123678HxCDD
  # 0.1,        # Dioxin_123789HxCDD
  0.01,       # Dioxin_1234678HpCDD
  0.0003,     # Dioxin_12346789OCDD
  # 0.1,        # Furan_2378TCDF
  # 0.03,       # Furan_12378PeCDF
  0.3,        # Furan_23478PeCDF
  0.1,        # Furan_123478HxCDF
  0.1,        # Furan_123678HxCDF
  # 0.1,        # Furan_123789HxCDF
  # 0.1,        # Furan_234678HxCDF
  0.01 #,       # Furan_1234678HpCDF
  # 0.01       # Furan_1234789HpCDF
)


TEQ_wts <- TEQ_coefs*xsd3 ## scaling by standard deviation to get appropriate TEQ weights
TEQ_index <- dat$bsmim$X[[3]]%*%TEQ_wts ## construct fixed TEQ index
w_prior <- TEQ_wts/sum(TEQ_wts) ## standardize for dirichlet (make sum to 1 so that it's easy to manipulate)


# w_prior <- c(0.00003, ##PCB118
#              0.1, ##1,2,3,6,7,8-hxcdd Lipid Adj (pg/g)
#              0.01, ##1,2,3,4,6,7,8-hpcdd Lipid Adj (pg/g)
#              0.0003, ##1,2,3,4,6,7,8,9-ocdd Lipid Adj (pg/g)
#              0.3, ##2,3,4,7,8-pncdf Lipid Adj (pg/g)
#              0.1, ##1,2,3,4,7,8-hxcdf Lipid Adj (pg/g)
#              0.1, ##1,2,3,6,7,8-hxcdf Lipid Adj (pg/g)
#              0.01 ##1,2,3,4,6,7,8-hxcdf Lipid Adj (pg/g)
# )



################################################
###             Format Exposures             ###
################################################
## choose exposure data (and constraints length)
if(dat_version=="index3"){  ## use only index 3 (8 exposures)
  X_list <- list(dat$bsmim$X[[3]])
  cnstrs <- NULL  #(only 1 index)
  prior_alphas_list <- list(10*w_prior) # for dirichlet
  prior_slabpos_list <- list(5*w_prior) # for dirichlet+varsel
  if(mod_version=="TEQ"){
    X_list[[1]] <- as.matrix(TEQ_index)
  }
}else{                      ## all 3 indices (18 exposures)
  X_list <- dat$bsmim$X
  cnstrs <- c(0,0)  #(only constraining 3rd index)
  prior_alphas_list <- list(NA,NA,10*w_prior) # for dirichlet
  prior_slabpos_list <- list(NA,NA,5*w_prior) # for dirichlet+varsel
  if(mod_version=="TEQ"){
    X_list[[3]] <- as.matrix(TEQ_index)
  }
}

################################################
###         Fit Models                       ###
################################################

if(mod_version=="unconstrained"){ ## fit unconstrained SIM on index 3
  
  fit <- bsmim2(y=y,x=X_list,z=dat$covariates,niter=R,nburn=R*burn,nthin=thin,prior_sigma=c(0.001,0.001),prior_lambda_shaperate=c(1,0.1),gaussian=TRUE,spike_slab=TRUE,gauss_prior=TRUE,prior_theta_slab_sd=0.25,stepsize_theta=jump,basis.opts=NULL,draw_h=FALSE,prior_pi=c(1,prior_b0))
  
}else if(mod_version=="constrained"){ ## fit constrained SIM on index 3
  
  fit <- bsmim2(y=y,x=X_list,z=dat$covariates,niter=R,nburn=R*burn,nthin=thin,prior_sigma=c(0.001,0.001),prior_lambda_shaperate=c(1,0.1),gaussian=TRUE,spike_slab=TRUE,gauss_prior=TRUE,prior_theta_slab_sd=0.25,stepsize_theta=jump,basis.opts=NULL,draw_h=FALSE,constraints=c(cnstrs,1),prior_slabpos=c(1.6,8),prior_pi=c(1,prior_b0)) ## shape and rate for gamma on thetastar
  
}else if(mod_version=="ordered"){ ## not currently running (NOTE: exposures not in correct order)
  
  # fit <- bsmim2(y=y,x=X_list,z=dat$covariates,niter=R,nburn=R*burn,nthin=thin,prior_sigma=c(0.001,0.001),prior_lambda_shaperate=c(1,0.1),gaussian=TRUE,spike_slab=TRUE,gauss_prior=TRUE,prior_theta_slab_sd=0.25,stepsize_theta=jump,draw_h=FALSE,constraints=c(cnstrs,1),prior_slabpos=c(1.6,8),prior_pi=c(1,prior_b0),basis.opts.list=list(list(type="RANKED"))) ## shape and rate for gamma on thetastar
  
}else if(mod_version=="dirichlet"){ ## fit dirichlet SIM on index 3

  fit <- bsmim2(y=y,x=X_list,z=dat$covariates,niter=R,nburn=R*burn,nthin=thin,prior_sigma=c(0.001,0.001),prior_lambda_shaperate=c(1,0.1),gaussian=TRUE,spike_slab=TRUE,gauss_prior=TRUE,prior_theta_slab_sd=0.25,stepsize_theta=jump,basis.opts=NULL,draw_h=FALSE,constraints=c(cnstrs,2),prior_alphas=prior_alphas_list,prior_slabrho=c(1,1)) 
  
}else if(mod_version=="dirichlet_varsel"){ ## fit SIM with dirichlet AND variable selection via informative gammas
  
  fit <- bsmim2(y=y,x=X_list,z=dat$covariates,niter=R,nburn=R*burn,nthin=thin,prior_sigma=c(0.001,0.001),prior_lambda_shaperate=c(1,0.1),gaussian=TRUE,spike_slab=TRUE,gauss_prior=TRUE,prior_theta_slab_sd=0.25,stepsize_theta=jump,basis.opts=NULL,draw_h=FALSE,constraints=c(cnstrs,1),prior_slabpos=c(1.6,8),prior_slabpos_shape_inf=prior_slabpos_list,prior_pi=c(1,prior_b0)) # old version had ,prior_pi=c(1,1) (uniform. but investigating priors shows less variable selecion makes for more reasonably shaped priors) 
  
}else if(mod_version=="TEQ"){ ## fit SIM fixed TEQ weights
  
  fit <- bsmim2(y=y,x=X_list,z=dat$covariates,niter=R,nburn=R*burn,nthin=thin,prior_sigma=c(0.001,0.001),prior_lambda_shaperate=c(1,0.1),gaussian=TRUE,spike_slab=TRUE,gauss_prior=TRUE,prior_theta_slab_sd=0.25,stepsize_theta=jump,basis.opts=NULL,draw_h=FALSE,prior_pi=c(1,prior_b0)) 
  
}


pred_assoc <- predict_hnew_assoc2(fit)
pred_overall <- predict_hnew_assoc2(fit,overall = TRUE)
pred_ind <- predict_hnew_indexwise2(fit)
if(dat_version=="index3"){
  pred_inter <- NA
}else{
  pred_inter <- pred_twoway(fit)
}
cv <- bsmim_crossval2(fit,kfolds=4) ## 4 fold cross validation

res_list <- list(fit=fit,y=y,pred_assoc=pred_assoc,pred_overall=pred_overall,pred_ind=pred_ind,pred_inter=pred_inter,cv=cv)
save(res_list, file=paste0(path,"res_list_",suffix,".RData") )







