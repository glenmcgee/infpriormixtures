#################################
###        Simulations        ###
#################################
## Run informative prior simulations
### Generate data from SIM
### Using only 8 exposures (M=1, L_1=8)
### Fit under 6 different SIM specifications (and bkmr for comparison)
#### 1 Unconstrained +varsel
#### 2 Constrained +varsel
#### 3 Dirichlet +varsel
#### 4 Dirichlet 
#### 5 Rank-ordered + varsel
#### 6 TEQ 
## Allows 3 different true index weight structures
### All positive, in decreasing ordered
### All positive, incorrect order
### Some negative, some positive

## set to TRUE to run locally ## FALSE is on cluster
runLOCAL=FALSE 

## params
R <- 100000            ## no. of iterations
burn <- 0.4            ## percent burn-in
thin <- 40             ## thinning number
doLog <- FALSE         ## dont log transform the exposures
swapX <- TRUE          ## permute exposures (within groups) to avoid artifacts as a result of unique correlation structure
folds <- 4             ## no. of folds for CV
jump <- 0.35           ## sd of random walk for theta* bsmim
sel <- seq(burn*R+1,R,by=thin) 
dat_names <- c("POS","POSwrong","NEG") ## different index settings
mod_names <- c("bkmr","SIM","SIMpos","SIMrank","SIMdir","SIMdirequal","SIMdirvarsel","TEQ") ## different models to fit
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



#### if runLOCAL=FALSE -->set up for batch jobs on compute cluster (via slurm)
## submit batch jobs using "submitJob_NHANES_infprior_sim.sh"
## combine results via function in "combineResults.R"
## may need to change paths depending on where files are saved
#### otherwise run single simulation locally
## "studypop.csv" should be saved in the NHANES_analysis folder
if(runLOCAL==TRUE){
  path <- "Results/" ## path for results
  n <- 300          ##  sample size
  sd <- 0.5         ##  of errors 
  low_corr <- 0     ## high correlation
  dat_version <- "POS" ## positive indices
  mod_version <- 1  ##  model version
  iter_no <- 0      ##  iteration number
  suffix <- paste0("_infprior","_n",n,"_sd0",sub("\\.","",sd),"_iter") 
  nhanes <- na.omit(read.csv("../NHANES_analysis/studypop.csv"))

} else{
  ## arguments passed from slurm for batch jobs on cluster
  path <- "Results/" ## path for results
  args <- commandArgs(trailingOnly = TRUE) ## collect arguments from shell script
  n <- as.integer(args[1])            ## get sample size via slurm
  sd <- as.integer(args[2])/10        ## sd of errors via slurm
  low_corr <- as.integer(args[3])==1  ## 1=low correlation (index 3), 0=high correlation (index 1)
  dat_version <- dat_names[as.integer(args[4])]
  mod_version <- mod_names[as.integer(args[5])]  ## get model version via slurm
  iter_no <- as.integer(args[6])      ## get iteration number
  if(low_corr==TRUE){
    low_corr_suff <- "lowcorr"
  }else{
    low_corr_suff <- ""
  }
  suffix <- paste0("_infprior","_n",n,"_",dat_version,low_corr_suff,"_sd0",sub("\\.","",sd),"_iter",iter_no) ## append output file names with the iteration number, to be combined later
  print(paste0("mod",mod_version,suffix))
  nhanes <- na.omit(read.csv("studypop.csv")) ## read in complete data only 

}


##########################
###  Pre-process data  ###  
##########################
## drop some large outliers for simulations
nhanes <- nhanes[-(1:10),]
## draw set of covariates for use in simulations
set.seed(0) 
# nhanes <- nhanes[1:500,]
nhanes <- nhanes[sample(nrow(nhanes))[1:300],]
set.seed(2000+iter_no) ## set,seed to permute exposures
print(2000+iter_no)
source("NHANES_cleandat_infprior.R")
set.seed(1000+iter_no) ## different outcomes/errors for each dataset
resample_ids <- 1:n  ## 
## standardize relative to actual sample only (and then standardize test set accordingly)
xmeans <- apply(as.matrix(X)[resample_ids,],2,mean)
xsds <- apply(as.matrix(X)[resample_ids,],2,sd)
X <- (as.matrix(X)-matrix(xmeans,byrow=TRUE,ncol=ncol(X),nrow=nrow(X)))/matrix(xsds,byrow=TRUE,ncol=ncol(X),nrow=nrow(X))
## format data
dat <- prep_data_split(resample_ids)



################################################
###         Data-Generating Mechanism        ###
################################################
wA <- c(0.50, 0.25, 0.10, 0.05, 0.05, 0.02, 0.02, 0.01) 
wB <- c(0.10, 0.25, 0.50, 0.05, 0.05, 0.02, 0.02, 0.01) ## swapped w1 and w3
wC <- c(0.50, -0.25, 0.10, 0.05, 0.05, 0.02, 0.02, 0.01) ## set w2 to be negative
## make true weight function
if(dat_version=="POS"){
  ww <- wA
}else if(dat_version=="POSwrong"){
  ww <- wB
}else if(dat_version=="NEG"){
  ww <- wC
}# else if(dat_version=="ZERO"){
#   ww <- c(4:1,0,0,0,0)
# }
wwnorm <- ww/sqrt(sum(ww^2))

if(low_corr==TRUE){ ## in low-correlation setting, use index 3 (median pairwise correlation of 0.46)
  w1 <- c(rep(0,8),rep(0,2),ww); w1 <- w1/sqrt(sum(w1^2))
}else{ ## otherwise use index 1 (median pairwise correlation of 0.86)
  w1 <- c(ww,rep(0,2),rep(0,8)); w1 <- w1/sqrt(sum(w1^2))
}


# index
mn1 <- mean(X%*%w1);sd1 <- sd(X%*%w1)

## exposure response functions
hfun1 <- function(z){ # the 1-dimensional single index function
  5*dnorm((z-mn1)/sd1)
}
hfun <- function(x){ ## function which takes vector and turns it into h
  hfun1(x%*%w1)
}
h <- hfun(X)
h_TRAIN <- h[resample_ids]
h_TEST <- h[-resample_ids]

## covariate effects 
gamma <- c(-0.43,0.00,-0.25,0.12,0.08)#,-0.06,0.23,-0.12,-0.02,0.02,0.08,-0.07,0.01,0.02,-0.04,-0.18,-0.28,-0.15)

## generate outcomes for entire dataset
yfull <- h + covariates%*%gamma + rnorm(nrow(nhanes),0,1)*sd
y <- dat$y <- yfull[resample_ids]
y_TEST <- dat$y_TEST <- yfull[-resample_ids]

#############################
###     true values       ###
#############################
## note that for overall effect contrasts, we look at an evenly spaced grid of quantiles, i.e. 25th, 30th, 35th percentiles....
## for component/index curves we take an evenly spaced grid of exposure values x between 25th and 75th percentiles (this is distinct from the overall contrast grid)

### true overall association (Note: these are contrasts)
overall_true <- hfun(apply(dat$SIM$X[[1]],2,function(x) quantile(x,seq(0.25,0.75,by=0.05)))) #overall_true <- hfun((apply(dat$SIM$X[[1]],2,function(x) quantile(x,seq(0.25,0.75,by=0.05)))%*%w1-mn1)/sd1)
overall_true <- overall_true-overall_true[6] ## subtract median

### true index curves
index_vals <- dat$SIM$X[[1]][,10*low_corr+1:8]%*%wwnorm ## 
grid_index <- seq(quantile(index_vals,0.25),quantile(index_vals,0.75),length=11)
index_true <- hfun1(grid_index)


### true component curves 
Xmedian <- apply(dat$SIM$X[[1]],2,function(x) quantile(x,0.5))
Xgrid <- matrix(Xmedian,nrow=11*8,ncol=length(Xmedian),byrow=TRUE) ## 11 points between 0.25 and 0.75, 8 exposures in our index
for (ii in 1:8){
  Xgrid[((ii-1)*11+1):(ii*11),(low_corr*10)+ii] <- seq(quantile(dat$SIM$X[[1]][,(low_corr*10)+ii],0.25),quantile(dat$SIM$X[[1]][,(low_corr*10)+ii],0.75),length=11) # old version was grid of quantiles, but now we've changed it to match what predict_hnew2() does by default #quantile(dat$SIM$X[[1]][,(low_corr*10)+ii],seq(0.25,0.75,length=11))
}
component_true <- hfun(Xgrid) #component_true <- hfun((Xgrid%*%w1-mn1)/sd1) 


##################################
###     data for fitting       ###
##################################

# consider only 8 exposures
w_prior <- wA
w_prior <- wA/sum(wA) # as proportions

# get data for TEQ model
XTEQ <- dat$SIM$X[[1]][,10*low_corr+1:8]%*%w_prior             ## for fitting
XTEQ_TEST <- dat$SIM$X_TEST[[1]][,10*low_corr+1:8]%*%w_prior   ## for predicting new h
XTEQ_component <- Xgrid[,10*low_corr+1:8]%*%w_prior
XTEQ_overall <- apply(dat$SIM$X[[1]][,10*low_corr+1:8],2,function(x) quantile(x,seq(0.25,0.75,by=0.05)))%*%w_prior
XTEQ_overall <- rbind(XTEQ_overall,XTEQ_overall[6]) ## final row is the median for comparisons




################################################
###         Fit Models                       ###
################################################
if(mod_version=="bkmr"){ ## bkmr (unstructured)

  fit <- bsmim2(y=y,x=dat$bkmr$X[10*low_corr+(1:8)],z=dat$covariates,niter=R,nburn=R*burn,nthin=thin,centering=F,scaling=F,prior_sigma=c(0.001,0.001),prior_lambda_shaperate=c(1,0.1),gaussian=TRUE,spike_slab=TRUE,gauss_prior=TRUE,prior_theta_slab_sd=0.25,stepsize_theta=jump,draw_h=FALSE)
  
}else if(mod_version=="SIM"){ ## SIM with no constraints
  
  fit <- bsmim2(y=y,x=list(dat$SIM$X[[1]][,10*low_corr+1:8]),z=dat$covariates,niter=R,nburn=R*burn,nthin=thin,centering=F,scaling=F,prior_sigma=c(0.001,0.001),prior_lambda_shaperate=c(1,0.1),gaussian=TRUE,spike_slab=TRUE,gauss_prior=TRUE,prior_theta_slab_sd=0.25,stepsize_theta=jump,draw_h=FALSE)

}else if(mod_version=="SIMpos"){ ## SIM with positive constraints
  
  fit <- bsmim2(y=y,x=list(dat$SIM$X[[1]][,10*low_corr+1:8]),z=dat$covariates,niter=R,nburn=R*burn,nthin=thin,centering=F,scaling=F,prior_sigma=c(0.001,0.001),prior_lambda_shaperate=c(1,0.1),gaussian=TRUE,spike_slab=TRUE,gauss_prior=TRUE,prior_theta_slab_sd=0.25,stepsize_theta=jump,draw_h=FALSE,constraints=c(1),prior_slabpos=c(1.6,8)) ## shape and rate for gamma on thetastar
 
}else if(mod_version=="SIMrank"){ ## SIM with positive constraints and ordered weights
  
  fit <- bsmim2(y=y,x=list(dat$SIM$X[[1]][,10*low_corr+1:8]),z=dat$covariates,niter=R,nburn=R*burn,nthin=thin,centering=F,scaling=F,prior_sigma=c(0.001,0.001),prior_lambda_shaperate=c(1,0.1),gaussian=TRUE,spike_slab=TRUE,gauss_prior=TRUE,prior_theta_slab_sd=0.25,stepsize_theta=jump,draw_h=FALSE,constraints=c(1),prior_slabpos=c(1.6,8),basis.opts.list=list(list(type="RANKED"))) ##checking new version c(1.6,24)

}else if(mod_version=="SIMdir"){ ## SIM with dirichlet constraints 
  
  # 10: # fit <- bsmim2(y=y,x=list(dat$SIM$X[[1]][,10*low_corr+1:8]),z=dat$covariates,niter=R,nburn=R*burn,nthin=thin,centering=F,scaling=F,prior_sigma=c(0.001,0.001),prior_lambda_shaperate=c(1,0.1),gaussian=TRUE,spike_slab=TRUE,gauss_prior=TRUE,prior_theta_slab_sd=0.25,stepsize_theta=jump,draw_h=FALSE,constraints=c(2),prior_alphas=list(10*w_prior),prior_slabrho=c(1,1))
  fit <- bsmim2(y=y,x=list(dat$SIM$X[[1]][,10*low_corr+1:8]),z=dat$covariates,niter=R,nburn=R*burn,nthin=thin,centering=F,scaling=F,prior_sigma=c(0.001,0.001),prior_lambda_shaperate=c(1,0.1),gaussian=TRUE,spike_slab=TRUE,gauss_prior=TRUE,prior_theta_slab_sd=0.25,stepsize_theta=jump,draw_h=FALSE,constraints=c(2),prior_alphas=list(5*w_prior),prior_slabrho=c(1,1))
  
}else if(mod_version=="SIMdirequal"){ ## SIM with equal dirichlet constraints 
  
  # 10: # fit <- bsmim2(y=y,x=list(dat$SIM$X[[1]][,10*low_corr+1:8]),z=dat$covariates,niter=R,nburn=R*burn,nthin=thin,centering=F,scaling=F,prior_sigma=c(0.001,0.001),prior_lambda_shaperate=c(1,0.1),gaussian=TRUE,spike_slab=TRUE,gauss_prior=TRUE,prior_theta_slab_sd=0.25,stepsize_theta=jump,draw_h=FALSE,constraints=c(2),prior_alphas=list(10*rep(1/8,8)),prior_slabrho=c(1,1))
  fit <- bsmim2(y=y,x=list(dat$SIM$X[[1]][,10*low_corr+1:8]),z=dat$covariates,niter=R,nburn=R*burn,nthin=thin,centering=F,scaling=F,prior_sigma=c(0.001,0.001),prior_lambda_shaperate=c(1,0.1),gaussian=TRUE,spike_slab=TRUE,gauss_prior=TRUE,prior_theta_slab_sd=0.25,stepsize_theta=jump,draw_h=FALSE,constraints=c(2),prior_alphas=list(5*rep(1/8,8)),prior_slabrho=c(1,1))
  
}else if(mod_version=="SIMdirvarsel"){ ## SIM with dirichlet constraints AND variable selection via informative gammas
  
  # 5: # fit <- bsmim2(y=y,x=list(dat$SIM$X[[1]][,10*low_corr+1:8]),z=dat$covariates,niter=R,nburn=R*burn,nthin=thin,centering=F,scaling=F,prior_sigma=c(0.001,0.001),prior_lambda_shaperate=c(1,0.1),gaussian=TRUE,spike_slab=TRUE,gauss_prior=TRUE,prior_theta_slab_sd=0.25,stepsize_theta=jump,draw_h=FALSE,constraints=c(1),prior_slabpos=c(1.6,8),prior_slabpos_shape_inf=list(5*w_prior)) ## using 10* since it gives roughly same sum of shape parameters as the positive constraint above
  fit <- bsmim2(y=y,x=list(dat$SIM$X[[1]][,10*low_corr+1:8]),z=dat$covariates,niter=R,nburn=R*burn,nthin=thin,centering=F,scaling=F,prior_sigma=c(0.001,0.001),prior_lambda_shaperate=c(1,0.1),gaussian=TRUE,spike_slab=TRUE,gauss_prior=TRUE,prior_theta_slab_sd=0.25,stepsize_theta=jump,draw_h=FALSE,constraints=c(1),prior_slabpos=c(1.6,8),prior_slabpos_shape_inf=list(10*w_prior)) ## using 10* since it gives roughly same sum of shape parameters as the positive constraint above
  
}else if(mod_version=="TEQ"){ ## SIM with fixed TEF weights
  
  fit <- bsmim2(y=y,x=list(XTEQ),z=dat$covariates,niter=R,nburn=R*burn,nthin=thin,centering=TRUE,scaling=TRUE,prior_sigma=c(0.001,0.001),prior_lambda_shaperate=c(1,0.1),gaussian=TRUE,spike_slab=TRUE,gauss_prior=TRUE,prior_theta_slab_sd=0.25,stepsize_theta=jump,draw_h=FALSE)
  
  
}


############################
### Summarize Results

if(mod_version=="bkmr"){ ## different for bkmr
  ## new h values
  pred_TEST <- predict_hnew_X2(fit,newX=dat$bkmr$X_TEST[10*low_corr+1:8])#,newY=y_TEST,newZ=dat$covariates_TEST)
  h_MSE <- mean((pred_TEST$fits$mean-h_TEST)^2)
  h_CVG <- mean(pred_TEST$fits$lower<=h_TEST & h_TEST<=pred_TEST$fits$upper)
  h_width <- mean(pred_TEST$fits$upper-pred_TEST$fits$lower)
  ## overall
  pred_overall <- predict_hnew_assoc2(fit,overall=TRUE,qtl_lims=c(0.25,0.75),points=11)
  overall_MSE <- mean(((pred_overall$contrasts$mean-c(overall_true))^2)[-6])
  overall_CVG <- mean((pred_overall$contrasts$lower<=c(overall_true) & c(overall_true)<=pred_overall$contrasts$upper)[-6])
  overall_width <- mean((pred_overall$contrasts$upper-pred_overall$contrasts$lower)[-6])
  ## indexwise curves
  index_MSE <- NA ## no index for bkmr
  index_CVG <- NA
  index_width <- NA
  index_grid <- c(NA,NA)
  ## componentwise curves
  pred_component <- predict_hnew2(fit,qtl_lims=c(0.25,0.75),points=11)
  component_MSE <- mean((pred_component$centered$mean-c(component_true-mean(component_true)))^2)
  component_CVG <- mean(pred_component$centered$lower<=c(component_true-mean(component_true)) & c(component_true-mean(component_true))<=pred_component$centered$upper)
  component_width <- mean(pred_component$centered$upper-pred_component$centered$lower)
  ## component weights (sum constraint)
  weight_avg <- rep(NA,8)
  weight_MSE <- rep(NA,8)
  weight_CVG <- rep(NA,8)
  weight_width <- rep(NA,8)
  ## theta weights (norm constraint)
  theta_avg <- rep(NA,8)
  theta_MSE <- rep(NA,8)
  theta_CVG <- rep(NA,8)
  theta_width <- rep(NA,8)
  # PIPs (for bkmr it is for Rho)
  PIP <- summarize_thetas(fit,w=TRUE)
  PIP <- unlist(PIP)[names(unlist(PIP))=="PIP_RHO"] ##need to unlist and then only select the PIPs for the RHOs themselves, since weights are 1
  # investigate index/feature weight rho_1
  rho_mean <- NA
  rho_PIP <- NA
  rho_sd <- NA
  
  
}else if(mod_version=="TEQ"){ ## TEQ model with fixed weights
  ## new h values
  pred_TEST <- predict_hnew_X2(fit,newX=list(XTEQ_TEST))#,newY=y_TEST,newZ=dat$covariates_TEST)
  h_MSE <- mean((pred_TEST$fits$mean-h_TEST)^2)
  h_CVG <- mean(pred_TEST$fits$lower<=h_TEST & h_TEST<=pred_TEST$fits$upper)
  h_width <- mean(pred_TEST$fits$upper-pred_TEST$fits$lower)
  ## overall
  pred_overall <- predict_hnew_assoc2(fit,newXcontrast = list(XTEQ_overall))
  overall_MSE <- mean(((pred_overall$contrasts$mean-c(overall_true))^2)[-6])
  overall_CVG <- mean((pred_overall$contrasts$lower<=c(overall_true) & c(overall_true)<=pred_overall$contrasts$upper)[-6])
  overall_width <- mean((pred_overall$contrasts$upper-pred_overall$contrasts$lower)[-6])
  ## indexwise curves
  pred_index <- predict_hnew_indexwise2(fit,qtl_lims=c(0.25,0.75),points=11) ## for averaging widths
  index_MSE <- mean((pred_index$mean_centered-c(index_true-mean(index_true)))^2) #old version; now updated for centered curves #mean((pred_index$mean-c(index_true))^2)
  index_CVG <- mean(pred_index$lower_centered<=c(index_true-mean(index_true)) & c(index_true-mean(index_true))<=pred_index$upper_centered) #old version; now updated for centered curves #mean(pred_index$lower<=c(index_true) & c(index_true)<=pred_index$upper)
  index_width <- mean(pred_index$upper_centered-pred_index$lower_centered)
  index_grid <- range(pred_index$grid)
  ## componentwise curves
  pred_component <- predict_hnew_X2(fit,newX=list(XTEQ_component))
  component_MSE <- mean((pred_component$centered$mean-c(component_true-mean(component_true)))^2)
  component_CVG <- mean(pred_component$centered$lower<=c(component_true-mean(component_true)) & c(component_true-mean(component_true))<=pred_component$centered$upper)
  component_width <- mean(pred_component$centered$upper-pred_component$centered$lower)
  ## component weights (sum constraint) ## done for each exposure individually, can sum later
    weight_avg <- w_prior       ## treated as fixed
    weight_upper <- w_prior     ## treated as fixed
    weight_lower <- w_prior     ## treated as fixed
    weight_MSE <- ((weight_avg-ww)^2)
    weight_CVG <- as.numeric(weight_lower<=ww & ww<=weight_upper)
    weight_width <- (weight_upper-weight_lower)
  ## theta weights (norm constraint) ## done for each exposure individually, can sum later
  theta_avg <- w_prior/sqrt(sum(w_prior^2))      ## treated as fixed
  theta_upper <- theta_avg                       ## treated as fixed
  theta_lower <- theta_avg                       ## treated as fixed
  theta_MSE <- ((theta_avg-wwnorm)^2)
  theta_CVG <- as.numeric(theta_lower<=wwnorm & wwnorm<=theta_upper)
  theta_width <- (theta_upper-theta_lower)
  # marginal PIPs (PIP for full index times conditional PIP for theta)
  summaryPIPs <- summarize_thetas(fit,w=TRUE)[[1]]
  PIP <- rep(summaryPIPs$PIP_RHO,8) ## all included or excluded simultaneously
  # new: investigate index/feature weight rho_1
  rho_mean <- mean(fit$rho)
  rho_PIP <- mean(fit$rho!=0)
  rho_sd <- sd(fit$rho)
  
  
}else{ ## all the SIM models
  ## new h values
  pred_TEST <- predict_hnew_X2(fit,newX=list(dat$SIM$X_TEST[[1]][,10*low_corr+1:8]))#,newY=y_TEST,newZ=dat$covariates_TEST)
  h_MSE <- mean((pred_TEST$fits$mean-h_TEST)^2)
  h_CVG <- mean(pred_TEST$fits$lower<=h_TEST & h_TEST<=pred_TEST$fits$upper)
  h_width <- mean(pred_TEST$fits$upper-pred_TEST$fits$lower)
  ## overall
  pred_overall <- predict_hnew_assoc2(fit,overall=TRUE,qtl_lims=c(0.25,0.75),points=11)
  overall_MSE <- mean(((pred_overall$contrasts$mean-c(overall_true))^2)[-6])
  overall_CVG <- mean((pred_overall$contrasts$lower<=c(overall_true) & c(overall_true)<=pred_overall$contrasts$upper)[-6])
  overall_width <- mean((pred_overall$contrasts$upper-pred_overall$contrasts$lower)[-6])
  ## indexwise curves
  pred_index <- predict_hnew_indexwise2(fit,qtl_lims=c(0.25,0.75),points=11) ## for averaging widths
  index_MSE <- mean((pred_index$mean_centered-c(index_true-mean(index_true)))^2) #old version; now updated for centered curves #mean((pred_index$mean-c(index_true))^2)
  index_CVG <- mean(pred_index$lower_centered<=c(index_true-mean(index_true)) & c(index_true-mean(index_true))<=pred_index$upper_centered) #old version; now updated for centered curves #mean(pred_index$lower<=c(index_true) & c(index_true)<=pred_index$upper)
  index_width <- mean(pred_index$upper_centered-pred_index$lower_centered)  
  index_grid <- range(pred_index$grid)
  ## componentwise curves
  pred_component <- predict_hnew2(fit,qtl_lims=c(0.25,0.75),points=11)
  component_MSE <- mean((pred_component$centered$mean-c(component_true-mean(component_true)))^2) #old version; now updated for centered curves #mean((pred_component$fits$mean-c(component_true))^2)
  component_CVG <- mean(pred_component$centered$lower<=c(component_true-mean(component_true)) & c(component_true-mean(component_true))<=pred_component$centered$upper) #old version; now updated for centered curves #mean(pred_component$fits$lower<=c(component_true) & c(component_true)<=pred_component$fits$upper)
  component_width <- mean(pred_component$centered$upper-pred_component$centered$lower)
  ## component weights (sum constraint) ## done for each exposure individually, can sum later
  if(mod_version=="SIM"){ ## if unconstrained
    weight_avg <- rep(NA,8)
    weight_MSE <- rep(NA,8)
    weight_CVG <- rep(NA,8)
    weight_width <- rep(NA,8)
  }else{
    weight_avg <- apply(fit$wPOS[[1]],2,mean)
    weight_upper <- apply(fit$wPOS[[1]],2,function(x) quantile(x,0.975))
    weight_lower <- apply(fit$wPOS[[1]],2,function(x) quantile(x,0.025))
    weight_MSE <- ((weight_avg-ww)^2)
    weight_CVG <- as.numeric(weight_lower<=ww & ww<=weight_upper)
    weight_width <- (weight_upper-weight_lower)
  }
  ## theta weights (norm constraint) ## done for each exposure individually, can sum later
  theta_avg <- apply(fit$w[[1]],2,mean)/sqrt(sum(apply(fit$w[[1]],2,mean)^2))
  theta_upper <- apply(fit$w[[1]],2,function(x) quantile(x,0.975))
  theta_lower <- apply(fit$w[[1]],2,function(x) quantile(x,0.025))
  theta_MSE <- ((theta_avg-wwnorm)^2)
  theta_CVG <- as.numeric(theta_lower<=wwnorm & wwnorm<=theta_upper)
  theta_width <- (theta_upper-theta_lower)
  # marginal PIPs (PIP for full index times conditional PIP for theta)
  summaryPIPs <- summarize_thetas(fit,w=TRUE)[[1]]
  PIP <- summaryPIPs$PIP_RHO*summaryPIPs$condPIP
  # new: investigate index/feature weight rho_1
  rho_mean <- mean(fit$rho)
  rho_PIP <- mean(fit$rho!=0)
  rho_sd <- sd(fit$rho)
}







############################
### collect results

results <- unlist(c(data.frame(iter_no=iter_no,
                               # holdout h
                               h_MSE=h_MSE,
                               h_width=h_width,
                               h_CVG=h_CVG,
                               # overall contrasts
                               overall_MSE=overall_MSE,
                               overall_width=overall_width,
                               overall_CVG=overall_CVG,
                               # indexwise curve
                               index_MSE=index_MSE,
                               index_width=index_width,
                               index_CVG=index_CVG,
                               index_grid1=index_grid[1],
                               index_grid2=index_grid[2], ## otherwise everything gets doubled
                               # componentwise curves
                               component_MSE=component_MSE,
                               component_width=component_width,
                               component_CVG=component_CVG),
                    # weights (sum constraint) (each is a vector of length 8)
                    data.frame(weight_avg),
                    data.frame(weight_MSE),
                    data.frame(weight_width),
                    data.frame(weight_CVG),
                    # thetas (norm constraint) (each is a vector of length 8)
                    data.frame(theta_avg),
                    data.frame(theta_MSE),
                    data.frame(theta_width),
                    data.frame(theta_CVG),
                    # marginal PIPs
                    data.frame(PIP),
                    # rho 
                    data.frame(rho_mean=rho_mean,
                               rho_PIP=rho_PIP, 
                               rho_sd=rho_sd) ))
results <- data.frame(t(results))
write.csv(results,file=paste0(path,mod_version,suffix,".csv"),row.names=F)








# ## investigate priors
# prior_sim <- investigate_priors(8,R=8000,constraint=0,gauss_prior=TRUE,prior_theta_slab_sd=0.25)
# prior_pos <- investigate_priors(8,R=8000,constraint=1,gauss_prior=TRUE,prior_theta_slab_sd=0.25,prior_slabpos=c(1.6,8))
# prior_direqual <- investigate_priors(8,R=8000,constraint=2,gauss_prior=TRUE,prior_theta_slab_sd=0.25,prior_slabpos=c(1.6,8),prior_alphas=(10*rep(1/8,8)),prior_slabrho=c(1,1))
# prior_dir <- investigate_priors(8,R=8000,constraint=2,gauss_prior=TRUE,prior_theta_slab_sd=0.25,prior_slabpos=c(1.6,8),prior_alphas=(10*w_prior),prior_slabrho=c(1,1))
# prior_dirvarsel <- investigate_priors(8,R=8000,constraint=1,gauss_prior=TRUE,prior_theta_slab_sd=0.25,prior_slabpos=c(1.6,8),prior_slabpos_shape_inf=5*w_prior)
# 
# boxplot(prior_sim$theta,ylim=c(-1,1))
# boxplot(prior_constrained$theta,ylim=c(-1,1))
# boxplot(prior_direqual$theta,ylim=c(-1,1))
# boxplot(prior_dir$theta,ylim=c(-1,1))
# boxplot(prior_dirvarsel$theta,ylim=c(-1,1))
# 
# boxplot(prior_constrained$thetaPOS)
# boxplot(prior_direqual$thetaPOS)
# boxplot(prior_dir$thetaPOS)
# boxplot(prior_dirvarsel$thetaPOS)
