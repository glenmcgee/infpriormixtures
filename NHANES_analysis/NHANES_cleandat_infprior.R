########################################################
###       Clean Data for NHANES Sims+Analysis        ###
########################################################
## October 11, 2021
## Script to clean NHANES data 
## Also defines analysis functions


# doLog=FALSE in main script that calls this one
# swapX=FALSE in main script that calls this one
##########################
###  Pre-process data  ###  
##########################
## based on code from Katrina Devick for CU mixtures workshop

## center/scale continous covariates and create indicators for categorical covariates
nhanes$age_z         <- scale(nhanes$age_cent)         ## center and scale age
nhanes$agez_sq       <- nhanes$age_z^2                 ## square this age variable
nhanes$bmicat2       <- as.numeric(nhanes$bmi_cat3==2) ## 25 <= BMI < 30
nhanes$bmicat3       <- as.numeric(nhanes$bmi_cat3==3) ## BMI >= 30 (BMI < 25 is the reference)
nhanes$educat1       <- as.numeric(nhanes$edu_cat==1)  ## no high school diploma
nhanes$educat3       <- as.numeric(nhanes$edu_cat==3)  ## some college or AA degree
nhanes$educat4       <- as.numeric(nhanes$edu_cat==4)  ## college grad or above (reference is high schol grad/GED or equivalent)
nhanes$otherhispanic <- as.numeric(nhanes$race_cat==1) ## other Hispanic or other race - including multi-racial
nhanes$mexamerican   <- as.numeric(nhanes$race_cat==2) ## Mexican American 
nhanes$black         <- as.numeric(nhanes$race_cat==3) ## non-Hispanic Black (non-Hispanic White as reference group)
nhanes$wbcc_z        <- scale(nhanes$LBXWBCSI)
nhanes$lymphocytes_z <- scale(nhanes$LBXLYPCT)
nhanes$monocytes_z   <- scale(nhanes$LBXMOPCT)
nhanes$neutrophils_z <- scale(nhanes$LBXNEPCT)
nhanes$eosinophils_z <- scale(nhanes$LBXEOPCT)
nhanes$basophils_z   <- scale(nhanes$LBXBAPCT)
nhanes$lncotinine_z  <- scale(nhanes$ln_lbxcot)         ## to access smoking status, scaled ln cotinine levels


## exposure matrix
mixture <- with(nhanes, cbind(LBX074LA, LBX099LA, LBX118LA, LBX138LA, LBX153LA, LBX170LA, LBX180LA, LBX187LA,LBX194LA, LBXHXCLA, LBXPCBLA,
                              LBXD03LA, LBXD05LA, LBXD07LA,
                              LBXF03LA, LBXF04LA, LBXF05LA, LBXF08LA)) 
if(doLog==TRUE){
  lnmixture   <- apply(mixture, 2, log)
}else{
  lnmixture   <- mixture
}

#######
## save standard deviations for TEQs
xsd <- apply(lnmixture,2,sd,na.rm=T)
xsd1 <- xsd[c(1,2,4,5,6,7,8,9)]
xsd2 <- xsd[c(10,11)]
xsd3 <- xsd[c(3,12:18)]

#######
## standardize exposures
X <- scale(lnmixture)
colnames(X) <- c(paste0("PCB",c("074", "099", 118, 138, 153, 170, 180, 187, 194, 169, 126)),paste0("Dioxin",1:3), paste0("Furan",1:4)) 

## clean up exposure names and reorder
colnames(X)[c(1,2,4,5,6,7,8,9)] <- paste0("A-",colnames(X)[c(1,2,4,5,6,7,8,9)])
colnames(X)[c(10,11)] <- paste0("B-",colnames(X)[c(10,11)])
colnames(X)[c(3)] <- paste0("C--",colnames(X)[c(3)])
colnames(X)[c(12:18)] <- paste0("C-",colnames(X)[c(12:18)])
X <- X[,sort(colnames(X))]
colnames(X) <- substring(colnames(X),3); colnames(X)[substring(colnames(X),1,1)=="-"] <- substring(colnames(X)[substring(colnames(X),1,1)=="-"],2)

## if swapX==TRUE, permute exposures within groups to avoid artifacts in simulations due to unique correlation structure
if(exists("swapX")){
  if(swapX==TRUE){
    X[,1:8] <- X[,sample(1:8)]
    X[,9:10] <- X[,sample(9:10)]
    X[,11:18] <- X[,sample(11:18)]
  }
}

## covariates
covariates <- with(nhanes, cbind(age_z, agez_sq, male, bmicat2, bmicat3))

## outcome
lnLTL_z <- scale(log(nhanes$TELOMEAN))

################################################
###        Prepare data for predictions      ###
################################################

## Construct new grid points, quantiles and weights
getGrid <- function(qtl=0.5,pts=20,qtl_lims=c(0.05,0.95)){
  Xq <- c()
  for(ll in 1:ncol(X)){
    tempXq <- matrix(apply(X,2,function(x) quantile(x,qtl)),nrow=pts,ncol=ncol(X),byrow=T)
    tempXq[,ll] <- seq(quantile(X[,ll],qtl_lims[1]),quantile(X[,ll],qtl_lims[2]),length=pts)
    Xq <- rbind(Xq,tempXq)
  }
  return(Xq)
}
X25 <- getGrid(qtl=0.25)
X50 <- getGrid(qtl=0.50)
X75 <- getGrid(qtl=0.75)



################################################
###     prepare data for model fitting       ###
################################################

subset_data <- function(ids){
  
  ## exposure 
  X <- as.matrix(X[ids,])
  
  ##  exposures for model fits
  X_bkmr <- list(as.matrix(X[,1]),as.matrix(X[,2]),as.matrix(X[,3]),as.matrix(X[,4]),as.matrix(X[,5]),as.matrix(X[,6]),as.matrix(X[,7]),as.matrix(X[,8]),as.matrix(X[,9]),as.matrix(X[,10]),as.matrix(X[,11]),as.matrix(X[,12]),as.matrix(X[,13]),as.matrix(X[,14]),as.matrix(X[,15]),as.matrix(X[,16]),as.matrix(X[,17]),as.matrix(X[,18]))
  X_SIM <- list(X[,1:18])
  X_bsmim <- list(X[,1:8],X[,9:10],X[,11:18])
  
  ## our covariate matrix
  covariates <- covariates[ids,]
  
  Y <- lnLTL_z[ids]
  
  ## return data
  return(list(X=X,X_bkmr=X_bkmr,X_SIM=X_SIM,X_bsmim=X_bsmim,covariates=covariates,Y=Y))
}


prep_data_full <- function(){
  
  traindat <- subset_data(1:nrow(nhanes))

  ### exposures for predictions
  X25_bkmr <- list(as.matrix(X25[,1]),as.matrix(X25[,2]),as.matrix(X25[,3]),as.matrix(X25[,4]),as.matrix(X25[,5]),as.matrix(X25[,6]),as.matrix(X25[,7]),as.matrix(X25[,8]),as.matrix(X25[,9]),as.matrix(X25[,10]),as.matrix(X25[,11]),as.matrix(X25[,12]),as.matrix(X25[,13]),as.matrix(X25[,14]),as.matrix(X25[,15]),as.matrix(X25[,16]),as.matrix(X25[,17]),as.matrix(X25[,18]))
  X50_bkmr <- list(as.matrix(X50[,1]),as.matrix(X50[,2]),as.matrix(X50[,3]),as.matrix(X50[,4]),as.matrix(X50[,5]),as.matrix(X50[,6]),as.matrix(X50[,7]),as.matrix(X50[,8]),as.matrix(X50[,9]),as.matrix(X50[,10]),as.matrix(X50[,11]),as.matrix(X50[,12]),as.matrix(X50[,13]),as.matrix(X50[,14]),as.matrix(X50[,15]),as.matrix(X50[,16]),as.matrix(X50[,17]),as.matrix(X50[,18]))
  X75_bkmr <- list(as.matrix(X75[,1]),as.matrix(X75[,2]),as.matrix(X75[,3]),as.matrix(X75[,4]),as.matrix(X75[,5]),as.matrix(X75[,6]),as.matrix(X75[,7]),as.matrix(X75[,8]),as.matrix(X75[,9]),as.matrix(X75[,10]),as.matrix(X75[,11]),as.matrix(X75[,12]),as.matrix(X75[,13]),as.matrix(X75[,14]),as.matrix(X75[,15]),as.matrix(X75[,16]),as.matrix(X75[,17]),as.matrix(X75[,18]))
  X25_SIM <- list(X25[,1:18])
  X50_SIM <- list(X50[,1:18])
  X75_SIM <- list(X75[,1:18])
  X25_bsmim <- list(X25[,1:8],X25[,9:10],X25[,11:18])
  X50_bsmim <- list(X50[,1:8],X50[,9:10],X50[,11:18])
  X75_bsmim <- list(X75[,1:8],X75[,9:10],X75[,11:18])
  
  ## return data
  SIM <- list(X=traindat$X_SIM,     X_TEST=NULL, X25=X25_SIM,   X50=X50_SIM,   X75=X75_SIM)
  bsmim <- list(X=traindat$X_bsmim, X_TEST=NULL, X25=X25_bsmim, X50=X50_bsmim, X75=X75_bsmim)
  bkmr <- list(X=traindat$X_bkmr,   X_TEST=NULL, X25=X25_bkmr,  X50=X50_bkmr,  X75=X75_bkmr)
  return(list(SIM=SIM,bsmim=bsmim,bkmr=bkmr,covariates=traindat$covariates,covariates_TEST=NULL,y=traindat$Y,y_TEST=NULL))
  
}

prep_data_split <- function(train_ids){
  
  ## split data
  trainset <- train_ids
  testset <- (1:length(lnLTL_z))[-train_ids]
  
  traindat <- subset_data(trainset)
  testdat <- subset_data(testset)
  
  ### exposures for predictions
  X25_bkmr <- list(as.matrix(X25[,1]),as.matrix(X25[,2]),as.matrix(X25[,3]),as.matrix(X25[,4]),as.matrix(X25[,5]),as.matrix(X25[,6]),as.matrix(X25[,7]),as.matrix(X25[,8]),as.matrix(X25[,9]),as.matrix(X25[,10]),as.matrix(X25[,11]),as.matrix(X25[,12]),as.matrix(X25[,13]),as.matrix(X25[,14]),as.matrix(X25[,15]),as.matrix(X25[,16]),as.matrix(X25[,17]),as.matrix(X25[,18]))
  X50_bkmr <- list(as.matrix(X50[,1]),as.matrix(X50[,2]),as.matrix(X50[,3]),as.matrix(X50[,4]),as.matrix(X50[,5]),as.matrix(X50[,6]),as.matrix(X50[,7]),as.matrix(X50[,8]),as.matrix(X50[,9]),as.matrix(X50[,10]),as.matrix(X50[,11]),as.matrix(X50[,12]),as.matrix(X50[,13]),as.matrix(X50[,14]),as.matrix(X50[,15]),as.matrix(X50[,16]),as.matrix(X50[,17]),as.matrix(X50[,18]))
  X75_bkmr <- list(as.matrix(X75[,1]),as.matrix(X75[,2]),as.matrix(X75[,3]),as.matrix(X75[,4]),as.matrix(X75[,5]),as.matrix(X75[,6]),as.matrix(X75[,7]),as.matrix(X75[,8]),as.matrix(X75[,9]),as.matrix(X75[,10]),as.matrix(X75[,11]),as.matrix(X75[,12]),as.matrix(X75[,13]),as.matrix(X75[,14]),as.matrix(X75[,15]),as.matrix(X75[,16]),as.matrix(X75[,17]),as.matrix(X75[,18]))
  X25_SIM <- list(X25[,1:18])
  X50_SIM <- list(X50[,1:18])
  X75_SIM <- list(X75[,1:18])
  X25_bsmim <- list(X25[,1:8],X25[,9:10],X25[,11:18])
  X50_bsmim <- list(X50[,1:8],X50[,9:10],X50[,11:18])
  X75_bsmim <- list(X75[,1:8],X75[,9:10],X75[,11:18])
  
  ## return data
  SIM <- list(X=traindat$X_SIM,     X_TEST=testdat$X_SIM,   X25=X25_SIM,   X50=X50_SIM,   X75=X75_SIM)
  bsmim <- list(X=traindat$X_bsmim, X_TEST=testdat$X_bsmim, X25=X25_bsmim, X50=X50_bsmim, X75=X75_bsmim)
  bkmr <- list(X=traindat$X_bkmr,   X_TEST=testdat$X_bkmr, X25=X25_bkmr,  X50=X50_bkmr,  X75=X75_bkmr)
  return(list(SIM=SIM,bsmim=bsmim,bkmr=bkmr,covariates=traindat$covariates,covariates_TEST=testdat$covariates,y=traindat$Y,y_TEST=testdat$Y))
  
}



################################################
###        functions for model fitting       ###
################################################


## compute overall effect from qgcomp fit
getOverall_qgcomp <- function(obj,qq=10,q2=6,q1=5,degree=2){ 
  
  quants <- (1:qq)-1
  qmat <- c()
  for(dd in 1:degree){
    qmat <- cbind(qmat,quants^dd)
  }
  newmat <- cbind(1,matrix(qmat,nrow=length(quants),ncol=degree*length(obj$expnms)),matrix(0,nrow=length(quants),ncol=ncol(covariates)))
  
  fitted <- newmat%*%obj$fit$coefficients
  Vmat <- newmat%*%vcov(obj$fit)%*%t(newmat)
  
  contrastvec <- rep(0,length(fitted))
  contrastvec[q2] <- 1    ## (q2,q1)=(6,5)--> compare 60th to 50th percentile
  contrastvec[q1] <- -1
  est_diff <- contrastvec%*%fitted
  SE_diff <- sqrt(contrastvec%*%Vmat%*%contrastvec)
  lci <- est_diff-1.96*SE_diff
  uci <- est_diff+1.96*SE_diff
  res <- c(est_diff,SE_diff,lci,uci)
  
  names(res) <- c("mean","sd","lower","upper")
  return(res)
}

## get quantiles
getQuants <- function(obj,newX,newCov){
  
  Xtemp <- newX #
  for(jj in 1:length(obj$expnms)){
    Xtemp[,jj] <- cut(Xtemp[,jj], breaks = obj$breaks[[jj]], labels = FALSE,include.lowest = TRUE) - 1
  }  
  colnames(Xtemp) <- obj$expnms
  res <- data.frame(cbind(1,Xtemp,newCov))        #
  
  return(res)
}

pred_twoway <- function(obj,qtls=c(0.1,0.9),qtl_lims=c(0.05,0.95),pts=20,include_median=TRUE){
  
  ## get predictions at each level (skip 0.5 since it is implicitly computed)
  res <- list()
  for(mm in 1:ncol(obj$rho)){
    res_mm <- list()
    for(qq in 1:length(qtls)){
      res_mm[[qq]] <- predict_hnew_indexwise2(obj,crossM=mm,qtl=qtls[[qq]],qtl_lims=qtl_lims,points=pts)
    }
    names(res_mm) <- qtls  ## label
    res[[mm]] <- res_mm
  }
  
  ## combine predictions into dataframe for plotting
  df_var1 <- df_var2 <- df_grid <- df_quantile <- df_est <- df_est_centered <- df_sd <- df_sd_centered <- c()
  for(xx in 1:ncol(obj$rho)){
    for(yy in 1:ncol(obj$rho)){
      if(xx==yy){
        next
      }
      for(qq in 1:length(qtls)){
        df_var1 <- c(df_var1,rep(paste0("Index ",xx),pts))
        df_var2 <- c(df_var2,rep(paste0("Index ",yy),pts))
        df_grid <- c(df_grid,res[[yy]][[qq]]$grid[(xx-1)*pts+1:pts])
        df_est <- c(df_est,res[[yy]][[qq]]$mean[(xx-1)*pts+1:pts])
        df_sd <- c(df_sd,res[[yy]][[qq]]$sd[(xx-1)*pts+1:pts])
        df_est_centered <- c(df_est_centered,res[[yy]][[qq]]$mean_centered[(xx-1)*pts+1:pts])
        df_sd_centered <- c(df_sd_centered,res[[yy]][[qq]]$sd_centered[(xx-1)*pts+1:pts])
        df_quantile <- c(df_quantile,rep(qtls[qq],pts))
      }
      if(include_median==TRUE){ ## implicitly predicted above
        df_var1 <- c(df_var1,rep(paste0("Index ",xx),pts))
        df_var2 <- c(df_var2,rep(paste0("Index ",yy),pts))
        df_grid <- c(df_grid,res[[xx]][[qq]]$grid[(xx-1)*pts+1:pts]) ## set yy (i.e. crossM) to xx (the index being predicted), which means we implicitly set everything else at the median
        df_est <- c(df_est,res[[xx]][[qq]]$mean[(xx-1)*pts+1:pts])   ## set yy (i.e. crossM) to xx (the index being predicted), which means we implicitly set everything else at the median
        df_sd <- c(df_sd,res[[xx]][[qq]]$sd[(xx-1)*pts+1:pts])   ## set yy (i.e. crossM) to xx (the index being predicted), which means we implicitly set everything else at the median
        df_est_centered <- c(df_est_centered,res[[xx]][[qq]]$mean_centered[(xx-1)*pts+1:pts])   ## set yy (i.e. crossM) to xx (the index being predicted), which means we implicitly set everything else at the median
        df_sd_centered <- c(df_sd_centered,res[[xx]][[qq]]$sd_centered[(xx-1)*pts+1:pts])   ## set yy (i.e. crossM) to xx (the index being predicted), which means we implicitly set everything else at the median
        df_quantile <- c(df_quantile,rep(0.5,pts))
      }
      
    }
  }
  
  pred_df <- data.frame(var1=df_var1,var2=df_var2,grid=df_grid,quantile=df_quantile,
                        est=df_est,est_centered=df_est_centered,
                        sd=df_sd,sd_centered=df_sd_centered)
  
  return(pred_df)
}

## get 2-way interactions for kmbayes (bkmr packages)
predict_km_inter<- function(obj){
  pred.resp.bivar <- PredictorResponseBivar(fit=obj,sel=sel,center=FALSE)
  pred.resp.bivar.levels <- PredictorResponseBivarLevels(pred.resp.df = pred.resp.bivar, qs = c(0.1, 0.5, 0.9), Z = obj$Z)
  return(pred.resp.bivar.levels)
}


## k-fold cross-validation for qgcomp
gqcomp_crossval <- function(obj,kfolds=4){
  
  id_k <- sample(rep(1:kfolds,ceiling(length(obj$fit$data$y)/kfolds))[1:length(obj$fit$data$y)],replace=F)
  df <- fit$fit$data
  
  test_out <- pred_out <- c()
  for(vv in 1:kfolds){
    
    fit_cv <- lm(formula=paste0("y~",paste(obj$fit$formula)[3]),data=df[id_k!=vv,])
    pred_out <- c(pred_out,predict(fit_cv,newdata=df[id_k==vv,]))
    test_out <- c(test_out,df$y[id_k==vv])
    
  }
  mse_out <- mean((pred_out-test_out)^2)
  bias_out <- mean((pred_out-test_out))
  
  return(list(mse_out=mse_out,bias_out=bias_out))
}

## estimate h and h+Zgamma (fitted vals) under QGC
pred_QGC <- function(fit,Xnew,Znew,fitted=FALSE,exp_names){
  ## mean of h component
  df_TEST <- getQuants(fit,Xnew,0*Znew)  
  df_NEW <- model.matrix(delete.response(terms(fit$fit)), df_TEST)
  h <- data.frame(mean=as.matrix(df_NEW)%*%fit$fit$coefficients, #predict(fit$fit,newdata=df_TEST), ## same as: as.matrix(df_TEST)%*%fit$fit$coefficients
                  sd=sqrt(diag(as.matrix(df_NEW)%*%vcov(fit$fit)%*%t(as.matrix(df_NEW)))) )

  ## mean fitted vals
  mean_fitted=NULL
  if(fitted==TRUE){
    df_TEST <- getQuants(fit,Xnew,Znew)   
    # mean_fitted <- predict(fit$fit,expnms=exp_names,newdata=df_TEST)
    df_NEW <- model.matrix(delete.response(terms(fit$fit)), df_TEST)
    mean_fitted <- data.frame(mean=as.matrix(df_NEW)%*%fit$fit$coefficients, #predict(fit$fit,newdata=df_TEST), ## same as: as.matrix(df_TEST)%*%fit$fit$coefficients
                    sd=sqrt(diag(as.matrix(df_NEW)%*%vcov(fit$fit)%*%t(as.matrix(df_NEW)))) )
    
  }

  return(list(fits=h,mean_fitted=mean_fitted))
}





## k-fold cross-validation for kmbayes
kmbayes_crossval <- function(obj,sel,groups_vec=NULL,kfolds=4){
  
  id_k <- sample(rep(1:kfolds,ceiling(length(obj$y)/kfolds))[1:length(obj$y)],replace=F)
  
  test_out <- pred_out <- c()
  for(vv in 1:kfolds){
    
    ## fit model
    if(!is.null(groups_vec)){
      fit_cv <-  kmbayes(y=obj$y[id_k!=vv], Z=obj$Z[id_k!=vv,], X=obj$X[id_k!=vv,], iter=R, verbose=FALSE, varsel=TRUE,groups=groups_vec)
    }else{
      fit_cv <-  kmbayes(y=obj$y[id_k!=vv], Z=obj$Z[id_k!=vv,], X=obj$X[id_k!=vv,], iter=R, verbose=FALSE, varsel=TRUE)
    }
    
    ## predict
    pred_TEST <- SamplePred(fit_cv,Znew=obj$Z[id_k==vv,],Xnew=obj$X[id_k==vv,],sel=sel)
    pred_out <- c(pred_out,apply(pred_TEST,2,mean))
    test_out <- c(test_out,obj$y[id_k==vv])
    
  }
  mse_out <- mean((pred_out-test_out)^2)
  bias_out <- mean((pred_out-test_out))
  
  return(list(mse_out=mse_out,bias_out=bias_out))
}



##function to get CIs from posterior draws 
pred_CI<- function(pred){
  pred_CI <- cbind(apply(pred,2,mean),apply(pred,2,sd)) ## get mean and sd
  pred_CI <- cbind(pred_CI,pred_CI[,1]-1.96*pred_CI[,2],pred_CI[,1]+1.96*pred_CI[,2]) ##
  colnames(pred_CI) <- c("mean","sd","lower","upper")
  return(pred_CI)
}

## get CIs from ComputePostmeanHnew() results 
get_CI<- function(hest){
  hest_CI <- data.frame(mean=hest$postmean,sd=sqrt(diag(hest$postvar)))
  hest_CI$lower <- hest_CI$mean-1.96*hest_CI$sd
  hest_CI$upper <- hest_CI$mean+1.96*hest_CI$sd
  return(hest_CI)
}

## get CIs from OverallRiskSummaries
ovrl_CI<- function(ovrl){
  ovrl_CI <- ovrl
  names(ovrl_CI)[2] <- "mean"
  ovrl_CI$lower <- ovrl_CI$mean-1.96*ovrl_CI$sd
  ovrl_CI$upper <- ovrl_CI$mean+1.96*ovrl_CI$sd  
  return(ovrl_CI)
}

## function to get posterior mean of h, post mean of h+Zgamma, and posterior summaries of draws of h+Zgamma
get_predkm <- function(fit,Znew,Xnew,sel,fitted=FALSE,samples=FALSE){
  ## get posterior means and sds of h
  pred_TEST_fits <- ComputePostmeanHnew(fit, Znew = Znew, sel = sel)#, method = "exact" ## using approx method as exact method errors out due to memory
  pred_TEST_fits <- get_CI(pred_TEST_fits)
  ## get posterior mean for h+Zgamma
  pred_TEST_fitted=NULL
  if(fitted==TRUE){
    pred_TEST_fitted <- pred_TEST_fits$mean+Xnew%*%apply(fit$beta[sel,],2,mean) 
  }
  ## summarize posterior for draws of h+Zgamma
  pred_TEST_samples=NULL
  if(samples==TRUE){
    pred_TEST_samples <- SamplePred(fit,Znew=Znew,Xnew=Xnew,sel=sel)
    pred_TEST_samples <- pred_CI(pred_TEST_samples) ## make CIs
  }
  ## collect results
  pred_TEST <-  list(fits=pred_TEST_fits,mean_fitted=pred_TEST_fitted,pred_out=pred_TEST_samples)
  return(pred_TEST)
}


