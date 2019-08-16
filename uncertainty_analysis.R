#### Packages ####
#install.packages("pacman")
library(pacman)
pacman::p_load("ggplot2","reshape2","drtmle","grf","foreach","uplift","data.table","ModelMetrics", "plyr", "parallel", "doParallel", "SuperLearner")

source("t_logit.R")
source("Qini.R")
source("data_generating_process.R")
source("supervised_randomization.R")
source("costs.R")
calc_ATE <- function(y, w, prop_score){
  #if(is.null(prop_score)){
  #  prop_score = rep(0.5, length(y))
  #}
  return( (sum(y*w/prop_score) - sum(y*(1-w)/(1-prop_score)) ) /length(y) )
}

#### Simulation Classification ####
Y = NULL
W = NULL

N_VAR=10
N_CUSTOMER=20000

set.seed(1234)
N_ITER = 100

temp <- foreach(iter=1:N_ITER, .multicombine = FALSE, .combine='+', .final=function(x)x/N_ITER)%do%{
  expCtrl <- expControl(n_var = N_VAR, mode = "regression", 
                        beta_sd = 1 , beta_zero = 0,  # >0 indicates more than 50% purchasers (for linear:  for 5%, for nonlinear: -3.12)
                        beta_tau_sd = 0.2, tau_zero = 0, # >0 indicates positive treatment effect (for linear: 0.645 for 5%, for nonlinear: 0.7)
                        DGP="nonlinear")
  
  X  <- data.frame(sapply(1:N_VAR, function(x)  rnorm(N_CUSTOMER, 0, 10)))  #runif(N_CUSTOMER, min = -10, max=10)
  
  
  balanced   <- do_experiment(X, expControl = expCtrl, prop_score = 0.5)
  
  prop_score <- map_propensity_logistic(balanced$tau, min_score = quantile(balanced$tau, p=0.05), max_score = quantile(balanced$tau, p=0.95)) 
  supervised <- do_experiment(X, expControl = expCtrl, prop_score = prop_score)
  
  
  targeting_model_balanced   <-           T_Regression(X=balanced$X,   y=balanced$y,   W=balanced$w)#, prop_score = rep(0.5, times=nrow(balanced$X)))
  targeting_model_supervised <-           T_Regression(X=supervised$X, y=supervised$y, W=supervised$w)#, prop_score = supervised$prop_score)
  targeting_model_supervised_corrected <- T_Regression(X=supervised$X, y=supervised$y, W=supervised$w, prop_score = supervised$prop_score)
  
  #targeting_model_balanced   <- T_rpart(X=balanced$X,   y=balanced$y,   W=balanced$w, method='anova')#
  #targeting_model_supervised <- T_rpart(X=supervised$X, y=supervised$y, W=supervised$w, method='anova')
  #targeting_model_supervised_corrected <- T_rpart(X=supervised$X, y=supervised$y, W=supervised$w, prop_score = supervised$prop_score, method='anova')
  
  tau_balanced   <- predict(targeting_model_balanced,   X)
  tau_supervised <- predict(targeting_model_supervised, X)
  tau_supervised_corrected <- predict(targeting_model_supervised_corrected, X)
  
  pred0_balanced   <- predict(targeting_model_balanced$model0,   X)
  pred0_supervised <- predict(targeting_model_supervised$model0, X)
  pred0_supervised_corrected <- predict(targeting_model_supervised_corrected$model0, X)
  
  pred1_balanced   <- predict(targeting_model_balanced$model1,   X)
  pred1_supervised <- predict(targeting_model_supervised$model1, X)
  pred1_supervised_corrected <- predict(targeting_model_supervised_corrected$model1, X)
  
  mapping <- factor(cut(prop_score, breaks=seq(0.05,0.95,length.out = 10)))
  #mapping <- factor(cut(prop_score, include_lowest=TRUE, breaks=quantile(prop_score, probs = seq(0,1,0.1)))) 
  
  rbind(
    "AB balanced (tau)"              = round(tapply((tau_balanced  - balanced$tau)^2,   mapping, mean), 3),
    "supervised randomization (tau)" = round(tapply((tau_supervised- supervised$tau)^2, mapping, mean), 3),
    "supervised randomization corrected (tau)" = round(tapply((tau_supervised_corrected- supervised$tau)^2, mapping, mean), 3),
    
    # "AB balanced"              = round(tapply( ((pred0_balanced   - balanced$y)^2)[balanced$w==0], mapping[balanced$w==0], mean), 3),
    # #"AB imbalanced"            = round(tapply( ((pred0_imbalanced - balanced$y)^2)[balanced$w==0], mapping[balanced$w==0], mean), 3),
    # "supervised randomization" = round(tapply( ((pred0_supervised - balanced$y)^2)[balanced$w==0], mapping[balanced$w==0], mean), 3),
    # 
    # "AB balanced"              = round(tapply( ((pred1_balanced   - balanced$y )^2)[balanced$w==1], mapping[balanced$w==1], mean), 3),
    # #"AB imbalanced"            = round(tapply( ((pred1_imbalanced - balanced$y )^2)[balanced$w==1], mapping[balanced$w==1], mean), 3),
    # "supervised randomization" = round(tapply( ((pred1_supervised - balanced$y )^2)[balanced$w==1], mapping[balanced$w==1], mean), 3)
    
    # Should be good for lower propensity scores and worse for higher prop. scores
    "AB balanced (model 0)"              =          round(tapply( ((pred0_balanced             - balanced$y)^2)[balanced$w==0], mapping[balanced$w==0], mean), 3),
    "supervised randomization (model 0)" =          round(tapply( ((pred0_supervised           - balanced$y)^2)[balanced$w==0], mapping[balanced$w==0], mean), 3),
    "supervised randomization corrected (model 0)"= round(tapply( ((pred0_supervised_corrected - balanced$y)^2)[balanced$w==0], mapping[balanced$w==0], mean), 3),
    
    # Should be worse for lower propensity scores and good/better for higher prop. scores
    "AB balanced (model 1)"              =          round(tapply( ((pred1_balanced             - balanced$y )^2)[balanced$w==1], mapping[balanced$w==1], mean), 3),
    "supervised randomization (model 1)" =          round(tapply( ((pred1_supervised           - balanced$y )^2)[balanced$w==1], mapping[balanced$w==1], mean), 3),
    "supervised randomization corrected (model 1)"= round(tapply( ((pred1_supervised_corrected - balanced$y )^2)[balanced$w==1], mapping[balanced$w==1], mean), 3)
  )
}
print(round(temp,2))



# ------------- Classification -----------
temp <- foreach(iter=1:N_ITER, .multicombine = FALSE, .combine='+', .final=function(x)x/N_ITER)%do%{
  expCtrl <- expControl(n_var = N_VAR, mode = "classification", 
                        beta_sd = 1, beta_zero = 0,  # >0 indicates more than 50% purchasers (for linear:  for 5%, for nonlinear: -3.12)
                        beta_tau_sd = 0.5, tau_zero = 0, # >0 indicates positive treatment effect (for linear: 0.645 for 5%, for nonlinear: 0.7)
                        DGP="nonlinear")
  
  X  <- data.frame(sapply(1:N_VAR, function(x)  rnorm(N_CUSTOMER, 0, 10)))  #runif(N_CUSTOMER, min = -10, max=10)
  
  
  balanced   <- do_experiment(X, expControl = expCtrl, prop_score = 0.5)
  
  prop_score <- map_propensity_logistic(balanced$tau, min_score = quantile(balanced$tau, p=0.05), max_score = quantile(balanced$tau, p=0.95)) 
  supervised <- do_experiment(X, expControl = expCtrl, prop_score = prop_score)
  
  
  targeting_model_balanced   <-           T_Logit(X=balanced$X,   y=balanced$y,   W=balanced$w)#, prop_score = rep(0.5, times=nrow(balanced$X)))
  targeting_model_supervised <-           T_Logit(X=supervised$X, y=supervised$y, W=supervised$w)#, prop_score = supervised$prop_score)
  targeting_model_supervised_corrected <- T_Logit(X=supervised$X, y=supervised$y, W=supervised$w, prop_score = supervised$prop_score)
  
  #targeting_model_balanced   <- T_rpart(X=balanced$X,   y=balanced$y,   W=balanced$w, method='anova')#
  #targeting_model_supervised <- T_rpart(X=supervised$X, y=supervised$y, W=supervised$w, prop_score = supervised$prop_score, method='anova')
  
  tau_balanced   <- predict(targeting_model_balanced,   X, type="response")
  tau_supervised <- predict(targeting_model_supervised, X, type="response")
  tau_supervised_corrected <- predict(targeting_model_supervised_corrected, X, type="response")
  
  pred0_balanced   <- predict(targeting_model_balanced$model0,   X, type="response")
  pred0_supervised <- predict(targeting_model_supervised$model0, X, type="response")
  pred0_supervised_corrected <- predict(targeting_model_supervised_corrected$model0, X, type="response")
  
  pred1_balanced   <- predict(targeting_model_balanced$model1,   X, type="response")
  pred1_supervised <- predict(targeting_model_supervised$model1, X, type="response")
  pred1_supervised_corrected <- predict(targeting_model_supervised_corrected$model1, X, type="response")
  
  mapping <- factor(cut(prop_score, breaks=seq(0.05,0.95,length.out = 10)))
  #mapping <- factor(cut(prop_score, include_lowest=TRUE, breaks=quantile(prop_score, probs = seq(0,1,0.1)))) 
  
  rbind(
    "AB balanced (tau)"              = round(tapply((tau_balanced  - balanced$tau)^2,   mapping, mean), 3),
    "supervised randomization (tau)" = round(tapply((tau_supervised- supervised$tau)^2, mapping, mean), 3),
    "supervised randomization corrected (tau)" = round(tapply((tau_supervised_corrected- supervised$tau)^2, mapping, mean), 3),
    
    # "AB balanced"              = round(tapply( ((pred0_balanced   - balanced$y)^2)[balanced$w==0], mapping[balanced$w==0], mean), 3),
    # #"AB imbalanced"            = round(tapply( ((pred0_imbalanced - balanced$y)^2)[balanced$w==0], mapping[balanced$w==0], mean), 3),
    # "supervised randomization" = round(tapply( ((pred0_supervised - balanced$y)^2)[balanced$w==0], mapping[balanced$w==0], mean), 3),
    # 
    # "AB balanced"              = round(tapply( ((pred1_balanced   - balanced$y )^2)[balanced$w==1], mapping[balanced$w==1], mean), 3),
    # #"AB imbalanced"            = round(tapply( ((pred1_imbalanced - balanced$y )^2)[balanced$w==1], mapping[balanced$w==1], mean), 3),
    # "supervised randomization" = round(tapply( ((pred1_supervised - balanced$y )^2)[balanced$w==1], mapping[balanced$w==1], mean), 3)
    
    # Should be good for lower propensity scores and worse for higher prop. scores
    "AB balanced (model 0)"              =          round(tapply( ((pred0_balanced             - balanced$y)^2)[balanced$w==0], mapping[balanced$w==0], mean), 3),
    "supervised randomization (model 0)" =          round(tapply( ((pred0_supervised           - balanced$y)^2)[balanced$w==0], mapping[balanced$w==0], mean), 3),
    "supervised randomization corrected (model 0)"= round(tapply( ((pred0_supervised_corrected - balanced$y)^2)[balanced$w==0], mapping[balanced$w==0], mean), 3),
    
    # Should be worse for lower propensity scores and good/better for higher prop. scores
    "AB balanced (model 1)"              =          round(tapply( ((pred1_balanced             - balanced$y )^2)[balanced$w==1], mapping[balanced$w==1], mean), 3),
    "supervised randomization (model 1)" =          round(tapply( ((pred1_supervised           - balanced$y )^2)[balanced$w==1], mapping[balanced$w==1], mean), 3),
    "supervised randomization corrected (model 1)"= round(tapply( ((pred1_supervised_corrected - balanced$y )^2)[balanced$w==1], mapping[balanced$w==1], mean), 3)
  )
}
print(round(temp,3))
