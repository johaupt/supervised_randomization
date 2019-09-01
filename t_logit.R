#### Two-Learner Approach with Logit as its base Learner ####
# Includes the standard version of the logit two-learner and 
# the double robust version

library(rpart)


#### Two-model Modified Outcome Logit ####
MOM_Logit <- function(X, y, W, prop_score){
  data <- data.frame(X,".y"=y)
  
  # # Train model for W0
  # logit0 <- glm(.y~., family=binomial(link='logit'),
  #               data = data[W==0,])
  # 
  # # Train model for W1
  # logit1 <- glm(.y~., family=binomial(link='logit'),
  #               data = data[W==1,])
  # 
  # mu0 <- predict(logit0, X, type='response')
  # mu1 <- predict(logit1, X, type='response')
  # Train model for W0
  logit0 <- lm(.y~., data = data[W==0,])
  
  # Train model for W1
  logit1 <- lm(.y~., data = data[W==1,])
  
  mu0 <- predict(logit0, X)
  mu1 <- predict(logit1, X)
  
  y_dr <- mu1 - mu0 + (W/prop_score)*(y-mu1) - ((1-W)/(1-prop_score))*(y-mu0)
  
  tau_reg <- lm(.y~., data=data.frame(".y"=y_dr,X))
  
  return(tau_reg)
}



#### Two-model Logit ####
T_Logit <- function(X, y, W, prop_score=NULL){
  data <- data.frame(X,".y"=y)
  if(is.null(prop_score)){
    weights0 = NULL
    weights1 = NULL
  }else{
    weights0 = (1/(1-prop_score[W==0])) 
    weights1 = (1/prop_score[W==1])
  }

  # Train model for W0
  logit0 <- glm(.y~., family=binomial(link='logit'),
                  weights = weights0,
                  data = data[W==0,])
  
  # Train model for W1
  logit1 <- glm(.y~., family=binomial(link='logit'),
                weights = weights1,
                data = data[W==1,])
  
  t_logit <- list(model0=logit0, model1=logit1)
  class(t_logit) <- c("tlearner")
  return(t_logit)
}

#### Two-model Reg ####
T_Regression <- function(X, y, W, prop_score=NULL){
  data <- data.frame(X,".y"=y)
  if(is.null(prop_score)){
    weights0 = NULL
    weights1 = NULL
  }else{
    weights0 = (1/(1-prop_score[W==0])) 
    weights1 = (1/prop_score[W==1])
  }
  
  # Train model for W0
  reg0 <- lm(.y~., 
                weights = weights0,
                data = data[W==0,])
  
  # Train model for W1
  reg1 <- lm(.y~., 
                weights = weights1,
                data = data[W==1,])
  
  t_reg <- list(model0=reg0, model1=reg1)
  class(t_reg) <- c("tlearner")
  return(t_reg)
}

# T-Tree
T_rpart <- function(X, y, W, prop_score=NULL, method='class'){
  data <- data.frame(X,".y"=y)
  if(is.null(prop_score)){
    weights0 = NULL
    weights1 = NULL
  }else{
    weights0 = (1/(1-prop_score[W==0])) 
    weights1 = (1/prop_score[W==1])
  }
  
  # Train model for W0
  model0 <- rpart(.y~., data = data[W==0,], 
                  weights = weights0, method=method)
  
  # Train model for W1
  model1 <- rpart(.y~., data = data[W==1,],
                weights = weights1,
                method=method)
  
  t_learner <- list(model0=model0, model1=model1)
  class(t_learner) <- c("tlearner")
  return(t_learner)
}


# Two-Model nnet
T_NNet <- function(X, y, W, prop_score=NULL, trace=FALSE){
  data <- data.frame(X,".y"=y)
  if(is.null(prop_score)){
    weights0 = NULL
    weights1 = NULL
  }else{
    weights0 = (1/(1-prop_score[W==0])) 
    weights1 = (1/prop_score[W==1])
  }
  
  # Train model for W0
  model0 <- nnet(.y~., size=ncol(X),
                weights = weights0,
                data = data[W==0,], 
                entropy=TRUE, maxit=500, MaxNWts=5000, decay=0.001, reltol = 1e-6, trace=trace)
  
  # Train model for W1
  model1 <- nnet(.y~., size=ncol(X),
                weights = weights1, 
                data = data[W==1,],
                entropy=TRUE, maxit=500, MaxNWts=5000, decay=0.001, reltol= 1e-6, trace=trace)
  
  t_learner <- list(model0=model0, model1=model1)
  class(t_learner) <- c("tlearner")
  return(t_learner)
}

# Prediction function for two-model class 
predict.tlearner <- function(object, newdata, ...){
  pred_diff <- predict(object$model1, newdata, ...) - 
               predict(object$model0, newdata, ...)
  
  if(!is.null(ncol(pred_diff))){
    pred_diff <- pred_diff[,ncol(pred_diff)]
  }
  return(pred_diff)
}



# Prediction function for two-model class 
predict.tlearner <- function(object, newdata, ...){
  pred_diff <- predict(object$model1, newdata, ...) - 
    predict(object$model0, newdata, ...)
  
  if(!is.null(ncol(pred_diff))){
    pred_diff <- pred_diff[,ncol(pred_diff)]
  }
  return(pred_diff)
}


#### Two-Model Logit Double Robust ####
T_Logit_DR <- function(X, y, W, prop_score=NULL){
  data <- data.frame(X,"y"=y)
  if(is.null(prop_score)){
    prop_score = rep(0.5, nrow(data))
  }

  y0 <- glm(y~., family=binomial(link='logit'),
            data = data[W==0,])
  
  y0_hat <- predict(y0,X, type="response")
  
  y1 <- glm(y~., family=binomial(link='logit'),
            data = data[W==1,])
  
  y1_hat <- predict(y1,X,type="response")
  
  y_mo <- (y1_hat - y0_hat) + (W*(y-y1_hat))/prop_score - ((1-W)*(y-y0_hat)/(1-prop_score))
  
  t_logit_dr <- glm(y_mo~.,data= cbind(X,"y_mo"=y_mo) )
  
  class(t_logit) <- c("tlearner")
  return(t_logit_dr)
}

