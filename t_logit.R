#### Two-Learner Approach with Logit as its base Learner ####
# Includes the standard version of the logit two-learner and 
# the double robust version

#### Two-model Logit ####
T_Logit <- function(X, y, W, prop_score=NULL){
  data <- data.frame(X,".y"=y)
  if(is.null(prop_score)){
    prop_score = rep(1, nrow(data))
  }

  # Train model for W0
  logit0 <- glm(.y~., family=binomial(link='logit'),
                  weights = 1/(1-prop_score[W==0]),
                  data = data[W==0,])
  
  # Train model for W1
  logit1 <- glm(.y~., family=binomial(link='logit'),
                weights = 1/prop_score[W==1],
                data = data[W==1,])
  
  t_logit <- list(model0=logit0, model1=logit1)
  class(t_logit) <- c("tlearner")
  return(t_logit)
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
  
  return(t_logit_dr)
}

