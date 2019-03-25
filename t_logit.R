T_Logit <- function(X, y, W, prop_score=NULL){
  data <- data.frame(X,".y"=y)
  if(is.null(prop_score)){
    prop_score = rep(1, nrow(data))
  }

  
  logit0 <- glm(.y~., family=binomial(link='logit'),
                  weights = 1/(1-prop_score[W==0]),
                  data = data[W==0,])
  
  logit1 <- glm(.y~., family=binomial(link='logit'),
                weights = 1/prop_score[W==1],
                data = data[W==1,])
  
  t_logit <- list(model0=logit0, model1=logit1)
  class(t_logit) <- c("tlearner")
  return(t_logit)
}

predict.tlearner <- function(object, newdata, ...){
  predict(object$model1, newdata, type="response", ...) - 
    predict(object$model0, newdata,type="response", ...)
}
