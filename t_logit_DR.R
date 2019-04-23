T_Logit_DR <- function(X, y, W, prop_score=NULL){
  data <- data.frame(X,"y"=y)
  if(is.null(prop_score)){
    pr = rep(1, nrow(data))
  }
  
  
  
  #pr <- glm(W~., family=binomial(link='logit'),
              #data = data[ , -which(names(data) %in% c(".y"))])
  
  y0 <- glm(y~., family=binomial(link='logit'),
                data = data[W==0,])
  
  y0_hat <- predict(y0,X, type="response")
  
  y1 <- glm(y~., family=binomial(link='logit'),
                data = data[W==1,])
  
  y1_hat <- predict(y1,X,type="response")
  
  y_mo <- (y1_hat - y0_hat) + (W*(y-y1_hat))/prop_score - ((1-W)*(y-y0_hat)/(1-prop_score))
  
  t_logit_dr <- glm(y_mo~.,data=data)
    
  
  
  return(t_logit_dr)
  
  
}




        