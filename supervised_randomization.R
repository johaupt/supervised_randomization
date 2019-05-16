## Individual biased random treatment (based on propensity)
map_propensity <- function(model_score, target_ratio, groups=9){
  # Platt scaling
  # platt_scaler <- glm(y~prob, family=binomial(link='logit'), 
  #                     weights = ifelse(exp$none$y==1, 1/(mean(exp$none$y==1)), 1),
  #                     data = data.frame(cbind('y'=exp$none$y,'prob'=churn_pred))
  #                     )
  # treat_prob <- predict(platt_scaler, newdata=data.frame(prob=churn_pred), type='response')
  # Cut into groups based on the score quantiles
  score_quantiles <- quantile(model_score,seq(0,1,1/groups))
  
  if(length(score_quantiles) != length(unique(score_quantiles))){
    score_quantiles <- unique(score_quantiles)
    warning(paste("Reduced number of groups to", length(score_quantiles)-1, "because scores between groups were constant"))
  }
  
  model_score <- cut(model_score,breaks=score_quantiles,labels=FALSE, include.lowest = TRUE)-1
  
  # Adjust to expected target ratio by shifting min or max
  if(target_ratio<=0.5){
    new_max <- 2*target_ratio - 0.05
    model_score <- 0.05 + model_score* (new_max-0.05)/(groups-1)
  }
  if(target_ratio>0.5){
    new_min <- 2*target_ratio - 0.95
    model_score <- new_min + model_score* (0.95-new_min)/(groups-1)
  }
  
  return(model_score)
}