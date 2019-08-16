# #### Static Quantile-based Mapping ####
# map_propensity <- function(model_score, target_ratio, groups=9){
#   # Platt scaling
#   # platt_scaler <- glm(y~prob, family=binomial(link='logit'), 
#   #                     weights = ifelse(exp$none$y==1, 1/(mean(exp$none$y==1)), 1),
#   #                     data = data.frame(cbind('y'=exp$none$y,'prob'=churn_pred))
#   #                     )
#   # treat_prob <- predict(platt_scaler, newdata=data.frame(prob=churn_pred), type='response')
#   # Cut into groups based on the score quantiles
#   score_quantiles <- quantile(model_score,seq(0,1,1/groups))
#   
#   if(length(score_quantiles) != length(unique(score_quantiles))){
#     score_quantiles <- unique(score_quantiles)
#     warning(paste("Reduced number of groups to", length(score_quantiles)-1, "because scores between groups were constant"))
#   }
#   
#   model_score <- cut(model_score,breaks=score_quantiles,labels=FALSE, include.lowest = TRUE)-1
#   
#   # Adjust to expected target ratio by shifting min or max
#   if(target_ratio<=0.5){
#     new_max <- 2*target_ratio - 0.05
#     model_score <- 0.05 + model_score* (new_max-0.05)/(groups-1)
#   }
#   if(target_ratio>0.5){
#     new_min <- 2*target_ratio - 0.95
#     model_score <- new_min + model_score* (0.95-new_min)/(groups-1)
#   }
#   
#   return(model_score)
# }


#### Static Quantile-based Mapping ####
get_score_quantiles_ <- function(model_score, groups=9){
  # Cut into groups based on the score quantiles
  score_quantiles <- quantile(model_score,seq(0,1,1/groups))
  
  if(length(score_quantiles) != length(unique(score_quantiles))){
    score_quantiles <- unique(score_quantiles)
    warning(paste("Reduced number of groups to", length(score_quantiles)-1, "because scores between groups were constant"))
  }
  
  return(score_quantiles)
}

map_propensity_quantiles <- function(model_score, target_ratio, score_breaks){
  # Platt scaling
  # platt_scaler <- glm(y~prob, family=binomial(link='logit'), 
  #                     weights = ifelse(exp$none$y==1, 1/(mean(exp$none$y==1)), 1),
  #                     data = data.frame(cbind('y'=exp$none$y,'prob'=churn_pred))
  #                     )
  # treat_prob <- predict(platt_scaler, newdata=data.frame(prob=churn_pred), type='response')

  model_score <- cut(model_score, breaks=score_breaks,labels=FALSE, include.lowest = TRUE)-1
  
  n_groups <- length(unique(model_score))
  
  # Adjust to expected target ratio by shifting min or max
  if(target_ratio<=0.5){
    new_max <- 2*target_ratio - 0.05
    model_score <- 0.05 + model_score* (new_max-0.05)/(n_groups-1)
  }
  if(target_ratio>0.5){
    new_min <- 2*target_ratio - 0.95
    model_score <- new_min + model_score* (0.95-new_min)/(n_groups-1)
  }
  
  return(model_score)
}


#### Generalized Logistic Mapping ####

generalized_logistic_ <- function(x, A=0, K=1, B=1, nu=1, Q=1, C=1){
  # B: growth rate
  # nu: skew
  return(A + (K-A)/((C+Q*exp(-B*x))^(1/nu)))
}

map_propensity_logistic <- function(model_score, min_score, max_score, target_ratio=NULL){
  if(!is.null(target_ratio)){
    stop("target_ratio not implemented")
  }
  
  # Min-Max scaling based on e.g. training data model scores
  model_score_normalized <- -1 + (model_score - min_score)*(2)/(max_score-min_score)
  
  prop_score <- generalized_logistic_(model_score_normalized, A=0.05, K=0.95, B=5)
  
  return(prop_score)
}
