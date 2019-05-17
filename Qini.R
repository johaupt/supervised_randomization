# # Test case (standard)
# W <- rbinom(1000, 1, 0.5)
# y <- rnorm(1000, 0,1)
# tau <- rnorm(1000, 0.5,0.1)
# Y <- as.numeric(1/(1+exp(-y - W*tau))>0.5)
# scores <- tau + rnorm(1000,0,0.2)
# groups <- 10
# p_treatment <- 0.5
# 
# # Test case (propensity scores)
# y <- rnorm(1000, 0,1)
# tau <- rnorm(1000, 0.5,0.1)
# scores <- tau + rnorm(1000,0,0.2)
# 
# # Treatment with 0.8 if tau>0.5 and 0.2 if tau<0.5
# W <- foreach(i=1:length(y), .combine=c)%do%{ 
#   if(tau[i]>0.5){
#     rbinom(1, 1, 0.8)
#   }else{
#     rbinom(1, 1, 0.2)
#   }
# }
# p_treatment <- ifelse(tau>0.5, 0.8,0.2)
# Y <- as.numeric(1/(1+exp(-y - W*tau))>0.5)
# groups <- 10
# 
# # Debugging test
#  uplift::qini(uplift::performance(scores, rep(0, length(scores)), Y, W))
#  qini_score(scores, Y, W, p_treatment = p_treatment, plotit=TRUE)

 
# DEPRECATED
#calc_ATE <- function(y, g, prop_score){
#  if(length(prop_score)==1){
#    prop_score = rep(prop_score, length(y))
#  }
#  return( sum(y*g/prop_score)/sum(g*1/prop_score) - sum(y*(1-g)/(1-prop_score))/sum((1-g)*1/prop_score)  )
#  #return( (sum(y*g/prop_score) - sum(y*(1-g)/(1-prop_score)) ) /length(y) )
#}

library(data.table)

### QINI SCORE (propensity corrected)
qini_score <- function(scores, Y, W, p_treatment=0.5, groups = 10, plotit=FALSE){
  if(length(unique(lengths(list(scores, Y, W)))) != 1){
    stop("input scores, Y, W must have same length")
  }
  
  # mm matrix to contain scores, deciles, observed response and experiment group indicator
  mm <- cbind(tau_hat = scores, y = Y, ct = W, tau_hat_rank = rank(-scores), prop_score = p_treatment)
  bk <- unique(quantile(mm[, "tau_hat_rank"], probs = seq(0, 1, 1/groups)))
  if ((length(bk) - 1) != groups) 
    warning("uplift: due to ties in uplift predictions, the number of groups is less than ", 
            groups)
  mm <- cbind(mm, decile = cut(mm[, "tau_hat_rank"], breaks = bk, labels = NULL, 
                               include.lowest = TRUE))
  
  mm <- data.table(mm)
  setorder(mm, ct,decile)
  
  # ATE per decile
  #uplift <- mm[, .(uplift=calc_ATE(y = y, g = ct, prop_score = prop_score)),by=decile]
  
  # No. of positive responses and group size
  #deciles <- mm[, .(y1=sum(y), n=.N),by=.(decile, ct)]
  # No. of positive responses and group size, corrected for propensity score
  deciles <- mm[, .(y1=sum(y/prop_score), n=sum(1/prop_score)),by=.(decile, ct)]
  # Ratio of positve responses in each decile relative to total number of positive responses in exp. group
  deciles[, inc.y1.ratio := y1/sum(n), by=ct]
  # Incremental gain in positive responses (as difference to control group)
  inc.gains <- cumsum(deciles[ct==1, inc.y1.ratio] - deciles[ct==0, inc.y1.ratio])
  
  # # No. of positive responses per decile per experiment group
  # y1 <- tapply(mm[,"y"], INDEX = list(mm[,"ct"],mm[,"decile"]),sum)
  # # No. of observations per decile per experiment group
  # n <- tapply(mm[,"y"], INDEX = list(mm[,"ct"],mm[,"decile"]),length)
  # # Success ratio per decile per experiment group
  # r.y1 <- y1/n
  # # Uplift as difference in success ratio per decile (ATE per decile)
  # uplift <- r.y1["1",] - r.y1["0",]
  # 
  # # Relative increase in gains incrementally including group after group
  # inc.gains <- cumsum( (y1["1",]/sum(n["1",])) - y1["0",]/sum(n["0",]) )
  
  ####
  overall.inc.gains <- inc.gains[length(inc.gains)]
  
  random.inc.gains <- cumsum(rep(overall.inc.gains/groups, 
                                 groups))
  
  x_axis <- seq(1/groups, 1, 1/groups)
  y_axis <- inc.gains
  
  auuc <- 0
  for (i in 2:length(x_axis)) {
    auuc <- auuc + 0.5 * (x_axis[i] - x_axis[i-1]) * (y_axis[i] + y_axis[i-1])
  }
  
  y_axis.rand <- random.inc.gains
  auuc.rand <- 0
  for (i in 2:length(x_axis)) {
    auuc.rand <- auuc.rand + 0.5 * (x_axis[i] - x_axis[i-1]) * (y_axis.rand[i] + y_axis.rand[i - 1])
  }
  qini <- auuc - auuc.rand
  
  if(plotit){
    miny <- 100 * min(c(random.inc.gains, inc.gains))
    maxy <- 100 * max(c(random.inc.gains, inc.gains))
    
    plot(inc.gains * 100 ~ seq(100/groups, 100, 100/groups), 
         type = "b", col = "blue", lty = 2, xlab = "Proportion of population targeted (%)", 
         ylab = "Cumulative incremental gains (pc pt)", ylim = c(miny, 
                                                                 maxy) )#, ...)
    lines(random.inc.gains * 100 ~ seq(100/groups, 100, 100/groups), 
          type = "l", col = "red", lty = 1)
    #legend("topright", c("Model", "Random"), col = c("blue", 
    #                                                 "red"), lty = c(2, 1))
  }
  
  return(unname(qini))
}

