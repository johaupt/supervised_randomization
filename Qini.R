# # Test case
# W <- rbinom(1000, 1, 0.5)
# y <- rnorm(1000, 0,1)
# tau <- rnorm(1000, 0.5,0.1)
# Y <- as.numeric(1/(1+exp(-y - W*tau))>0.5)
# scores <- tau + rnorm(1000,0,0.2)
# 
# # Debuggin test
# uplift::qini(uplift::performance(scores, rep(0, length(scores)), Y, W))
# qini_score(scores, Y, W, plotit=TRUE)


qini_score <- function(scores, Y, W, groups = 10, plotit=FALSE){
  if(length(unique(lengths(list(scores, Y, W)))) != 1){
    stop("input scores, Y, W must have same length")
  }
  
  mm <- cbind(dif.pred = scores, y = Y, ct = W, dif.pred_r = rank(-scores))
  bk <- unique(quantile(mm[, 4], probs = seq(0, 1, 1/groups)))
  if ((length(bk) - 1) != groups) 
    warning("uplift: due to ties in uplift predictions, the number of groups is less than ", 
            groups)
  mm <- cbind(mm, decile = cut(mm[, 4], breaks = bk, labels = NULL, 
                               include.lowest = TRUE))
  
  y1 <- tapply(mm[,2], INDEX = list(mm[,3],mm[,5]),sum)
  n <- tapply(mm[,2], INDEX = list(mm[,3],mm[,5]),length)
  r.y1 <- y1/n
  
  uplift <- r.y1["1",] - r.y1["0",]
  # Relative increase in gains incrementally including group after group
  inc.gains <- cumsum( (y1["1",]/sum(n["1",])) - y1["0",]/sum(n["0",]) )
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
                                                                 maxy), ...)
    lines(random.inc.gains * 100 ~ seq(100/groups, 100, 100/groups), 
          type = "l", col = "red", lty = 1)
    #legend("topright", c("Model", "Random"), col = c("blue", 
    #                                                 "red"), lty = c(2, 1))
  }
  
  return(unname(qini))
}

