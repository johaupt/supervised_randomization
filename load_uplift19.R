library(data.table)

load_uplift19 <- function(path){
  data <- fread(path, check.names = TRUE)
  # Keep only campaigns that assign a 15 absolute amount coupon for consistency
  data <- data[campaignValue==15 & campaignUnit=="CURRENCY", ]
  
  # Keep only view counts<=60
  data <- data[targetViewCount<=60, ]
  
  # Rename Treatment Indicator
  data$treatmentGroup <- abs(data$controlGroup-1)
  data$controlGroup <- NULL
  
  # Ignore unrelated discounts by adding them back to the basket
  # Simultanously add back the coupon cost of 15, where it applied
  data[, checkoutAmount := checkoutAmount+checkoutDiscount]
  data[, checkoutDiscount := NULL]
  
  # Drop variables
  dropVar <- c(
    # Index or and campaign descriptors
    "epochSecond","campaignId", "campaignUnit", "campaignValue", "campaignTags",
    # Transformed outcome variable
    "label",
    # Cancellation during conversion process?
    "dropOff",
    # Confirm or Abort on banner that treatment was seen?
    "confirmed", "aborted",
    # Duplicate variables         
    "HoursSinceOldestSession",
    # Baseline dummy for time of day
    "IsEarlyMorning"
    )
  data[, (dropVar) := NULL]
  
  # Zip Code Simplification
  data[,ZipCode := as.character(ZipCode)]
   # All German Zip Codes have 5 digits 
  data[nchar(ZipCode) != 5, ZipCode := "Foreign"]
  # The first digits indicates the most general zone
  data[ZipCode != "Foreign", ZipCode := substring(ZipCode,1,1)]
  
  # Factor variables
  data[, ZipCode := factor(ZipCode)]
  data[, DeviceCategory := factor(DeviceCategory)]

  # Drop binary variables with <0.01 of positive values
  dropVar <- c("HasConfirmedBefore",    "DidConfirmLastWeek",    "DidConfirmLastYear",    "DidConfirmLastMonth", "InitPageWas.sale.", "ChannelIs.EMAIL.", "ChannelIs.PAID.", "ChannelIs.SOCIAL.",
               "ScreenTypeIs.home.",    "ScreenTypeIs.sale.",    "ScreenTypeIs.account.", "ScreenTypeIs.about.")
  data[, (dropVar) := NULL]
  # Drop constant variables
  data <- data[, sapply(data, function(x) length(unique(x))) != 1, with=FALSE]
  
  # Drop variables with a correlation >0.9
  # Multicollinearity
  # metavars_left = c('converted', 'checkoutAmount', 'treatmentGroup')
  # data = data[c(metavars_left, setdiff(names(data),metavars_left))]
  # num <- c(4:ncol(data)) # disregard meta-variables
  # cor_num_P <- cor(data[,num],use="pairwise.complete.obs", method="pearson")
  # #corrplot(corr = cor_num_P, method="pie", tl.cex=0.3, diag=F,type="lower")
  # corr_var <- findCorrelation(cor_num_P, cutoff = .90, verbose= TRUE, names= TRUE, exact = TRUE)
  # data <- data[ , !(names(data) %in% corr_var)]
  dropVarColl <- c("DidConvertLastYear", "DidVisitLastYear", "ViewedBefore.overview.", 
               "log.of.SecondsSinceFirst.home.", "log.of.SecondsSinceFirst.account.", 
               "log.of.SecondsSinceFirst.about.", "log.of.SecondsSinceFirst.search.", 
               "log.of.SecondsSinceFirst.overview.", "log.of.SecondsSinceFirst.product.", 
               "log.of.SecondsSinceFirst.sale.", "log.of.ViewsOn.overview.", 
               "log.of.NumberOfDifferent.sale.", "log.of.NumberOfDifferent.search.", 
               "log.of.NumberOfDifferent.about.", "ClientKnown", "log.of.SecondsSinceOn.account.", 
               "log.of.SecondsSinceOn.about.", "log.of.SecondsSinceOn.sale.", 
               "log.of.SecondsSpentOn.account.", "log.of.SecondsSpentOn.about.", 
               "log.of.SecondsSpentOn.overview.", "log.of.SecondsSinceOn.home.", 
               "log.of.SecondsSpentOn.sale.", "ScreenTypeIs.product.", "ViewedBefore.account.", 
               "ViewedBefore.about.", "log.of.ViewsOn.product.", "log.of.ViewsOn.sale.", 
               "ViewedBefore.home.", "log.of.ViewsOn.home.", "log.of.ViewsOn.account."
  )
  data[, (dropVarColl) := NULL]
  
  data <- data.table(predict(caret::dummyVars(~., data, sep="_", fullRank=TRUE), data), keep.rownames = FALSE)
  
  # Missing Values
  #apply(data, 2, function(x) any(is.na(x))) # no missing values
  # replace negative values with zero since variables should be non-zero
  #data <- as.data.frame(lapply(data, function(x){replace(x, x <0,0)})) # no negative values

  # # Multivariate outliers
  # mod <- lm(converted ~ ., data=data)
  # cooksd <- cooks.distance(mod)
  # sample_size <- nrow(data)
  # plot(cooksd, pch="*", cex=2, main="Outliers as determined by Cooks distance")
  # abline(h = 4/sample_size, col="blue")
  # text(x=1:length(cooksd)+1, y=cooksd, labels=ifelse(cooksd>4/sample_size, names(cooksd),""), col="blue")
  # outl <- as.numeric(names(cooksd)[(cooksd > (4/sample_size))])
  # data <- data[-outl, ] # removed 5000 observations that had outliers
  
  TARGET = "converted"
  PROFIT = "checkoutAmount"
  W = "treatmentGroup"
  return(list("X"=data[,!c(TARGET, PROFIT, W), with=FALSE], "Y"=data[[TARGET]], "value"=data[[PROFIT]], "W"=data[[W]]))
}

#### TESTS ####
#data <- load_uplift19("../data/explore.csv")
#str(data,1)
#lm <- glm(w~., cbind(data$X, "w"=data$W), family = "binomial")
#ModelMetrics::auc(actual=data$W, predict(lm, type='response'))
