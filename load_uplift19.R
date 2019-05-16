library(data.table)
library(caret)

load_uplift19 <- function(path){
  data <- fread(path, check.names = TRUE)
  # Keep only campaigns that assign a 15 absolute amount coupon for consistency
  data <- data[campaignValue==15 & campaignUnit=="CURRENCY", ]
  
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
  
  data <- data.table(predict(caret::dummyVars(~., data, sep="_", fullRank=TRUE), data))
  
  TARGET = "converted"
  PROFIT = "checkoutAmount"
  W = "treatmentGroup"
  return(list("X"=data[,!c(TARGET, PROFIT, W), with=FALSE], "Y"=data[[TARGET]], "value"=data[[PROFIT]], "W"=data[[W]]))
}
