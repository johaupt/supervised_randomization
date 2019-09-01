# Test ------------------------------
#path <- "../Data/DMC01"
#temp <- load_dmc01(path)
# End Test --------------------------


library(data.table)


load_dmc01 <- function(folder){
  data_train <- fread(paste0(path,"/dmc2001-train.txt"))
  data_class <- fread(paste0(path,"/dmc2001-class.txt"))
  y_class <- fread(paste0(path,"/dmc2001-trueclass.txt"))
  
  data_class$AKTIV <- y_class$AKTIV[match(y_class$ID, data_class$ID)]
  data <- rbind(data_train, data_class)
  
  # Change target variable coding to 1 (purchase) and 0 (no future purchase)
  data$AKTIV <- ifelse(data$AKTIV==1,0,1)
  
  # Recode West/East variable
  data[WO=="F", WO := NA]

  # factor variables
  for (var in c("WO","Regiotyp","Bebautyp", "Strtyp")){
    set(data, j=var, value=paste("level",data[[var]],sep="_"))
    #set(data, which(is.na(data[[var]])), var, "UNKNOWN")
    set(data, j=var, value=factor(data[[var]]))
  }
  
  # Merge rare levels for robustness over cross validation
  data[Bebautyp=="level_0", Bebautyp:="level_1"]
  data[Bebautyp=="level_5", Bebautyp:="level_4"]
  
  data[Strtyp=="level_5", Strtyp:="level_4"]
  
  # Impute missing with median and include missing indicator
  # Transform to factor
  data[, "missing_habitation":=ifelse(data$Regiotyp == "level_NA", 1,0)]
  
  for (var in c("WO","Regiotyp","Bebautyp", "Strtyp")){
    set(data, i=which(data[[var]] == "level_NA"), j=var, value= names(which.max(table(data[[var]]))) )
    set(data, j=var, value= factor(data[[var]]))
  }
  
  # Length of customer relationship
  data[, customer_age := 2001-jahrstart]
  data[customer_age>10,customer_age := 0] # Correct for no history
  data[, jahrstart := NULL]
  
  # Drop ID
  data[, ID := NULL]
  
  # Median imputation for missing values
  for (var in colnames(data)[which(sapply(data, is.numeric))]){
    temp_median <- median(data[[var]], na.rm=TRUE)
    if(!is.integer(temp_median)){
      set(data, j=var, value=as.numeric(data[[var]]))
    }
    
    set(data, i=which(is.na(data[[var]])), j=var, value=median(data[[var]], na.rm=TRUE))
  }
  
  # Drop due to inconsistent scaling between regions
  data[, Kaufkraft := NULL]
  
  # Rename target variable
  setnames(data, "AKTIV", "y")
  
  # Normalize non-ordinal
  for(var in c("AnzHH", "AnzGew", "customer_age")){
    set(data, j=var, value=scale(data[[var]]))
  }
  
  # One-hot encode data
  data <- data.frame(model.matrix(~.-1,data))

  return(data)

}
