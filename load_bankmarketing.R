#### Debugging
#path <- "/Users/hauptjoh/Data/UCI/bank_marketing/bank-additional/bank-additional-full.csv"
####
library(data.table)

load_bankmarketing <- function(path){
    data <- fread(path, stringsAsFactors = FALSE)  
    
    # Duration is actual call duration during targeting
    # previous and pdays are only available for few customers
    dropList <- c("duration", "pdays", "previous")
    data[, (dropList) := NULL]
    
    # Drop very rare (<0.005) classes in 'marital', 'education' and 'default' to avoid issues during cross validation
    data <- data[marital!='unknown' & education != 'illiterate' & default != 'yes', ]
    
    # Transform string outcome
    data[, y:=ifelse(y=="yes",1,0)]
    
    # Standardize
    numVars <- c("age", "campaign", "emp.var.rate", "cons.price.idx", "cons.conf.idx", 
                 "euribor3m", "nr.employed")
    data <- data.frame(data)
    data[, numVars] <- sapply(data[, numVars], function(x) (x-mean(x))/sd(x) )
    
    # One-hot/dummy encode categorical variables
    data <- data.frame(model.matrix(~.-1, data))
    
    return(data)
}
