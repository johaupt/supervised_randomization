#### Packages ####
#install.packages("pacman")
library(pacman)
pacman::p_load("ggplot2","reshape2","caret","drtmle","grf","foreach","uplift","data.table","ModelMetrics", "plyr", "parallel", "doParallel", "SuperLearner", "tidyverse", "nnet", "xtable")

source("supervised_randomization.R")
source("t_logit.R")
source("Qini.R")
source("costs.R")

RNGkind("L'Ecuyer-CMRG")

calc_ATE <- function(y, w, prop_score){
  #if(is.null(prop_score)){
  #  prop_score = rep(0.5, length(y))
  #}
  return( (sum(y*w/prop_score) - sum(y*(1-w)/(1-prop_score)) ) /length(y) )
}


#### Load Real data ####

# Results path
results_path = "../results/bankmarketing/"

source("load_bankmarketing.R")
path = "../data/bank-additional-full.csv"
data <- load_bankmarketing(path)

# source("load_dmc01.R")
# path <- "../data/"
# data <- load_dmc01(path)


tau_model <- function(X, hidden_layer, ATE){
  X <- as.matrix(X)
  no_var <- ncol(X)
  W1 <- matrix(rnorm(no_var*hidden_layer, 0, 1), nrow=no_var, ncol=hidden_layer)
  W2 <- matrix(rnorm(hidden_layer,        0, 1), nrow=hidden_layer, ncol=1)
  
  h <- X %*% W1
  h <- dlogis(X %*% W1)
  
  o <- c(h %*% W2)
  o <- o / (25*sd(o))
  o <- o - mean(o) + ATE
  return(o)
}

X <- data[,!(colnames(data)=="y")]

# Simulation bankmarketing
X_tau <- X[,c(1:24, 43:48)]
set.seed(123456789)
TAU <- tau_model(X_tau, hidden_layer=ncol(X_tau), ATE=0.05)
X <- X[,!(colnames(X) %in% c("age", "maritalmarried", "maritalsingle"))]

# # Simulation DMC01 catalogue
# X_tau <- X
# TAU <- tau_model(X_tau, hidden_layer=ncol(X_tau), ATE=0.05)
# X <- X[,!grepl("Typ*|PHARM*", colnames(X))]

quantile(TAU, probs = seq(0,1,0.025))


R <- rbinom(nrow(data), size = 1, prob = abs(TAU))
Y <- cbind(Y0=data$y, Y1=data$y)

for(i in 1:nrow(Y)){
  if(R[i] == 1){
      if(TAU[i]>0) Y[i, ] <- c(0,1)
      if(TAU[i]<0) Y[i, ] <- c(1,0)
  }
}

# Sanity check: ATE should match target ATE
ATE <- mean(Y[,2]) - mean(Y[,1])
ATE

#### Existing Targeting Model ####
# Vary error variance to make the existing model better or worse
# for(e in seq(0,0.05, 0.005)){
#     tau_hat <- TAU + rnorm(length(tau), 0, e)
#     print(cor(tau_hat, tau, method = "spearman"))
# }


#### Create experiments ####
do_experiment <- function(Y1, Y0, prop_score){
     # Sample treatment group based on propensity
     W <- rbinom(length(Y1), 1, prop_score)
     Y <- ifelse(W==1, Y1, Y0)
     return(list('treatment'=W, 'outcome'=Y))
}


#### SET PARAMETERS ####

# Cross-Validation split
NO_FOLDS = 4
NO_EXPERIMENT_ITER = 50
IMBALANCED_EXP_RATIO = 0.666

DEFAULT_TARGETING_ERROR = 0.025

#set.seed(123456789)
test_idx_list <- caret::createFolds(Y[,2] - Y[,1], k = NO_FOLDS, list = TRUE, returnTrain = FALSE)

# Model error
UPLIFT_MODEL_ERROR <- c(0.025, 0.04, 0.05, 0.08) #c(0, 0.005, 0.01, 0.025, 0.04, 0.05, 0.08, 0.1, 0.15)
uplift_model_error_list <- foreach(e=UPLIFT_MODEL_ERROR)%do% rnorm(length(TAU), 0, e)

# Profit setting
CONTACT_COST = 1 # Contact costs
OFFER_COST = 0 # Price reduction

VALUE = c(seq(5, 80, 5))

exp_grid = data.frame(expand.grid(list(
  "uplift_model_error_idx" = 1:length(uplift_model_error_list),
  "testset"=1:length(test_idx_list),
  "iteration"=1:NO_EXPERIMENT_ITER,
  "randomization" = c("full", "imbalanced", "supervised")
)))
# The uplift model error is only relevant for supervised randomization
exp_grid <- exp_grid[ !((exp_grid$uplift_model_error_idx != 1) & (exp_grid$randomization != "supervised")), ]


#### Register cluster for parallel processing ####
cl <- makeCluster(25)
registerDoParallel(cl)
RNGkind("L'Ecuyer-CMRG") 
clusterSetRNGStream(cl,iseed = 123456789)


#### Experiment Cost Summary  ####
exp_grid_ABtest <- rbind(exp_grid, 
                         data.frame(expand.grid(list(
                           "uplift_model_error_idx" = 1,
                           "testset"=1:length(test_idx_list),
                           "iteration"=1:NO_EXPERIMENT_ITER,
                           "randomization" = c("none", "all")
                         )))
)

# For each quality level of existing model
result_exp <- foreach(exp_idx = 1:nrow(exp_grid_ABtest), .inorder=TRUE, .combine=rbind)%dopar%{
          # Select experiment from grid
          randomization          <- exp_grid_ABtest$randomization[exp_idx]
          uplift_model_error_idx <- exp_grid_ABtest$uplift_model_error_idx[exp_idx]
          testset                <- exp_grid_ABtest$testset[exp_idx]
          i                      <- exp_grid_ABtest$iteration[exp_idx]
          
          # Test set observations
          test_idx <- test_idx_list[[testset ]]
        
          # Simulate existing uplift model
          tau_hat <- TAU + uplift_model_error_list[[uplift_model_error_idx ]]
          # Get score quantiles from training set for mapping 
          score_quantiles_ <- get_score_quantiles_(tau_hat[-test_idx], groups=18)
          supervised_prop_score <- map_propensity_quantiles(model_score = tau_hat[-test_idx], score_breaks = score_quantiles_, target_ratio = 0.5)
          # supervised_prop_score <- map_propensity_logistic(model_score = tau_hat[-test_idx], 
          #                                                  min_score = quantile(tau_hat[-test_idx], probs=0.05), 
          #                                                  max_score = quantile(tau_hat[-test_idx], probs=0.95))
          
          if(randomization=="none") prop_score <- 0
          if(randomization=="all") prop_score <- 1
          if(randomization=="full") prop_score <- 0.5
          if(randomization=="imbalanced") prop_score <- IMBALANCED_EXP_RATIO
          if(randomization=="supervised") prop_score <- supervised_prop_score
          
          # Assign treatment
          Y1 <- Y[-test_idx, "Y1"]
          Y0 <- Y[-test_idx, "Y0"]
          experiment <- do_experiment(Y1=Y1, Y0=Y0, prop_score = prop_score)
          y <- experiment$outcome
          w <- experiment$treatment
          
          # Profit during experiment
          campaign_margin <- numeric()
          campaign_profit <- numeric()
          i = 1
          for(value in VALUE){
            campaign_margin[i] <- value
            campaign_profit[i] <- sum(w * (Y1 * value - CONTACT_COST) + (1-w)*(Y0 * value))
            i = i+1
          }
          
          # Collect results
          data.frame(list(
            "uplift_model_error"=UPLIFT_MODEL_ERROR[uplift_model_error_idx], "testset"=testset, "randomization_iteration" = i, "randomization"=randomization,
            "treatment_ratio"=mean(w), "outcome_ratio"=mean(y), "campaign_margin"=campaign_margin, "campaign_profit" = campaign_profit
            ), stringsAsFactors = FALSE) 
}

### Experiment summary statistics
library(tidyverse)
exp_summary <- result_exp[result_exp$campaign_margin==10,] %>% group_by(randomization, uplift_model_error) %>% summarize("outcome_ratio"=mean(outcome_ratio), "treatment_ratio"=mean(treatment_ratio))
row.names(exp_summary) <- NULL

print(exp_summary)

fwrite(data.table(exp_summary,keep.rownames = F), 
       paste0(results_path, "raw_experiment_summary_20190826.txt"))

randomization_values <- c("none", "full", "supervised","imbalanced", "all")
randomization_labels <- c("None", "Full", "Supervised","Full (Imb.)", "All")
exp_summary$randomization <- factor(exp_summary$randomization, randomization_values, randomization_labels)

exp_summary_table <- melt(exp_summary, id.vars = c("randomization", "uplift_model_error"))
exp_summary_table <- dcast( exp_summary_table[exp_summary$randomization!="Supervised" | exp_summary$uplift_model_error==DEFAULT_TARGETING_ERROR,], variable ~ randomization)

summary_values <- c("treatment_ratio", "outcome_ratio")
summary_labels <- c("Targeted Fraction of Customers", "Conversion Rate")
exp_summary_table$variable <- factor(mapvalues(exp_summary_table$variable, summary_values, summary_labels), summary_labels, summary_labels)
exp_summary_table <- exp_summary_table[order(exp_summary_table$variable),]

sink(paste0(results_path, "experiment_summary_table.txt"))
print(exp_summary_table, digits=3)
cat("\n\n\n")
print(xtable(exp_summary_table, digits=3), include.rownames=FALSE)
sink()

### Experiment profit
exp_profit <- result_exp %>% group_by(randomization, uplift_model_error, campaign_margin) %>% summarize("campaign_profit"=mean(campaign_profit))
row.names(exp_profit) <- NULL
exp_profit$randomization <- factor(exp_profit$randomization, randomization_values, randomization_labels)

exp_profit_table <- melt(exp_profit, id.vars = c("randomization", "uplift_model_error","campaign_margin"))
# Per person profit (divide by number of customer in training set)
exp_profit_table$value <- exp_profit_table$value / mean(sapply(test_idx_list, function(test_idx) nrow(X)-length(test_idx)))

exp_profit_table <- dcast( exp_profit_table[exp_profit$randomization!="Supervised" | exp_profit$uplift_model_error==DEFAULT_TARGETING_ERROR,], campaign_margin ~ randomization)

sink(paste0(results_path, "experiment_profit_table.txt"))
print(exp_profit_table, digits=2)
cat("\n\n\n")
print(xtable(exp_profit_table), include.rownames=FALSE)
sink()

#### ATE ####

result_ATE <- foreach(exp_idx = 1:nrow(exp_grid), .inorder=TRUE, .combine=rbind, .packages=c("drtmle", "SuperLearner"))%dopar%{
  # Select experiment from grid
  randomization          <- exp_grid$randomization[exp_idx]
  uplift_model_error_idx <- exp_grid$uplift_model_error_idx[exp_idx]
  testset                <- exp_grid$testset[exp_idx]
  i                      <- exp_grid$iteration[exp_idx]
  
  # Test set observations
  test_idx <- test_idx_list[[testset ]]
  
  # Simulate existing uplift model
  tau_hat <- TAU + uplift_model_error_list[[uplift_model_error_idx ]]
  
  # Get score quantiles from training set for mapping 
  score_quantiles_ <- get_score_quantiles_(tau_hat[test_idx], groups=18)
  supervised_prop_score <- map_propensity_quantiles(model_score = tau_hat[test_idx], score_breaks = score_quantiles_, target_ratio = 0.5)
  # supervised_prop_score <- map_propensity_logistic(model_score = tau_hat[-test_idx],
  #                                                  min_score = quantile(tau_hat[-test_idx], probs=0.05), 
  #                                                  max_score = quantile(tau_hat[-test_idx], probs=0.95))
  
  if(randomization=="none") prop_score <- rep(0, nrow(X[test_idx,]))
  if(randomization=="all") prop_score <- rep(1, nrow(X[test_idx,]))
  if(randomization=="full") prop_score <- rep(0.5, nrow(X[test_idx,]))
  if(randomization=="imbalanced") prop_score <- rep(IMBALANCED_EXP_RATIO, nrow(X[test_idx,]))
  if(randomization=="supervised") prop_score <- supervised_prop_score
  
  # Assign treatment
  experiment <- do_experiment(Y1=Y[test_idx, "Y1"], Y0=Y[test_idx, "Y0"], prop_score = prop_score)
  y <- experiment$outcome
  w <- experiment$treatment
  x <- X[test_idx, ]
    
    ### ATE ###
    ATE_temp <- list()
    
    # True ATE
    ATE_temp[["true_ATE"]] <- mean(tau_hat[test_idx])
    
    # IPW
    ATE_temp[["IPW"]] <- calc_ATE(y, w, prop_score = prop_score)

    output <- data.frame(list(
      "uplift_model_error"=UPLIFT_MODEL_ERROR[uplift_model_error_idx], "testset"=testset, "randomization_iteration" = i, "randomization"=randomization,
      "ATE_true" = ATE_temp$true_ATE, 
      "estimator"= "IPW", ATE_hat = ATE_temp$IPW
    ), stringsAsFactors = FALSE)
    
    #if(randomization=="supervised"){
    gn <- list()
    gn[[1]] <- prop_score
    gn[[2]] <- 1 - prop_score

    ATE_temp[["doubly_robust"]] <-  drtmle::ci(drtmle::drtmle(Y=y,A=w,W=x,a_0 = c(1,0),
                                            family=binomial(),
                                            stratify=TRUE,
                                            SL_Q = c("SL.glm"),
                                            SL_g = c("SL.glm"),
                                            SL_Qr = "SL.glm",
                                            SL_gr = "SL.glm",
                                            gn=gn,
                                            maxIter = 1),contrast=c(1,-1))$drtmle[1]
    
    output <- rbind(output,
          data.frame(list(
            "uplift_model_error"=UPLIFT_MODEL_ERROR[uplift_model_error_idx], "testset"=testset, "randomization_iteration" = i, "randomization"=randomization,
            "ATE_true" = ATE_temp$true_ATE, 
            "estimator"= "DR", ATE_hat = ATE_temp$doubly_robust
          ), stringsAsFactors = FALSE) 
    )
    
    #}

}
saveRDS(result_ATE, paste0(results_path, "ATE_predictions_20190821.rds"))

# Keep only one target model error rate
result_ATE_table <- result_ATE[result_ATE$randomization!="supervised" | result_ATE$uplift_model_error==DEFAULT_TARGETING_ERROR,]
# Sort the boxplots by sorting the factor variable
result_ATE_table$estimator <- factor(result_ATE_table$estimator, levels = c("IPW","DR"), labels=c("IPW", "DR"))

# Create the boxplot
ATE_plot <- ggplot(data = result_ATE[result_ATE_table$estimator=="DR"|result_ATE_table$randomization=="supervised",], aes(x=interaction(randomization, estimator, lex.order=TRUE), y=ATE_hat)) + 
  geom_boxplot() + 
  stat_summary(fun.y = "mean", geom = "point", colour = "blue", shape = 15) + #, size = 2
  geom_hline(aes(yintercept=ATE),linetype="dashed") +
  labs(x="Randomization", y = "ATE") +
  scale_x_discrete(labels=c("full.DR" = "full (balanced)", "imbalanced.DR" = "full (imbalanced)",
                            "supervised.IPW" = "supervised (IPW)","supervised.DR" = "supervised (DR)")) +
  theme_classic()+
  theme(axis.text.x = element_text(size=12, color="black"))+
  # theme(axis.text.x = element_text(size=14),
  #       axis.text.y = element_text(size=14),
  #       axis.title.x = element_text( size=16),
  #       axis.title.y = element_text(size=16)) +
  ylim(c(0.02,0.08))

ggsave(paste0("ATE_boxplot_",NO_FOLDS*NO_EXPERIMENT_ITER, ".pdf"), plot = ATE_plot, device = "pdf", path=results_path,
              width = 14, height = 8, units = "cm",
              dpi = 400)

##### CATE #####
# Set model parameters
# CF
MTRY= ceil(sqrt(ncol(X)))
NUM.TREES=500
CI.GROUP.SIZE=1
MIN.NODE.SIZE = 20
SAMPLE.FRACTION = 0.5


pred_CATE <- foreach(exp_idx = 1:nrow(exp_grid), .inorder=TRUE, .packages=c("grf", "rpart", "nnet"), .export=c("predict.tlearner"))%dopar%{
  # Select experiment from grid
  randomization          <- exp_grid$randomization[exp_idx]
  uplift_model_error_idx <- exp_grid$uplift_model_error_idx[exp_idx]
  testset                <- exp_grid$testset[exp_idx]
  i                      <- exp_grid$iteration[exp_idx]
  
  # Test set observations
  test_idx <- test_idx_list[[testset ]]
  X_train <- X[-test_idx, ]
  X_test <- X[test_idx, ]

  # Simulate existing uplift model
  tau_hat <- TAU + uplift_model_error_list[[uplift_model_error_idx ]]
  
  # Get score quantiles from training set for mapping 
  score_quantiles_ <- get_score_quantiles_(tau_hat[-test_idx], groups=18)
  supervised_prop_score <- map_propensity_quantiles(model_score = tau_hat[-test_idx], score_breaks = score_quantiles_, target_ratio = 0.5)
  # supervised_prop_score <- map_propensity_logistic(model_score = tau_hat[-test_idx], 
  #                                                  min_score = quantile(tau_hat[-test_idx], probs=0.05), 
  #                                                  max_score = quantile(tau_hat[-test_idx], probs=0.95))
  
  if(randomization=="none") prop_score <- rep(0, nrow(X_train))
  if(randomization=="all") prop_score <- rep(1, nrow(X_train))
  if(randomization=="full") prop_score <- rep(0.5, nrow(X_train))
  if(randomization=="imbalanced") prop_score <- rep(IMBALANCED_EXP_RATIO, nrow(X_train))
  if(randomization=="supervised") prop_score <- supervised_prop_score
  
  # Assign treatment
  experiment <- do_experiment(Y1=Y[-test_idx, "Y1"], Y0=Y[-test_idx, "Y0"], prop_score = prop_score)
  y <- experiment$outcome
  w <- experiment$treatment
  
  ### CATE ###
  pred <- list()
  
  ## ATE baseline
  ate_hat <- calc_ATE(y = y, w = w, prop_score = prop_score)
  pred[["ATE"]] <- rep(ate_hat, times = nrow(X_test) ) 
  
  ## Oracle regression
  pred[["true"]] <- TAU[test_idx]
  
  oracle_reg <- lm(.y~., data = data.frame(".y"=TAU[-test_idx], X_train), weights = ifelse(w==1, prop_score, 1-prop_score))
  pred[['oracle_regression']] <- unname(predict(oracle_reg, X_test))
  
  ## T Logit
  t_logit <- T_Logit(X = X_train, y = y, W = w, prop_score = prop_score)
  pred[['t-logit']] <- unname(predict(t_logit, X_test, type='response'))
  
  ## MOM Logit (see Knaus 2018)
  #mom_logit <- MOM_Logit(X = X_train, y = y, W = w, prop_score = prop_score)
  #pred[['mom-logit']] <- unname(predict(mom_logit, X_test))
  
  ## T Tree
  # t_rpart <- T_rpart(X = X_train, y = y, W = w, prop_score = prop_score, method="class")
  # pred[['t-tree']] <- unname(predict(t_rpart, X_test, type='prob'))
  
  ## T Neural Network
  #t_nnet <- T_NNet(X = X_train, y = y, W = w, prop_score = prop_score, trace=TRUE)
  #pred[['t-nnet']] <- unname(predict(t_nnet, X_test, type='raw'))
  
  #Causal Forest
  cf <- grf::causal_forest(X=X_train,Y=y, W=w,
                           W.hat = prop_score,
                           honesty=TRUE,
                           ci.group.size = CI.GROUP.SIZE,
                           num.trees=NUM.TREES, min.node.size=MIN.NODE.SIZE,
                           mtry=MTRY, sample.fraction = SAMPLE.FRACTION)

  pred[['CF']] <- predict(cf, X_test)[,1]
  
  # Collect results
  return(pred)
}

saveRDS(pred_CATE, "../CATE_predictions_20190826.rds")

# 
CATE_perf <- foreach(exp_idx=1:nrow(exp_grid), .combine=rbind, .packages = "data.table", .inorder=TRUE, .export=)%dopar%{
    # Select experiment from grid
    randomization          <- exp_grid$randomization[exp_idx]
    uplift_model_error_idx <- exp_grid$uplift_model_error_idx[exp_idx]
    testset                <- exp_grid$testset[exp_idx]
    i                      <- exp_grid$iteration[exp_idx]
    
    test_idx <- test_idx_list[[testset]]
    Y_test <-   Y[test_idx, ]
    
    tau_scores <- pred_CATE[[exp_idx]]
    
    # Results
    raw <- list()
    raw[["mae"]] <- sapply(tau_scores, function(x) mean(abs(x - TAU[test_idx])))
    raw[["qini"]] <- sapply(tau_scores, function(x) qini_score(scores = c(x, x), Y = c(Y_test[,"Y1"], Y_test[,"Y0"]), W = c(rep(1, length(Y_test[,"Y1"])), rep(0, length(Y_test[,"Y0"])))))
    
    for(value in VALUE){
    targeting_decisions <- sapply(tau_scores, function(x) targeting_policy(tau_hat=x, offer_cost = CONTACT_COST, customer_value = value), simplify = FALSE)
    raw[[paste0("targeting_ratio_", value)]] <- sapply(targeting_decisions, mean)
    raw[[paste0("profit_cutoff_", value)]]          <- sapply(targeting_decisions, function(x) sum(x * (Y_test[,"Y1"] * value - CONTACT_COST) + (1-x)*(Y_test[,"Y0"] * value)))
    
    targeting_decisions_top_decile          <- sapply(tau_scores, function(x) top_decile_targeting_policy(x), simplify=FALSE)
    raw[[paste0("profit_decile_", value)]]     <- sapply(targeting_decisions_top_decile, function(x) sum(x * (Y_test[,"Y1"] * value - CONTACT_COST) + (1-x)*(Y_test[,"Y0"] * value)))
    }
    
    temp <- data.frame(raw)
    temp <- cbind(temp, "model"=row.names(temp))
    temp <- melt(temp, id.vars = "model", variable.name=c("metric"))
    temp <- cbind(
      temp, 
      data.frame(list("uplift_model_error"=UPLIFT_MODEL_ERROR[uplift_model_error_idx], "testset"=testset, "randomization_iteration" = i, "randomization"=randomization), 
                 stringsAsFactors = FALSE, row.names = NULL) 
      )
    return(temp)
}

CATE_summary <- CATE_perf %>% group_by(randomization, uplift_model_error, metric, model) %>% summarize("value"=mean(value))
model_names <- c("true","ATE",  "oracle_regression", "t-logit", "t-tree", "mom-logit", "CF", "t-nnet")
CATE_summary$model <- factor(CATE_summary$model, levels = model_names, labels = model_names)

CATE_summary <- data.frame(CATE_summary)
row.names(CATE_summary) <- NULL

saveRDS(CATE_summary, "../CATE_summary_20190826.rds")

### Statistical Performance Table
stat_table <- dcast(CATE_summary[ CATE_summary$metric %in% c("mae","qini") & CATE_summary$model %in% c("ATE","t-logit","CF") & (CATE_summary$randomization!="supervised" | CATE_summary$uplift_model_error==DEFAULT_TARGETING_ERROR),], 
                    randomization ~ model+ metric, row.names=FALSE)
row.names(stat_table) <- NULL                    

sink(paste0(results_path, "CATE_stat_table.txt"))
print(stat_table, digits=4)
cat("\n\n\n")
print( 
  xtable(stat_table, digits=4),
  include.rownames=FALSE)
sink()

### Profit Performance Table
## Targeting ratio as a sanity check
stat_table <- dcast(CATE_summary[ grepl(pattern = "targeting_ratio", CATE_summary$metric) & CATE_summary$model %in% c("ATE","t-logit","CF") & (CATE_summary$randomization!="supervised" | CATE_summary$uplift_model_error==DEFAULT_TARGETING_ERROR),], 
                    metric ~ model + randomization, row.names=FALSE)
row.names(stat_table) <- NULL     

sink(paste0(results_path, "CATE_targeting-ratio_table.txt"))
print(stat_table)
cat("\n\n\n")
print( 
  xtable(stat_table, digits=4),
  include.rownames=FALSE)
sink()

## Profit based on profit decision
stat_table <- dcast(CATE_summary[ grepl(pattern = "profit_cutoff", CATE_summary$metric) & CATE_summary$model %in% c("t-logit","CF") & (CATE_summary$randomization!="supervised" | CATE_summary$uplift_model_error==DEFAULT_TARGETING_ERROR),], 
                    metric ~ model + randomization, row.names=FALSE)
row.names(stat_table) <- NULL   
stat_table$metric <-  str_remove(stat_table$metric, "profit_cutoff_")

# Per customer profit
stat_table[, 2:ncol(stat_table)] <- stat_table[, 2:ncol(stat_table)]/mean(sapply(test_idx_list, length))

sink(paste0(results_path, "CATE_profit_table.txt"))
print(stat_table, digits=2)
cat("\n\n\n")
print( 
  xtable(stat_table, digits=2),
  include.rownames=FALSE)
sink()

## Profit at fixed targeting rate
stat_table <- dcast(CATE_summary[ grepl(pattern = "profit_decile", CATE_summary$metric) & CATE_summary$model %in% c("t-logit","CF") & CATE_summary$uplift_model_error==DEFAULT_TARGETING_ERROR,], 
                    metric ~ model + randomization, row.names=FALSE)
row.names(stat_table) <- NULL 


### Targeting model performance robustness check ####

## Experiment profit
error_exp_profit <- result_exp %>% group_by(randomization, uplift_model_error, campaign_margin) %>% summarize("campaign_profit"=mean(campaign_profit))
row.names(error_exp_profit) <- NULL
error_exp_profit$randomization <- factor(error_exp_profit$randomization, randomization_values, randomization_labels)

error_exp_profit_table <- melt(error_exp_profit, id.vars = c("randomization", "uplift_model_error","campaign_margin"))
error_exp_profit_table$value <- error_exp_profit_table$value / mean(sapply(test_idx_list, function(test_idx) nrow(X)-length(test_idx)))

error_exp_profit_table <- dcast( error_exp_profit_table[error_exp_profit$randomization=="Supervised",], campaign_margin ~ uplift_model_error)

print(xtable(error_exp_profit_table), include.rownames=FALSE)

## Model profit
# Statistical performance
error_model_profit <- dcast(CATE_summary[ CATE_summary$metric %in% c("mae","qini") & CATE_summary$model %in% c("ATE","t-logit","CF") & CATE_summary$randomization=="supervised",], 
                     uplift_model_error ~ model+ metric, row.names=FALSE)

print(xtable(error_model_profit, digits=3), include.rownames=FALSE)

# CATE Profit
error_model_profit <- dcast(CATE_summary[ grepl(pattern = "profit_cutoff", CATE_summary$metric)  & CATE_summary$model %in% c("t-logit","CF") & CATE_summary$randomization=="supervised",],  #
                            metric ~ model + uplift_model_error, row.names=FALSE)

error_model_profit[, 2:ncol(error_model_profit)] <- error_model_profit[, 2:ncol(error_model_profit)]/mean(sapply(test_idx_list, length))

print(xtable(error_model_profit, digits=2), include.rownames=FALSE)

