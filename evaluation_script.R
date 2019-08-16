#### Packages ####
#install.packages("pacman")
library(pacman)
pacman::p_load("ggplot2","reshape2","drtmle","grf","foreach","uplift","data.table","ModelMetrics", "plyr", "parallel", "doParallel", "SuperLearner", "tidyverse")

source("load_bankmarketing.R")
#source("data_generating_process.R")
source("supervised_randomization.R")

source("t_logit.R")
source("Qini.R")
source("costs.R")

calc_ATE <- function(y, w, prop_score){
  #if(is.null(prop_score)){
  #  prop_score = rep(0.5, length(y))
  #}
  return( (sum(y*w/prop_score) - sum(y*(1-w)/(1-prop_score)) ) /length(y) )
}


#### Real data ####
source("load_bankmarketing.R")
path = "../data/bank-additional-full.csv"
data <- load_bankmarketing(path)


# TODO: Neural Network Process best?
tau_model <- function(X, hidden_layer, ATE){
  X <- as.matrix(X)
  no_var <- ncol(X)
  W1 <- matrix(rnorm(no_var*hidden_layer, 0, 1),    nrow=no_var, ncol=hidden_layer)
  W2 <- matrix(rnorm(hidden_layer,        0, 1), nrow=hidden_layer, ncol=1)

  h <- X %*% W1
  h <- dlogis(X %*% W1)
  
  o <- c(h %*% W2)
  o <- o / (25*sd(o))
  o <- o - mean(o) + ATE
  return(o)
}

set.seed(123456789)
X <- data[,!(colnames(data)=="y")]
X_tau <- X[,c(1:24, 43:48)]
TAU <- tau_model(X_tau, hidden_layer=ncol(X_tau), ATE=0.05)
X <- X[,!(colnames(X) %in% c("age", "maritalmarried", "maritalsingle"))]

summary(TAU)


R <- rbinom(nrow(data), size = 1, prob = abs(TAU))
Y <- cbind(Y0=data$y, Y1=data$y)

for(i in 1:nrow(Y)){
  if(R[i] == 1){
      if(TAU[i]>0) Y[i, ] <- c(0,1)
      if(TAU[i]<0) Y[i, ] <- c(1,0)
  }
}

# Sanity check: ATE should match target ATE
mean(Y[,2]) - mean(Y[,1])


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


#### Register cluster for parallel processing ####
cl <- makeCluster(3)
registerDoParallel(cl)
RNGkind("L'Ecuyer-CMRG")                                                          
clusterSetRNGStream(cl,iseed = 123456789)

#### SET PARAMETERS ####
# Cross-Validation split
NO_FOLDS = 4
NO_EXPERIMENT_ITER = 50
IMBALANCED_EXP_RATIO = 0.75

set.seed(1234)
test_idx_list <- caret::createFolds(Y[,2] - Y[,1], k = NO_FOLDS, list = TRUE, returnTrain = FALSE)

# Model error
UPLIFT_MODEL_ERROR <- c(0, 0.005, 0.01, 0.02, 0.05, 0.1)
uplift_model_error_list <- foreach(e=UPLIFT_MODEL_ERROR)%do% rnorm(length(TAU), 0, e)

# Profit setting
CONTACT_COST = 2 # Contact costs
OFFER_COST = 0 # Price reduction

VALUE = 0.5*c(20,40,60,80,100,120,140)

exp_grid = data.frame(expand.grid(list(
  "uplift_model_error_idx" = 1:length(uplift_model_error_list),
  "testset"=1:length(test_idx_list),
  "iteration"=1:NO_EXPERIMENT_ITER,
  "randomization" = c("full", "imbalanced", "supervised")
)))
# The uplift model error is only relevant for supervised randomization
exp_grid <- exp_grid[ !((exp_grid$uplift_model_error_idx != 1) & (exp_grid$randomization != "supervised")), ]

#### Cost Summary  ####
exp_grid_ABtest <- rbind(exp_grid, 
                         data.frame(expand.grid(list(
                           "uplift_model_error_idx" = 1,
                           "testset"=1:length(test_idx_list),
                           "iteration"=1:NO_EXPERIMENT_ITER,
                           "randomization" = c("none", "all")
                         )))
)

# For each quality level of existing model
result_ATE <- foreach(exp_idx = 1:nrow(exp_grid_ABtest), .inorder=TRUE, .combine=rbind)%dopar%{
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
          
          if(randomization=="none") prop_score <- 0
          if(randomization=="all") prop_score <- 1
          if(randomization=="full") prop_score <- 0.5
          if(randomization=="imbalanced") prop_score <- IMBALANCED_EXP_RATIO
          if(randomization=="supervised") prop_score <- supervised_prop_score
          
          # Assign treatment
          experiment <- do_experiment(Y1=Y[-test_idx, "Y1"], Y0=Y[-test_idx, "Y0"], prop_score = prop_score)
          y <- experiment$outcome
          w <- experiment$treatment
          
          # Collect results
          data.frame(list(
            "uplift_model_error"=UPLIFT_MODEL_ERROR[uplift_model_error_idx], "testset"=testset, "randomization_iteration" = i, "randomization"=randomization,
            "treatment_ratio"=mean(w), "outcome_ratio"=mean(y)
            ), stringsAsFactors = FALSE) 
}

library(tidyverse)
exp_summary <- result_ATE %>% group_by(randomization, uplift_model_error) %>% summarize("outcome_ratio"=mean(outcome_ratio), "treatment_ratio"=mean(treatment_ratio))
row.names(exp_summary) <- NULL

print(exp_summary)

fwrite(data.table(exp_summary,keep.rownames = F), 
       "../experiment_summary_real_20190815.txt")


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

result_ATE$estimator <- factor(result_ATE$estimator, levels = c("IPW","DR"), labels=c("IPW", "DR"))

ATE_plot <- ggplot(data = result_ATE[result_ATE$estimator=="DR",], aes(x=interaction(randomization, estimator, lex.order=FALSE), y=ATE_hat)) + 
  geom_boxplot() + 
  stat_summary(fun.y = "mean", geom = "point", colour = "blue", shape = 15) + #, size = 2
  geom_hline(aes(yintercept=mean(TAU)),linetype="dashed") +
  labs(x="Randomization", y = "ATE") +
  scale_x_discrete(labels=c("full.DR" = "full (balanced)", "imbalanced.DR" = "full (imbalanced)",
                            "supervised.DR" = "supervised")) +
  theme_classic()+
  theme(axis.text.x = element_text(size=12, color="black"))+
  # theme(axis.text.x = element_text(size=14),
  #       axis.text.y = element_text(size=14),
  #       axis.title.x = element_text( size=16),
  #       axis.title.y = element_text(size=16)) +
  ylim(c(0.02,0.08))

ggsave(paste0("ATE_boxplot_",NO_FOLDS*NO_EXPERIMENT_ITER, ".pdf"), plot = ATE_plot, device = "pdf", path="../ICIS19_revision",
              width = 14, height = 8, units = "cm",
              dpi = 400)

##### CATE #####
# Set model parameters
# CF
MTRY=7
NUM.TREES=500
CI.GROUP.SIZE=1
MIN.NODE.SIZE = 20
SAMPLE.FRACTION = 0.5


pred_CATE <- foreach(exp_idx = 1:nrow(exp_grid), .inorder=TRUE, .packages=c("grf", "rpart"), .export=c("predict.tlearner"))%dopar%{
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
  mom_logit <- MOM_Logit(X = X_train, y = y, W = w, prop_score = prop_score)
  pred[['mom-logit']] <- unname(predict(mom_logit, X_test))
  
  ## T Tree
  t_rpart <- T_rpart(X = X_train, y = y, W = w, prop_score = prop_score, method="class")
  pred[['t-tree']] <- unname(predict(t_rpart, X_test, type='prob'))
  
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

# 
CATE_perf <- foreach(exp_idx=1:nrow(exp_grid), .combine=rbind, .packages = "data.table", .inorder=TRUE)%dopar%{
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
    raw[["qini"]] <- sapply(tau_scores, function(x) qini_score(scores = c(x, x), Y = c(Y_test[,"Y1"], Y_test[,"Y0"]), W = c(rep(1, length(Y_test[,"Y1"])), rep(0, length(Y_test[,"Y0"]))), p_treatment = prop_score))
    
    for(value in VALUE){
    targeting_decisions <- sapply(tau_scores, function(x) targeting_policy(tau_hat=x, offer_cost = CONTACT_COST, customer_value = value), simplify = FALSE)
    raw[[paste0("targeting_ratio_", value)]] <- sapply(targeting_decisions, mean)
    raw[[paste0("profit_", value)]]          <- sapply(targeting_decisions, function(x) sum(x * (Y_test[,"Y1"] * value - CONTACT_COST) + (1-x)*(Y_test[,"Y0"] * value)))
    
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
row.names(CATE_summary) <- NULL

CATE_summary <- data.frame(CATE_summary)

print(CATE_summary[CATE_summary$metric=="mae",])
print(CATE_summary[CATE_summary$metric=="qini",])
print(CATE_summary[CATE_summary$metric=="profit_30",])
print(CATE_summary[CATE_summary$metric=="profit_decile_30",])
