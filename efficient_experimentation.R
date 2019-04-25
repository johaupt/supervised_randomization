#### Packages ####
#library(pacman)
#pacman::p_load("ggplot2","reshape","cowplot","car","drtmle","SuperLearner","stargazer","grf","foreach","uplift")

if(!require("ggplot2")) install.packages("ggplot2"); library("ggplot2")
if(!require("reshape")) install.packages("reshape"); library("reshape")
if(!require("cowplot")) install.packages("cowplot"); library("cowplot")
if(!require("car")) install.packages("car"); library("car")
if(!require("drtmle")) install.packages("drtmle"); library("drtmle")
if(!require("SuperLearner")) install.packages("SuperLearner"); library("SuperLearner")
if(!require("stargazer")) install.packages("stargazer"); library("stargazer")
if(!require("grf")) install.packages("grf"); library("grf")
if(!require("foreach")) install.packages("foreach"); library("foreach")
if(!require("uplift")) install.packages("uplift"); library("uplift")

source("data_generating_process.R")

N_VAR=20
N_CUSTOMER=1e5
#RATIO_SAMPLE=0.05

# Full set of customers
X_Customers <- sapply(1:N_VAR, function(x)rnorm(N_CUSTOMER, mean = 0, sd=1))
X <- data.frame(X_Customers)

expCtrl <- expControl(n_var = N_VAR, mode = "classification", beta_zero = -1.75,  # >0 indicates more than 50% purchasers
                      tau_zero =   0.75, # >0 indicates positive treatment effect)
                      DGP="nonlinear")


#### Run experiments ####
exp <- list()
exp$none <- do_experiment(X, expControl = expCtrl, prop_score = 0)
exp$all <- do_experiment(X, expControl = expCtrl, prop_score = 1)
# Balanced random treatment
exp$balanced <- do_experiment(X, expControl = expCtrl, prop_score = 0.5)
# Conservative random treatment
exp$imbalanced <- do_experiment(X, expControl = expCtrl, prop_score = 0.66)
# baseline
mean(exp$none$y)

#### Examplary basic churn reponse model
# This can be any model that maps some business value to a treatment probability
# For churn, the churn probability is an obvious candidate, but could be rescaled
response_model <- glm(y~., cbind(X, y=exp$none$y), family = binomial(link="logit"))
churn_pred <- predict(response_model, X, type = "response")
ModelMetrics::auc(exp$none$y, churn_pred)
## Individual biased random treatment (based on churn probability)
# Platt scaling
churn_pred <- pmin(pmax(churn_pred, 0.05), 0.95)
# platt_scaler <- glm(y~prob, family=binomial(link='logit'), 
#                     weights = ifelse(exp$none$y==1, 1/(mean(exp$none$y==1)), 1),
#                     data = data.frame(cbind('y'=exp$none$y,'prob'=churn_pred))
#                     )
# treat_prob <- predict(platt_scaler, newdata=data.frame(prob=churn_pred), type='response')
#treat_prob <- pmin(pmax(treat_prob, 0.1), 0.9)
treat_prob <- churn_pred
exp$individual <- do_experiment(X, expControl = expCtrl, prop_score = treat_prob)

# # Second response model (probit)
# response_model2 <- glm(y~., cbind(X, y=exp$none$y), family = binomial(link="probit"))
# churn_pred2 <- predict(response_model2, X, type = "response")
# ModelMetrics::auc(exp$none$y, churn_pred)
# churn_pred2 <- pmin(pmax(churn_pred, 0.05), 0.95)
# treat_prob2 <- churn_pred2
# exp$individual2 <- do_experiment(X, expControl = expCtrl, prop_score = treat_prob2)
# 
# # Third response model (cauchit)
# response_model3 <- glm(y~., cbind(X, y=exp$none$y), family = binomial(link="cauchit"))
# churn_pred3 <- predict(response_model3, X, type = "response")
# ModelMetrics::auc(exp$none$y, churn_pred)
# churn_pred3 <- pmin(pmax(churn_pred, 0.05), 0.95)
# treat_prob3 <- churn_pred3
# exp$individual3 <- do_experiment(X, expControl = expCtrl, prop_score = treat_prob3)
# 
# # Fourth response model (cloglog; complementary log-log)
# response_model4 <- glm(y~., cbind(X, y=exp$none$y), family = binomial(link="cloglog"))
# churn_pred4 <- predict(response_model4, X, type = "response")
# ModelMetrics::auc(exp$none$y, churn_pred)
# churn_pred4 <- pmin(pmax(churn_pred, 0.05), 0.95)
# treat_prob4 <- churn_pred4
# exp$individual4 <- do_experiment(X, expControl = expCtrl, prop_score = treat_prob4)

# # Output influence of covariates on target variable for response models in publication-ready quality
# stargazer(response_model, type="text") # logit model
# stargazer(response_model2, type="text") # probit model
# stargazer(response_model3, type="text") # cauchit model
# stargazer(response_model4, type="text") # cloglog model

### Experiment outcomes ####
EXPERIMENT_SIZE = 1e5 # Number of people in experiment
CONTACT_COST = 2 # Contact costs

CLV_matrix <- c(10, 50, 100, 200, 500, 1000, 5000, 10000, 50000)
cost_all <- matrix(NA,nrow=length(CLV_matrix),ncol=6)
colnames(cost_all) <- c("CLV","none","all","balanced","imbalanced","individual")

source("costs.R")
cost <- catalogue_profit

for(j in 1:length(CLV_matrix)) {
  
  CLV = CLV_matrix[j] # Customer lifetime value
  
  
  COST_TREATMENT_VAR = 1/20*CLV # Price reduction
  COST_CHURN = CLV # Foregone profit

  # Ratio of treated
  sapply(exp, function(x)mean(x$g))
  # Churn rate
  sapply(exp, function(x)mean(x$y))
  # Expected outcome per customer (max. 0, higher is better)
  sapply(exp[c("none","all","balanced","imbalanced","individual")], 
         function(A) cost(A$y, A$g, 
                                COST_TREATMENT_FIX, COST_TREATMENT_VAR, COST_CHURN))
  
  # Churn costs per scenario and unit of observation
  sapply(exp[c("none","all","balanced","imbalanced","individual")],
         function(B) cost(B$y, B$g, 
                                COST_TREATMENT_FIX, COST_TREATMENT_VAR, COST_CHURN) / EXPERIMENT_SIZE)
  
  churn_cost_scenario <- as.vector(sapply(exp[c("none","all","balanced","imbalanced","individual")],
                                          function(B) cost(B$y, B$g, 
                                                                 COST_TREATMENT_FIX, COST_TREATMENT_VAR, COST_CHURN) / EXPERIMENT_SIZE))
  
  
  cost_all[j,1] <- CLV_matrix[j]
  cost_all[j,c(2:ncol(cost_all))] <- churn_cost_scenario 
  
}

cost_all


#### Create experiments ####
# Fixed response model
response_model <- glm(y~., cbind(X, y=exp$none$y), family = binomial(link="logit"))
churn_pred <- predict(response_model, X, type = "response")
treat_prob <- pmin(pmax(churn_pred, 0.05), 0.95)

balanced <- list()
imbalanced <- list()
individual <- list()
# Repeat sampling n times



set.seed(123)

NO_EXPERIMENT_ITER = 50

for(i in 1:NO_EXPERIMENT_ITER){
  balanced[[i]] <- do_experiment(X, expControl = expCtrl, prop_score = 0.5)
  imbalanced[[i]] <- do_experiment(X, expControl = expCtrl, prop_score = 0.25)
  individual[[i]] <- do_experiment(X, expControl = expCtrl, prop_score = treat_prob)
}



### ATE Estimation ####
calc_ATE <- function(y, g, prop_score){
  #if(is.null(prop_score)){
  #  prop_score = rep(0.5, length(y))
  #}
  return( (sum(y*g/prop_score) - sum(y*(1-g)/(1-prop_score)) ) /length(y) )
}

ATE <- data.frame()

for(i in 1:NO_EXPERIMENT_ITER){
  ATE[i,"balanced"] <- calc_ATE(balanced[[i]]$y, balanced[[i]]$g, prop_score = 0.5)
  ATE[i,"balanced_dr"] <- ci(drtmle(Y=balanced[[i]]$y,A=balanced[[i]]$g,W=X,a_0 = c(1,0),
                                    family=binomial(),
                                    stratify=TRUE,
                                    SL_Q = c("SL.glm"),
                                    SL_g = c("SL.glm"),
                                    SL_Qr = "SL.glm",
                                    SL_gr = "SL.glm", maxIter = 1),contrast=c(1,-1))$drtmle[1]
  
  ATE[i,"imbalanced"] <- calc_ATE(imbalanced[[i]]$y, imbalanced[[i]]$g, prop_score = 0.25)
  
  ATE[i,"individual"] <- calc_ATE(individual[[i]]$y, individual[[i]]$g, individual[[i]]$prop_score)
  
  ATE[i,"individual_dr"] <- ci(drtmle(Y=individual[[i]]$y,A=individual[[i]]$g,W=X,a_0 = c(1,0),
                                      family=binomial(),
                                      stratify=TRUE,
                                      SL_Q = c("SL.glm"),
                                      SL_g = c("SL.glm"),
                                      SL_Qr = "SL.glm",
                                      SL_gr = "SL.glm", maxIter = 1),contrast=c(1,-1))$drtmle[1]
}



# True ATE
mean(exp$all$y) - mean(exp$none$y)
# Estimated ATE
ATE_hat <- apply(ATE,2,mean)
ATE_hat

ate_box <- melt(ATE)
 ggplot(data = ate_box, aes(x=variable, y=value)) + 
  geom_boxplot() + 
  stat_summary(fun.y = "mean", geom = "point", colour = "blue", shape = 15, size = 2) +
  geom_hline(aes(yintercept=mean(exp$all$y) - mean(exp$none$y)),colour="red") +
  labs(x="Experiment design", y = "ATE") +
  scale_x_discrete(labels=c("balanced" = "balanced", "imbalanced" = "imbalanced",
                            "individual" = "supervised (IPW)", "individual_dr"= "supervised (DR)"))
                                          
  
                                       
                                          
                                          
# T-Test for mean Difference (H0)
t.test(ATE$balanced,ATE$imbalanced) # iterations: 200, H0: diff in mean = 0 accepeted
t.test(ATE$balanced,ATE$individual) # iterations: 200, H0: diff in mean = 0 accepeted
t.test(ATE$individual,ATE$individual_dr) # iterations: 200, H0: diff in mean = 0 accepeted

# Test for Normal-distribution (H0)
#shapiro.test(ATE$balanced) # H0: not sig. different from normal distributon
#shapiro.test(ATE$individual_dr)

#plot(density(ATE$balanced))
#plot(density(ATE$individual))
#plot(density(ATE$individual_dr))

# Test for equal variance (H0)
# Group samples 
lev_sample <- c(ATE$balanced, ATE$individual,ATE$individual_dr)


lev_group <- as.factor(c(rep("b", length(ATE$balanced)), rep("ind", length(ATE$individual)),rep("ind_dr", length(ATE$individual_dr))))
leveneTest(lev_sample,lev_group) # H0: Homogeneity of Variance 

                                        
                                          
                                          
#### CATE Estimation####


# Check the MAE and especially Qini for the cost of
# training CATE models on biased and corrected instead 
# of 1:1 randomized data

# New metric can be specified here
performance_CATE <- function(tau_score, y_true=NULL, w=NULL, tau_true=NULL){
  #res <- data.frame("metric" = c("MAE","Qini"), "value"=NA)
  #if(!is.null(tau_true))             res[res$metric=="MAE","value"]  <- mean(abs(tau_true - tau_score))
  #if(!is.null(y_true) & !is.null(w)) res[res$metric=="Qini","value"] <- qini_score(tau_score, y_true, w)

  res <- list()
  if(!is.null(tau_true))  res[["MAE"]] <- mean(abs(tau_true - tau_score))
  if(!is.null(y_true) & !is.null(w)) res[["Qini"]] <- qini_score(tau_score, y_true, w)
  
  return(res)
}

MTRY=5
NUM.TREES=100


## Two-model approach
# TODO: Finish evaluation of two-model approach
source("t_logit.R")
source("Qini.R")



exp_list <- list("balanced"=balanced, "individual"=individual)

perf_CATE <- foreach(exp_setting = exp_list, .combine='list', .inorder=TRUE, .multicombine=TRUE) %:%
               foreach(exp=exp_setting, .combine = "rbind", .multicombine = TRUE) %dopar% {
                 #res <- data.frame()
                 # Predictions from each model are saved in a list 'pred'
                 pred <- list()
                 
                 # Split data
                 idx_train <- sample(length(exp$y), size = floor(length(exp$y)*0.75), replace=FALSE)
                 
                 ## ATE baseline
                 ate_hat <- calc_ATE(exp$y, exp$g, exp$prop_score)
                 
                 ## Train t-learner
                 t_logit <- T_Logit(X[idx_train,], exp$y[idx_train], exp$g[idx_train],
                                    exp$prop_score[idx_train])
          
                # Save predictions on test set in list pred
                 pred[['t_logit']] <- unname(predict(t_logit, X[-idx_train,]))
                 
                 #res <- rbind(res, 
                #       cbind("model"="t_logit", performance_CATE(tau_score=predict(t_logit, X[-idx_train,]),
                #                                                                                   y_true=exp$y[-idx_train], w = exp$g[-idx_train], tau_true = exp$tau[-idx_train]))
                 #      )
                 
                 
                 ## Train causal forest
                 cf <- grf::causal_forest(X=X[idx_train,],Y=exp$y[idx_train], W=exp$g[idx_train],
                                          W.hat = exp$prop_score[idx_train],
                                          honesty=TRUE,
                                          num.trees=NUM.TREES, min.node.size=100, 
                                          mtry=MTRY, sample.fraction = 0.2,
                                          seed=123)
               
                 pred[['CF']] <- predict(cf, X[-idx_train,])[,1]
                 
                 # Calculate performance for each model
                 res <- sapply(pred, simplify=FALSE, function(PRED) performance_CATE(tau_score= PRED,
                                                      y_true=exp$y[-idx_train], w = exp$g[-idx_train], tau_true = exp$tau[-idx_train]))
                 res <- melt(res)
                 colnames(res)[match(c("L1","L2"), names(res))] <- c("model","metric")
                 
                 return(res)
               }

# Calculate mean performance
names(perf_CATE) <- names(exp_list)

melt(sapply(perf_CATE, function(x) tapply(x$value, INDEX=list(x$model, x$metric), mean), simplify=FALSE))

sapply(perf_CATE, function(x) tapply(x$value, INDEX=list(x$model, x$metric), sd), simplify=FALSE)

t.test(perf_CATE$t_logit$balanced,perf_CATE$t_logit_DR$balanced)
t.test(perf_CATE$t_logit_DR$balanced,perf_CATE$CF$balanced)

t.test(perf_CATE$t_logit$individual,perf_CATE$t_logit_DR$individual)
t.test(perf_CATE$t_logit_DR$individual,perf_CATE$CF$individual)

