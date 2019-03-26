#### Experiment Functions ####
# Define the experiment to be able to run it several times with 
# the same coefficients
expControl <- function(n_var, mode="regression", tau_zero=NULL, beta_zero=NULL){
  ###
  # n_var (int): Number of variables
  # model (str): 'regression' or binary 'classification'
  # tau_zero (float): Homogeneous/baseline treatment effect
  # beta_zero (float): Response constant
  ###
  beta_zero = ifelse(is.null(beta_zero),0,beta_zero)
  beta = rnorm(n_var, mean = 0, sd = 0.2)
  beta_x2 = rnorm(n_var, mean = 0, sd = 0.2)
  beta_tau = rnorm(n_var, mean = 0, sd = 0.2)
  if(is.null(tau_zero)){
      tau_zero = rnorm(1, mean=-0.1, sd=0.01)
  }
  
  return(list("beta_zero"=beta_zero, "beta"=beta, "beta_x2"=beta_x2, 
              "beta_tau"=beta_tau, "tau_zero" = tau_zero,
              "mode"=mode))
  
}

# Create data given the data generating process defined in 
# the experiment control object
do_experiment <- function(X, expControl, g=NULL, prop_score=NULL, X_out=FALSE){
  ###
  # X (array):
  #    Matrix of observations
  # g (vector of int): 
  #    Binary treatment assignment (0: control, 1:treatment) 
  # prop_score (single float or vector of floats): 
  #    Propensity score for random treatment. If g is not NULL it will override prop_score
  #
  ###
  n_obs = nrow(X)
  beta_zero = expControl$beta_zero
  beta = expControl$beta
  beta_x2 = expControl$beta_x2
  beta_tau = expControl$beta_tau
  tau_zero = expControl$tau_zero
  mode = expControl$mode
  
  logit <- function(x) 1/(1+exp(-x))
  
  if(!"matrix" %in% class(X)){
    X <- as.matrix(X)
  }
  
  if(is.null(g)){
    if(is.null(prop_score)){
        g = rep(0, times=nrow(X))
    }else{
        g = rbinom(nrow(X), 1, prop_score)
        if(length(prop_score)==1){
          prop_score <- rep(prop_score, times=nrow(X))
        }
    }
  }
  
if(mode %in% c('regression','classification')){
  
  if(mode == "regression"){
      tau = tau_zero + X%*%beta_tau + rnorm(n_obs, 0, 0.01)
      y = beta_zero + X%*%beta + g*tau + rnorm(n_obs, 0, 0.5)
  }
  
  if(mode == "classification"){
    tau = tau_zero + X%*%beta_tau + rnorm(n_obs, 0, 0.01)
    
    y = beta_zero + X%*%beta + rnorm(n_obs, 0, 0.5)
    
    y0 = logit(y)
    y1 = logit(y+tau)
    tau = y1-y0
    
    y = logit(y+g*tau)
    y = as.numeric(y>=0.5)
  }
}else{
  stop("Only mode 'classification' or 'regression' currently implemented")
}
  
  if(X_out==FALSE){
    X = NULL
  }
  
  return(list("X"=X, "y"=y, "tau"=tau, "g"=g, "prop_score"=prop_score))
  
}

N_VAR=20
N_CUSTOMER=100000
RATIO_SAMPLE=0.2

# Full set of customers
X_Customers <- sapply(1:N_VAR, function(x)rnorm(N_CUSTOMER, mean = 0, sd=1))
X_Customers <- data.frame(X_Customers)
# Select subset for experimentation
# TODO: In practice, we could withhold treatment for selected customers only, does that work?
#       Would produce only the in-sample treatment effect (-> reject inference problem)
X <- X_Customers[sample(1:N_CUSTOMER, size = N_CUSTOMER*RATIO_SAMPLE, replace = FALSE),]
expCtrl <- expControl(n_var = N_VAR, mode = "classification", 
                      beta_zero=-0.6, tau_zero = -0.05)

#### Run experiments ####
exp <- list()
exp$none <- do_experiment(X, expControl = expCtrl, prop_score = 0)
exp$all <- do_experiment(X, expControl = expCtrl, prop_score = 1)
# Balanced random treatment
exp$balanced <- do_experiment(X, expControl = expCtrl, prop_score = 0.5)
# Conservative random treatment
exp$imbalanced <- do_experiment(X, expControl = expCtrl, prop_score = 0.25)
# Churn baseline
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

### Experiment outcomes ####
EXPERIMENT_SIZE = 150000 # Number of people in experiment
COST_TREATMENT_FIX = 0 # Contact costs
CLV = 50 # Customer lifetime value
COST_TREATMENT_VAR = 1/20*CLV # Price reduction
COST_CHURN = CLV # Foregone profit
# Expected churn costs per customer
churn_cost <- function(n_customer, y, g, cost_treatment_fix, 
                       cost_treatment_var, cost_churn){
    (mean(g)*               -cost_treatment_fix +
    mean(g)*(1-mean(y))*   -cost_treatment_var +
    mean(y)*               -cost_churn)*n_customer
}

# Ratio of treated
sapply(exp, function(x)mean(x$g))
# Churn rate
sapply(exp, function(x)mean(x$y))
# Expected outcome per customer (max. 0, higher is better)
sapply(exp[c("none","all","balanced","imbalanced","individual")], 
       function(A) churn_cost(EXPERIMENT_SIZE, A$y, A$g, 
                     COST_TREATMENT_FIX, COST_TREATMENT_VAR, COST_CHURN))

### ATE Estimation ####
calc_ATE <- function(y, g, prop_score){
  #if(is.null(prop_score)){
  #  prop_score = rep(0.5, length(y))
  #}
  return( (sum(y*g/prop_score) - sum(y*(1-g)/(1-prop_score)) ) /length(y) )
}

# Fixed response model
response_model <- glm(y~., cbind(X, y=exp$none$y), family = binomial(link="logit"))
churn_pred <- predict(response_model, X, type = "response")
treat_prob <- pmin(pmax(churn_pred, 0.05), 0.95)

ATE <- data.frame()
balanced <- list()
imbalanced <- list()
individual <- list()
# Repeat sampling n times
for(i in 1:200){
 balanced[[i]] <- do_experiment(X, expControl = expCtrl, prop_score = 0.5)
 ATE[i,"balanced"] <- calc_ATE(balanced[[i]]$y, balanced[[i]]$g)
 imbalanced[[i]] <- do_experiment(X, expControl = expCtrl, prop_score = 0.25)
 ATE[i,"imbalanced"] <- calc_ATE(imbalanced[[i]]$y, imbalanced[[i]]$g, prop_score = 0.25)
 individual[[i]] <- do_experiment(X, expControl = expCtrl, prop_score = treat_prob)
 ATE[i,"individual"] <- calc_ATE(individual[[i]]$y, individual[[i]]$g, individual[[i]]$prop_score)
 }

# True ATE
mean(exp$all$y) - mean(exp$none$y)
# Estimated ATE
ATE_hat <- apply(ATE,2,mean)
ATE_hat
boxplot(ATE)
abline(h=mean(exp$all$y) - mean(exp$none$y),col="red")

#### CATE Estimation####
library(foreach)
library(uplift)
# Check the MAE and especially Qini for the cost of
# training CATE models on biased and corrected instead 
# of 1:1 randomized data

## Two-model approach
# TODO: Finish evaluation of two-model approach
source("t_logit.R")
perf_CATE <- list()

perf_CATE[["t_logit"]][["balanced"]] <- foreach(exp=balanced[1:10],
                                   .combine = "rbind", .multicombine = TRUE) %do% {
   t_logit <- T_Logit(X,exp$y, exp$g,
                     exp$prop_score)
   tau_hat <- predict(t_logit, X)
   MAE <- mean(abs(exp$tau - tau_hat))
   perf <- uplift::performance(tau_hat, rep(0, times = length(tau_hat)), exp$y, exp$g,
                               direction =1) # Negative treatment effect
   Qini <- qini(perf, direction = 1, plotit = FALSE)$Qini
   c("MAE" = MAE, "Qini" = Qini)
                                   }

perf_CATE[["t_logit"]][["individual"]] <- foreach(exp=individual[1:10],
                                                .combine = "rbind", .multicombine = TRUE) %do% {
                                                  t_logit <- T_Logit(X,exp$y, exp$g,
                                                                     exp$prop_score)
                                                  tau_hat <- predict(t_logit, X)
                                                  MAE <- mean(abs(exp$tau - tau_hat))
                                                  perf <- uplift::performance(tau_hat, rep(0, times = length(tau_hat)), exp$y, exp$g,
                                                                              direction =1) # Negative treatment effect
                                                  Qini <- qini(perf, direction = 1, plotit = FALSE)$Qini
                                                  c("MAE" = MAE, "Qini" = Qini)
                                                }


## Causal Forest ####
library(causalTree)
library(grf)
# Build causal trees based on balanced and efficient experiments
# and compare mean absolute error and Qini score
# TODO: The ATE is again a competitive predictor. That's 
#       weird, find out why! 

MTRY=5
NUM.TREES=1000

perf_CATE[["CF"]][["balanced"]] <- foreach(exp=balanced[1:10],
                                 .combine="rbind",.multicombine=TRUE) %do%{
cf <- grf::causal_forest(X=X,Y=exp$y, W=exp$g,
                 W.hat = exp$prop_score,
                 honesty=TRUE,
                 num.trees=NUM.TREES, min.node.size=10, 
                 mtry=MTRY, sample.fraction = 0.5,
                 seed=123)
tau_hat <- predict(cf, X)
MAE <- mean(abs(exp$tau - tau_hat)[,1])
perf <- uplift::performance(tau_hat[,1],rep(0,times=length(tau_hat[,1])), exp$y, exp$g,
                            direction=1) # Negative treatment effect
Qini <- qini(perf, direction=1, plotit=FALSE)$Qini
c("MAE"=MAE,"Qini"=Qini)
}

perf_CATE[["CF"]][["individual"]] <- foreach(exp=individual[1:10],
                                   .combine="rbind",.multicombine=TRUE) %do%{
  cf <- grf::causal_forest(X=X,Y=exp$y, W=exp$g,
                           W.hat = exp$prop_score,
                           honesty=TRUE,
                           num.trees=NUM.TREES, min.node.size=10, 
                           mtry=MTRY, sample.fraction=0.5,
                           seed=123)
  tau_hat <- predict(cf, X)
  MAE <- mean(abs(exp$tau - tau_hat)[,1])
  perf <- uplift::performance(tau_hat[,1],rep(0,times=length(tau_hat[,1])), exp$y, exp$g,
                      direction=1) # Negative treatment effect
  Qini <- qini(perf, direction=1, plotit=FALSE)$Qini
  c("MAE"=MAE,"Qini"=Qini)
}

plot.default(tau_hat[,1], exp$tau[,1])

mean(abs(exp$tau - ATE_hat["balanced"])[,1])
lapply(perf_CATE, lapply, function(x) colMeans(x))



