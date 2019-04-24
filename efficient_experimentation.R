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
EXPERIMENT_SIZE = 100000 # Number of people in experiment
COST_TREATMENT_FIX = 1 # Contact costs

CLV_matrix <- c(10, 50, 100, 200, 500, 1000, 5000, 10000, 50000)
cost_all <- matrix(NA,nrow=length(CLV_matrix),ncol=9)
colnames(cost_all) <- c("CLV","none","all","balanced","imbalanced","individual", "individual2", "individual3", "individual4")

source("costs.R")

for(j in 1:length(CLV_matrix)) {
  
  CLV = CLV_matrix[j] # Customer lifetime value
  
  
  COST_TREATMENT_VAR = 1/20*CLV # Price reduction
  COST_CHURN = CLV # Foregone profit

  # Ratio of treated
  sapply(exp, function(x)mean(x$g))
  # Churn rate
  sapply(exp, function(x)mean(x$y))
  # Expected outcome per customer (max. 0, higher is better)
  sapply(exp[c("none","all","balanced","imbalanced","individual", "individual2", "individual3", "individual4")], 
         function(A) churn_cost(A$y, A$g, 
                                COST_TREATMENT_FIX, COST_TREATMENT_VAR, COST_CHURN))
  
  # Churn costs per scenario and unit of observation
  sapply(exp[c("none","all","balanced","imbalanced","individual", "individual2", "individual3", "individual4")],
         function(B) churn_cost(B$y, B$g, 
                                COST_TREATMENT_FIX, COST_TREATMENT_VAR, COST_CHURN) / EXPERIMENT_SIZE)
  
  churn_cost_scenario <- as.vector(sapply(exp[c("none","all","balanced","imbalanced","individual", "individual2", "individual3", "individual4")],
                                          function(B) churn_cost(B$y, B$g, 
                                                                 COST_TREATMENT_FIX, COST_TREATMENT_VAR, COST_CHURN) / EXPERIMENT_SIZE))
  
  
  cost_all[j,1] <- CLV_matrix[j]
  cost_all[j,c(2:ncol(cost_all))] <- churn_cost_scenario 
  
}

cost_all

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
 ATE[i,"balanced"] <- calc_ATE(balanced[[i]]$y, balanced[[i]]$g, prop_score = 0.5)
 imbalanced[[i]] <- do_experiment(X, expControl = expCtrl, prop_score = 0.25)
 ATE[i,"imbalanced"] <- calc_ATE(imbalanced[[i]]$y, imbalanced[[i]]$g, prop_score = 0.25)
 individual[[i]] <- do_experiment(X, expControl = expCtrl, prop_score = treat_prob)
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
ate_plot <- ggplot(data = ate_box, aes(x=variable, y=value)) + 
  geom_boxplot() + 
  stat_summary(fun.y = "mean", geom = "point", colour = "blue", shape = 15, size = 2) +
  geom_hline(aes(yintercept=mean(exp$all$y) - mean(exp$none$y)),colour="red") +
  labs(x="Experiment design", y = "ATE") +
  scale_x_discrete(labels=c("balanced" = "balanced", "imbalanced" = "imbalanced",
                            "individual" = "supervised (IPW)", "individual_dr"= "supervised (DR)"))
                                          
  plot_grid(ate_plot) + theme(plot.background = element_rect(color = "black", size=0.7))
                                       
                                          
                                          
# T-Test for mean Difference (H0)
t.test(ATE$balanced,ATE$imbalanced) # iterations: 200, H0: diff in mean = 0 accepeted
t.test(ATE$balanced,ATE$individual) # iterations: 200, H0: diff in mean = 0 accepeted
t.test(ATE$individual,ATE$individual_dr) # iterations: 200, H0: diff in mean = 0 accepeted

# Test for Normal-distribution (H0)
shapiro.test(ATE$balanced) # H0: not sig. different from normal distributon
shapiro.test(ATE$individual_dr)

plot(density(ATE$balanced))
plot(density(ATE$individual))
plot(density(ATE$individual_dr))

# Test for equal variance (H0)
# Group samples 
lev_sample <- c(ATE$balanced, ATE$individual,ATE$individual_dr)


lev_group <- as.factor(c(rep("b", length(ATE$balanced)), rep("ind", length(ATE$individual)),rep("ind_dr", length(ATE$individual_dr))))
leveneTest(lev_sample,lev_group) # H0: Homogeneity of Variance -> rejected -> ind_dr has smaller variance 

                                        
                                          
                                          
#### CATE Estimation####


# Check the MAE and especially Qini for the cost of
# training CATE models on biased and corrected instead 
# of 1:1 randomized data

## Two-model approach
# TODO: Finish evaluation of two-model approach
source("t_logit.R")
source("Qini.R")
perf_CATE <- list()

perf_CATE[["t_logit"]][["balanced"]] <- foreach(exp=balanced[1:20],
                                              .combine = "rbind", .multicombine = TRUE) %do% {
                                               t_logit <- T_Logit(X,exp$y, exp$g,
                                                                  exp$prop_score)
                                               tau_hat <- predict(t_logit, X)
                                               MAE <- mean(abs(exp$tau - tau_hat))
                                               Qini <- qini_score(tau_hat, exp$y, exp$g)
                                               c("MAE" = MAE, "Qini" = Qini)
                                   }

perf_CATE[["t_logit"]][["individual"]] <- foreach(exp=individual[1:20],
                                                .combine = "rbind", .multicombine = TRUE) %do% {
                                                  t_logit <- T_Logit(X,exp$y, exp$g,
                                                                     exp$prop_score)
                                                  tau_hat <- predict(t_logit, X)
                                                  MAE <- mean(abs(exp$tau - tau_hat))
                                                  Qini <- qini_score(tau_hat, exp$y, exp$g)
                                                  c("MAE" = MAE, "Qini" = Qini)
                                                }

# Double Robust for CATE 
source("t_logit_DR.R")
perf_CATE[["t_logit_DR"]][["balanced"]] <- foreach(exp=balanced[1:20],
                                                   .combine = "rbind", .multicombine = TRUE) %do% {
                                                     t_logit_dr <- T_Logit_DR(X,exp$y, exp$g,
                                                                           exp$prop_score)
                                                     tau_hat <- predict(t_logit_dr, X)
                                                     MAE <- mean(abs(exp$tau - tau_hat))
                                                     Qini <- qini_score(tau_hat, exp$y, exp$g)
                                                     c("MAE" = MAE, "Qini" = Qini)
                                                   }



perf_CATE[["t_logit_DR"]][["individual"]] <- foreach(exp=individual[1:20],
                                                     .combine = "rbind", .multicombine = TRUE) %do% {
                                                       t_logit_dr <- T_Logit_DR(X,exp$y, exp$g,
                                                                             exp$prop_score)
                                                       tau_hat <- predict(t_logit_dr, X)
                                                       MAE <- mean(abs(exp$tau - tau_hat))
                                                       Qini <- qini_score(tau_hat, exp$y, exp$g)
                                                       c("MAE" = MAE, "Qini" = Qini)
                                                     }


## Causal Forest ####

# Build causal trees based on balanced and efficient experiments
# and compare mean absolute error and Qini score
# TODO: The ATE is again a competitive predictor. That's 
#       weird, find out why! 

MTRY=5
NUM.TREES=1000

perf_CATE[["CF"]][["balanced"]] <- foreach(exp=balanced[1:20],
                                 .combine="rbind",.multicombine=TRUE) %do%{
cf <- grf::causal_forest(X=X,Y=exp$y, W=exp$g,
                 W.hat = exp$prop_score,
                 honesty=TRUE,
                 num.trees=NUM.TREES, min.node.size=10, 
                 mtry=MTRY, sample.fraction = 0.2,
                 seed=123)
tau_hat <- predict(cf, X)
c("MAE"=mean(abs(exp$tau - tau_hat)[,1]),
  "Qini"=qini_score(tau_hat[,1], exp$y, exp$g))
}

perf_CATE[["CF"]][["individual"]] <- foreach(exp=individual[1:20],
                                   .combine="rbind",.multicombine=TRUE) %do%{
  cf <- grf::causal_forest(X=X,Y=exp$y, W=exp$g,
                           W.hat = exp$prop_score,
                           honesty=TRUE,
                           num.trees=NUM.TREES, min.node.size=10, 
                           mtry=MTRY, sample.fraction=0.5,
                           seed=123)
  tau_hat <- predict(cf, X)
  c("MAE"=mean(abs(exp$tau - tau_hat)[,1]),
    "Qini"=qini_score(tau_hat[,1], exp$y, exp$g))
}

plot.default(tau_hat[,1], exp$tau[,1])

mean(abs(exp$tau - ATE_hat["balanced"])[,1])
res <- lapply(perf_CATE, lapply, function(x) colMeans(x))
data.frame(res)

