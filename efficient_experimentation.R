#### Packages ####
#library(pacman)

pacman::p_load("ggplot2","reshape","reshape2","cowplot","car","drtmle","grf","foreach","uplift","data.table")

source("data_generating_process.R")

N_VAR=20
N_CUSTOMER=100000
EXPERIMENT_SIZE=N_CUSTOMER
#RATIO_SAMPLE=0.05

expCtrl <- expControl(n_var = N_VAR, mode = "classification", beta_zero = -3,  # >0 indicates more than 50% purchasers
                      tau_zero =   0.8, # >0 indicates positive treatment effect)
                      DGP="nonlinear")

#### Examplary basic churn reponse model
# This can be any model that maps some business value to a treatment probability
# For churn, the churn probability is an obvious candidate, but could be rescaled
X_train  <- make_customers(N_CUSTOMER, N_VAR)
X_test  <- make_customers(N_CUSTOMER, N_VAR)
response_model <- glm(y~., family = binomial(link="logit"),
                      cbind(X_train, y=do_experiment(X_train, expControl = expCtrl, prop_score = 0)$y) 
)
response_pred <- predict(response_model, X_test, type = "response")
ModelMetrics::auc(do_experiment(X_test, expControl = expCtrl, prop_score = 0)$y, 
                  response_pred)

## Individual biased random treatment (based on propensity)
map_propensity <- function(model_score, target_ratio, groups=9){
# Platt scaling
# platt_scaler <- glm(y~prob, family=binomial(link='logit'), 
#                     weights = ifelse(exp$none$y==1, 1/(mean(exp$none$y==1)), 1),
#                     data = data.frame(cbind('y'=exp$none$y,'prob'=churn_pred))
#                     )
# treat_prob <- predict(platt_scaler, newdata=data.frame(prob=churn_pred), type='response')
  # Cut into groups based on the score quantiles
  model_score <- cut(model_score,breaks=quantile(model_score,seq(0,1,1/groups)),labels=FALSE, include.lowest = TRUE)-1

  # Adjust to expected target ratio by shifting min or max
  if(target_ratio<=0.5){
    new_max <- 2*target_ratio - 0.05
    model_score <- 0.05 + model_score* (new_max-0.05)/(groups-1)
  }
  if(target_ratio>0.5){
    new_min <- 2*target_ratio - 0.95
    model_score <- new_min + model_score* (0.95-new_min)/(groups-1)
  }

 return(model_score)
}


#### Create experiments ####
X <- list()
none <- list()
all <- list()
balanced <- list()
imbalanced <- list()
individual <- list()
test <- list()

set.seed(123)
# Repeat sampling n times  
NO_EXPERIMENT_ITER = 100
IMBALANCED_EXP_RATIO = 0.66

for(i in 1:NO_EXPERIMENT_ITER){
  X_iter <- make_customers(EXPERIMENT_SIZE, 20)
  X[[i]] <- X_iter
  balanced[[i]] <- do_experiment(X_iter, expControl = expCtrl, prop_score = 0.5)
  imbalanced[[i]] <- do_experiment(X_iter, expControl = expCtrl, prop_score = IMBALANCED_EXP_RATIO)
  individual[[i]] <- do_experiment(X_iter, expControl = expCtrl, 
                                   prop_score = map_propensity(predict(response_model, X_iter, type="response"), 
                                                               target_ratio=IMBALANCED_EXP_RATIO, groups=20))
  none[[i]] <- do_experiment(X_iter, expControl = expCtrl, prop_score = 0)
  all[[i]] <- do_experiment(X_iter, expControl = expCtrl, prop_score = 1)
}


#### Experiment Summary statistics ####
exp_summary <- foreach(exp=list(none, all,balanced, imbalanced, individual),.combine="rbind", .inorder = TRUE)%do%{
  temp <- sapply(exp, function(x)   c("treatment_ratio" = mean(x$g), "response_ratio" = mean(x$y)) )
  c(rowMeans(temp), "response_ratio_sd"=sd(temp["response_ratio",]) )
}
row.names(exp_summary) <- c("none","all","balanced","imbalanced","supervised")
round(t(exp_summary),3)

#### Experiment outcomes ####
CONTACT_COST = 2 # Contact costs
OFFER_COST = 0 # Price reduction

VALUE_matrix <- c(20, 40, 60, 80, 100, 120, 140)
profit_all <- matrix(NA,nrow=length(VALUE_matrix),ncol=6)
colnames(profit_all) <- c("basket","none","all","balanced","imbalanced","individual")

source("costs.R")
profit <- catalogue_profit

### TODO: Should not only relate to one experiment, but repeat for different X
exp = list("balanced"=balanced[[1]], "imbalanced"=imbalanced[[1]], 
           "individual"=individual[[1]], "none"=none[[1]],"all"=all[[1]])

for(j in 1:length(VALUE_matrix)) {
  
  VALUE = VALUE_matrix[j] # Basket value
  
  
  # Ratio of treated
  sapply(exp, function(x)mean(x$g))
  # Churn rate
  sapply(exp, function(x)mean(x$y))
  # Expected outcome per customer (max. 0, higher is better)
  sapply(exp[c("none","all","balanced","imbalanced","individual")], 
         function(A) profit(A$y, A$g, 
                            contact_cost = CONTACT_COST, offer_cost = OFFER_COST, value=VALUE))
  
  # Churn costs per scenario and unit of observation
  sapply(exp[c("none","all","balanced","imbalanced","individual")],
         function(B) profit(B$y, B$g, 
                            contact_cost = CONTACT_COST, offer_cost = OFFER_COST, value=VALUE) / EXPERIMENT_SIZE)
  
  profit_scenario <- as.vector(sapply(exp[c("none","all","balanced","imbalanced","individual")],
                                      function(B) profit(B$y, B$g, 
                                                         contact_cost = CONTACT_COST, offer_cost = OFFER_COST, value=VALUE) / EXPERIMENT_SIZE))
  
  
  profit_all[j,1] <- VALUE_matrix[j]
  profit_all[j,c(2:ncol(profit_all))] <- profit_scenario 
  
}

df_profit <- as.data.frame(profit_all)
df_profit <- df_profit[,-c(2)] # remove scenario "none"

# Robustness check: Profit per customer for different basket values
d <- melt(df_profit, id.vars="basket")
colnames(d) <- c("Basket_value","Scenario","Profit")

Fig_BasketProfit <- ggplot(d, aes(Basket_value,Profit, col=Scenario)) + 
  geom_line() +
  #ggtitle("Average profit per customer for different basket values") +
  ylab("Profit per customer") + xlab("Basket value")

#df_profit

write.csv(df_profit, file = "TableBasketValuesScenarios2.csv")

ggsave("../ICIS19/FigureBasketProfit.pdf")

df_profit <- cbind(df_profit, df_profit[,c("imbalanced","individual")]-df_profit[,"balanced"])
print(df_profit)

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

  ATE[i,"imbalanced"] <- calc_ATE(imbalanced[[i]]$y, imbalanced[[i]]$g, prop_score = IMBALANCED_EXP_RATIO)
  
  ATE[i,"individual"] <- calc_ATE(individual[[i]]$y, individual[[i]]$g, individual[[i]]$prop_score)
  
  ATE[i,"balanced_dr"] <- ci(drtmle(Y=balanced[[i]]$y,A=balanced[[i]]$g,W=X[[i]],a_0 = c(1,0),
                                    family=binomial(),
                                    stratify=TRUE,
                                    SL_Q = c("SL.glm"),
                                    SL_g = c("SL.glm"),
                                    SL_Qr = "SL.glm",
                                    SL_gr = "SL.glm", maxIter = 1),contrast=c(1,-1))$drtmle[1]

  ATE[i,"individual_dr"] <- ci(drtmle(Y=individual[[i]]$y,A=individual[[i]]$g,W=X[[i]],a_0 = c(1,0),
                                      family=binomial(),
                                      stratify=TRUE,
                                      SL_Q = c("SL.glm"),
                                      SL_g = c("SL.glm"),
                                      SL_Qr = "SL.glm",
                                      SL_gr = "SL.glm", maxIter = 1),contrast=c(1,-1))$drtmle[1]
  # True ATE
  ATE[i,"true_ATE"] <- mean(all[[i]]$y) - mean(none[[i]]$y)
}





# Estimated ATE
ATE_hat <- apply(ATE,2,mean)
ATE_hat

ate_box <- melt(ATE)

ggplot(data = ate_box, aes(x=variable, y=value)) + 
  geom_boxplot() + 
  stat_summary(fun.y = "mean", geom = "point", colour = "blue", shape = 15, size = 2) +
  geom_hline(aes(yintercept=ATE_hat[6]),colour="red") +
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

### Model evaluation ####
source("t_logit.R")
source("Qini.R")

library("parallel")
library("doParallel")
cl <- makeCluster( max(1,detectCores()-1))
registerDoParallel(cl)
RNGkind("L'Ecuyer-CMRG")  
clusterSetRNGStream(cl,iseed = 1234567)

N_ITER=100

## Model specifics
# CF
MTRY=4
NUM.TREES=500
CI.GROUP.SIZE=1
MIN.NODE.SIZE = 100
SAMPLE.FRACTION = 0.5

# Return predictions and true values for each of N sets of customers X
model_library <- foreach(i=1:N_ITER, .combine="list", .multicombine=TRUE, .export = c("predict.tlearner"), .errorhandling = "pass", .packages = "data.table")%dopar%{
  X <- make_customers(N_CUSTOMER, N_VAR)
  idx_split <- ceiling(N_CUSTOMER/2)
  X_train <- X[1:idx_split,]
  X_test <- X[(idx_split+1):nrow(X),]
  
  exp <- list(
  "balanced" = do_experiment(X_train, expControl = expCtrl, prop_score = 0.5),
  "imbalanced" = do_experiment(X_train, expControl = expCtrl, prop_score = IMBALANCED_EXP_RATIO),
  "individual" = do_experiment(X_train, expControl = expCtrl, prop_score = map_propensity(predict(response_model, X_train, type="response"), target_ratio=IMBALANCED_EXP_RATIO))
  )
  
  # Predictions from each model are saved in a list 'pred'
  pred <- list()
  
  for(rand_scheme in names(exp)){
  ## ATE baseline
  ate_hat <- calc_ATE(exp[[rand_scheme]]$y, exp[[rand_scheme]]$g, exp[[rand_scheme]]$prop_score)
  pred[[paste0("ATE_",rand_scheme)]] <- rep(ate_hat, times = nrow(X_test) ) 
  
  ## Train t-learner
  t_logit <- T_Logit(X_train, exp[[rand_scheme]]$y, exp[[rand_scheme]]$g,
                     exp[[rand_scheme]]$prop_score)
  
  # Save predictions on test set in list pred
  pred[[paste0('t-logit_',rand_scheme)]] <- unname(predict(t_logit, X_test))
  
  ## Train causal forest
  cf <- grf::causal_forest(X=X_train,Y=exp[[rand_scheme]]$y, W=exp[[rand_scheme]]$g,
                           W.hat = exp[[rand_scheme]]$prop_score,
                           honesty=TRUE,
                           ci.group.size = CI.GROUP.SIZE,
                           num.trees=NUM.TREES, min.node.size=MIN.NODE.SIZE, 
                           mtry=MTRY, sample.fraction = SAMPLE.FRACTION)
  
  pred[[paste0('CF_', rand_scheme)]] <- predict(cf, X_test)[,1]
  }
  pred <- as.data.table(pred)
  
  test <- do_experiment(X_test, expControl = expCtrl, prop_score = 0.5)
  
  return(list("pred"=pred, "true"=test))
}

#saveRDS(model_library, "../ICIS19/model_library_20190428.rds")
model_library <- readRDS("../ICIS19/model_library_20190428.rds")

#### Calculate performance for each model ####
# New metric can be specified here
performance_CATE <- function(tau_score, y_true=NULL, w=NULL, prop_score=NULL, tau_true=NULL){
  res <- list()
  if(!is.null(tau_true))  res[["MAE"]] <- mean(abs(tau_true - tau_score))
  if(!is.null(y_true) & !is.null(w)) res[["Qini"]] <- qini_score(scores = tau_score, Y = y_true, W = w, p_treatment = prop_score)

  for(basket_value in c(20,40,60,80,100,120,140)){
    res[[paste0("profit_",basket_value)]] <- catalogue_profit(y=y_true, contact_cost = 2, offer_cost = 0, value = basket_value,
                                 g=targeting_policy(
                                   tau_hat = tau_score, offer_cost = 2, customer_value = basket_value)
                                 )
  }
  
  return(as.data.table(res))
}

# Calculate performance per iteration
res <- lapply(model_library, function(ITER) lapply(ITER$pred, function(PRED){
  performance_CATE(tau_score= PRED,
                   y_true=ITER$true$y, w = ITER$true$g, prop_score = ITER$true$prop_score, tau_true = ITER$true$tau)
} ))

# Combine models into one data table per iteration
res <- lapply(res, function(x) rbindlist(x,idcol="model"))#, keep.rownames="model"))
# Combine iterations into one data table
res <- rbindlist(res, idcol="iter")
# Restructure to long format
res <- melt(res, id.vars = c("model","iter"), variable.name = "metric")
# Calculate summary statistics over iterations
res <- res[,.("perf_mean" = mean(value), "perf_sd" = sd(value)),by=.(model, metric)]
setorder(res, metric,model)

# CATE model performance for statistical performance measures/KPIs
res_KPI <- dcast(res[!grepl("profit",metric),], metric~model, value.var = "perf_mean")

# CATE model performance for profit
res_profit <- dcast(res[grepl("profit",metric),], metric~model, value.var = "perf_mean")
res_profit[,"metric"] <- gsub(res_profit[,"metric"],pattern = "profit_",replacement = "")
res_profit[,-1] <- res_profit[,-1]/(0.5*N_CUSTOMER)
for(model in c("ATE","CF","logit")){
res_profit[,grepl(model, colnames(res_profit))] <- res_profit[,grepl(model, colnames(res_profit))] - 
                                                   unlist(res_profit[grepl(paste0(model,"_balanced"), colnames(res_profit))])
}

# Final tables
res_KPI
res_profit

# t.test(perf_CATE$t_logit$balanced,perf_CATE$t_logit_DR$balanced)
# t.test(perf_CATE$t_logit_DR$balanced,perf_CATE$CF$balanced)
# 
# t.test(perf_CATE$t_logit$individual,perf_CATE$t_logit_DR$individual)
# t.test(perf_CATE$t_logit_DR$individual,perf_CATE$CF$individual)

