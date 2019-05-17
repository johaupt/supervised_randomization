#### Packages ####
#install.packages("pacman")
library(pacman)
pacman::p_load("ggplot2","reshape2","drtmle","grf","foreach","uplift","data.table","ModelMetrics", "plyr", "parallel", "doParallel", "SuperLearner")

source("t_logit.R")
source("Qini.R")
source("data_generating_process.R")
source("supervised_randomization.R")
source("costs.R")

#### Simulation ####
Y = NULL
W = NULL

N_VAR=20
N_CUSTOMER=100000
EXPERIMENT_SIZE=N_CUSTOMER
#RATIO_SAMPLE=0.05

set.seed(1234)
expCtrl <- expControl(n_var = N_VAR, mode = "classification", beta_zero = -3,  # >0 indicates more than 50% purchasers
                      tau_zero =   0.425, # >0 indicates positive treatment effect)
                      DGP="nonlinear")

X  <- make_customers(N_CUSTOMER, N_VAR)
exp_temp <- do_experiment(X, expControl = expCtrl, prop_score = 0.5)
Y <- exp_temp$y
W <- exp_temp$w
TAU <- exp_temp$tau

####
expCtrl <- NULL
source("load_uplift19.R")
path = "../data/explore.csv"
data <- load_uplift19(path)

# Sample the amount of control observations from the indices of the treatment group
idx_balanced <- c( which(data[["W"]] == 1)[sample(sum(data[["W"]] == 1), size=sum(data[["W"]] == 0), replace=FALSE)], 
                   which(data[["W"]] == 0))

X <- data[["X"]][idx_balanced,]
Y <- data[["Y"]][idx_balanced]
W <- data[["W"]][idx_balanced]
X_VALUE <- data[["value"]][idx_balanced]


#### Examplary basic churn reponse model
# This can be any model that maps some business value to a treatment probability
# For churn, the churn probability is an obvious candidate, but could be rescaled
set.seed(123)
# Remove training data for targeting model from dataset
idx_targeting_model <- sample(nrow(X), size=10000)
X_targeting <- X[idx_targeting_model,]
Y_targeting <- Y[idx_targeting_model]
W_targeting <- W[idx_targeting_model]
X <- X[-idx_targeting_model,]
Y <- Y[-idx_targeting_model]
W <- W[-idx_targeting_model]

# Train targeting model: We assume that a company would have an existing model
#targeting_model <- glm(y~., family = binomial(link="logit"), cbind(X, "y"=Y) )
targeting_model <- T_Logit(X_targeting, Y_targeting, W_targeting, prop_score = rep(0.5, length(idx_targeting_model)))
#targeting_model <- T_Logit_DR(X_targeting, Y_targeting, W_targeting)
# targeting_model <- causalTree::causalTree(y~., cbind("y"=Y[idx_targeting_model], X[idx_targeting_model,]), treatment = W[idx_targeting_model], 
#                                           split.Rule='CT', split.Honest = TRUE, cv.option="CT", split.Bucket=FALSE,
#                                           cp=0,minsize=40)
 # targeting_model <- grf::causal_forest(X= X_targeting, Y = Y_targeting, W = W_targeting,
 #                         #W.hat = 0.5,
 #                         honesty=TRUE,
 #                         ci.group.size = 1,
 #                         num.trees=2000, min.node.size=12,
 #                         mtry=5, sample.fraction = 0.5)

## Training data performance
ModelMetrics::auc(actual=Y_targeting[W_targeting==0], predicted=predict(targeting_model$model0, X_targeting, type = "response")[W_targeting==0])
ModelMetrics::auc(actual=Y_targeting[W_targeting==1], predicted=predict(targeting_model$model0, X_targeting, type = "response")[W_targeting==1])

pred_targeting <- predict(targeting_model, X_targeting, type = "response")
#pred_targeting <- predict(targeting_model, X_targeting)[,1]
cor(pred_targeting[W_targeting==1], Y_targeting[W_targeting==1],method = "spearman")
qini_score(scores = pred_targeting, Y_targeting, W_targeting, plotit=TRUE)
transformed_outcome_loss(pred_targeting, Y_targeting, W_targeting, p_treatment = 0.5)

# Test data performance
ModelMetrics::auc(actual=Y[W==0], predicted=predict(targeting_model$model0, X[W==0,], type = "response"))
ModelMetrics::auc(actual=Y[W==1], predicted=predict(targeting_model$model1, X[W==1,], type = "response"))
pred <- predict(targeting_model, X, type = "response")
#pred <- predict(targeting_model, X)[,1]
cor(pred[W==1], Y[W==1], method = "spearman")
qini_score(scores = pred, Y, W, plotit=TRUE)
transformed_outcome_loss(pred, Y, W, p_treatment = 0.5)


#### Create experiments ####
#X <- list()
none <- list()
all <- list()
balanced <- list()
imbalanced <- list()
individual <- list()
test <- list()

set.seed(1234)
# Repeat sampling n times  
NO_EXPERIMENT_ITER = 100
IMBALANCED_EXP_RATIO = 0.75

# ICIS submission:
# Repeat experiment NO_EXPERIMENT_ITER times. For each iteration:
# - Fix the customer data X_iter and save it to list X
# - Run experiment for each randomization procedure, save experiment in respective list
# NOW:
# Fix customer data X. For each iteration:
# - Run experiment for each randomization procedure, save experiment in respective list
X_iter <- X
for(i in 1:NO_EXPERIMENT_ITER){
  #X_iter <- make_customers(EXPERIMENT_SIZE, 20)
  #X[[i]] <- X_iter
  balanced[[i]] <- do_experiment(X_iter, expControl = expCtrl, y=Y, w=W, prop_score = 0.5)
  imbalanced[[i]] <- do_experiment(X_iter, expControl = expCtrl, y=Y, w=W, prop_score = IMBALANCED_EXP_RATIO)
  individual[[i]] <- do_experiment(X_iter, expControl = expCtrl, y=Y, w=W,
                                   prop_score = map_propensity(predict(targeting_model, X_iter), 
                                                               target_ratio=IMBALANCED_EXP_RATIO, groups=20))
  none[[i]] <- do_experiment(X_iter, expControl = expCtrl, y=Y, w=W, prop_score = 0)
  all[[i]] <- do_experiment(X_iter, expControl = expCtrl, y=Y, w=W, prop_score = 1)
}

#### Experiment Summary statistics ####
exp_summary <- foreach(exp=list(none, all,balanced, imbalanced, individual),.combine="rbind", .inorder = TRUE)%do%{
  temp <- sapply(exp, function(x)   c("treatment_ratio" = mean(x$w), "response_ratio" = mean(x$y)) )
  c(apply(temp,1,median), "treatment_ratio_sd"=sd(temp["treatment_ratio",]) ,"response_ratio_sd"=sd(temp["response_ratio",]) )
}
row.names(exp_summary) <- c("none","all","balanced","imbalanced","supervised")
round(t(exp_summary),3)

fwrite(data.table(round(t(exp_summary),3),keep.rownames = T), 
       "../experiment_summary_real_20190515.txt")

#### Experiment outcomes ####
CONTACT_COST = 1 # Contact costs
OFFER_COST = 0 # Price reduction

VALUE = 0.5*c(20,40,60,80,100,120,140)

VALUE_matrix <- rep(VALUE, NO_EXPERIMENT_ITER) # 100 iterations x 7 values
profit_all <- matrix(NA,nrow=length(VALUE_matrix),ncol=5+1) # create empty profit matrix
colnames(profit_all) <- c("basket","balanced","imbalanced","individual","none","all")

source("costs.R")
profit <- catalogue_profit

for(i in 1:NO_EXPERIMENT_ITER){

exp = list("balanced"=balanced[[i]], "imbalanced"=imbalanced[[i]], 
           "individual"=individual[[i]], "none"=none[[i]],"all"=all[[i]])

for(j in 1:length(VALUE)) {
  
  VALUE_iter = VALUE_matrix[j] # Basket value
  
  profit_scenario <- as.vector(sapply(
    exp,
    #exp[c("none","all","balanced","imbalanced","individual")],
      function(B) profit(B$y, B$w, contact_cost = CONTACT_COST, offer_cost = OFFER_COST, value=VALUE_iter) / length(B$y)))
  
profit_all

  profit_all[(i-1)*length(VALUE)+j,1] <- VALUE_matrix[j]
  profit_all[(i-1)*length(VALUE)+j,c(2:ncol(profit_all))] <- profit_scenario 
  
}
}
profit_scenario
  
df_profit <- as.data.frame(profit_all)

#df_profit <- df_profit[,-c(5)] # remove procedure "none"
#df_profit <- as.data.table(df_profit)
#df_profit[, colMeans(.SD), by=basket]
df_exhaustive <- df_profit # store all profit values from x experiments for the different profit margins

df_profit <- aggregate(.~basket, FUN=mean, data=df_profit)

# Plot
#d <- melt(df_profit, id.vars="basket")
#d$variable <- revalue(d$variable, c("individual"="supervised"))
#colnames(d) <- c("BasketMargin","Procedure","Profit")

#Fig_ProfitBasketMargin <- ggplot(d, aes(BasketMargin,Profit, col=Procedure)) + 
#  geom_line() +
#  ylab("Profit per customer") + xlab("Basket margin")

write.csv(df_profit, file = "RawProfit.csv")

ggsave("../FigProfitBasketMargin.pdf")

Fig_ProfitBasketMargin

df_profit <- cbind(df_profit, df_profit[,c("imbalanced","individual")]-df_profit[,"balanced"])
colnames(df_profit) <- c("BasketMargin", "profit_balanced", "profit_imbalanced","profit_supervised", "profit_none", "profit_all", "savings_imbalanced", "savings_supervised")
print(df_profit)

write.csv2(df_profit, file = "Table5_ProfitsSavings.csv")

### ATE Estimation ####
calc_ATE <- function(y, w, prop_score){
  #if(is.null(prop_score)){
  #  prop_score = rep(0.5, length(y))
  #}
  return( (sum(y*w/prop_score) - sum(y*(1-w)/(1-prop_score)) ) /length(y) )
}


# Register cluster for parallel processing
cl <- makeCluster(3)
registerDoParallel(cl)
RNGkind("L'Ecuyer-CMRG")  
clusterSetRNGStream(cl,iseed = 1234567)

# ATE for each experimental setting
ATE <- foreach(i=1:NO_EXPERIMENT_ITER,.combine = "rbind", .multicombine = TRUE, .packages = c("drtmle","SuperLearner"))%dopar%{
  ATE_temp <- c()
  ATE_temp["balanced"] <- calc_ATE(balanced[[i]]$y, balanced[[i]]$w, prop_score = 0.5)

  ATE_temp["imbalanced"] <- calc_ATE(imbalanced[[i]]$y, imbalanced[[i]]$w, prop_score = IMBALANCED_EXP_RATIO)

  
  ATE_temp["individual"] <- calc_ATE(individual[[i]]$y, individual[[i]]$w, individual[[i]]$prop_score)
  
  gn <- list()
  gn[[1]] <- individual[[i]]$prop_score
  gn[[2]] <- 1 - individual[[i]]$prop_score

  ATE_temp["individual_dr"] <-  ci(drtmle(Y=individual[[i]]$y,A=individual[[i]]$w,W=individual[[i]]$X,a_0 = c(1,0),
                                       family=binomial(),
                                       stratify=TRUE,
                                       SL_Q = c("SL.glm"),
                                       SL_g = c("SL.glm"),
                                       SL_Qr = "SL.glm",
                                       SL_gr = "SL.glm",
                                       gn=gn,
                                       maxIter = 1),contrast=c(1,-1))$drtmle[1]

  # True ATE
  ATE_temp["true_ATE"] <- mean(all[[i]]$y) - mean(none[[i]]$y)
  
  return(ATE_temp)
}


# Estimated ATE
ATE_hat <- apply(ATE,2,mean)
ATE_hat

ATE2 <- ATE[,c(5,1,2,3,4)]
ate_box <- melt(ATE2)

ggplot(data = ate_box, aes(x=Var2, y=value)) + 
  geom_boxplot() + 
  stat_summary(fun.y = "mean", geom = "point", colour = "blue", shape = 15, size = 2) +
  geom_hline(aes(yintercept=ATE_hat[5]),linetype="dashed") +
  labs(x="Experiment design", y = "ATE") +
  scale_x_discrete(labels=c("balanced" = "balanced", "imbalanced" = "imbalanced",
                            "individual" = "supervised (IPW1)", "individual_dr"= "supervised (DR)", "true_ATE"="true ATE")) +
  theme(axis.text.x = element_text(size=14), 
        axis.text.y = element_text(size=14),
        axis.title.x = element_text( size=16),
        axis.title.y = element_text(size=16))+ 
  ylim(c(-0.05,0.05))


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
lev_sample <- c( ATE$individual,ATE$individual_dr)

lev_group <- as.factor(c(rep("ind", length(ATE$individual)),rep("ind_dr", length(ATE$individual_dr))))
leveneTest(lev_sample,lev_group) # H0: Homogeneity of Variance 

                                        
                                    
#### CATE Estimation####

### Model evaluation ####

#cl <- makeCluster( max(1,detectCores()-1))
cl <- makeCluster(30)
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
  "individual" = do_experiment(X_train, expControl = expCtrl, prop_score = map_propensity(predict(targeting_model, X_train, type="response"), target_ratio=IMBALANCED_EXP_RATIO))
  )
  
  # Predictions from each model are saved in a list 'pred'
  pred <- list()
  
  for(rand_scheme in names(exp)){
  ## ATE baseline
  ate_hat <- calc_ATE(y = exp[[rand_scheme]]$y, w = exp[[rand_scheme]]$w, prop_score = exp[[rand_scheme]]$prop_score)
  pred[[paste0("ATE_",rand_scheme)]] <- rep(ate_hat, times = nrow(X_test) ) 
  
  ## Train t-learner
  t_logit <- T_Logit(X_train, exp[[rand_scheme]]$y, exp[[rand_scheme]]$w,
                     exp[[rand_scheme]]$prop_score)
  
  # Save predictions on test set in list pred
  pred[[paste0('t-logit_',rand_scheme)]] <- unname(predict(t_logit, X_test))
  
  ## Train causal forest
  cf <- grf::causal_forest(X=X_train,Y=exp[[rand_scheme]]$y, W=exp[[rand_scheme]]$w,
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

saveRDS(model_library, "../model_library_20190501.rds")
#model_library <- readRDS("../model_library_20190430.rds")

#### Calculate performance for each model ####
# New metric can be specified here
performance_CATE <- function(tau_score, y_true=NULL, w=NULL, prop_score=NULL, tau_true=NULL, value){
  res <- list()
  if(!is.null(tau_true))  res[["MAE"]] <- mean(abs(tau_true - tau_score))
  if(!is.null(y_true) & !is.null(w)) res[["Qini"]] <- qini_score(scores = tau_score, Y = y_true, W = w, p_treatment = prop_score)

  for(basket_value in value){
    res[[paste0("targeted_ratio",basket_value)]] <- mean(targeting_policy(
      tau_hat = tau_score, offer_cost = CONTACT_COST, customer_value = basket_value))
    res[[paste0("profit_",basket_value)]] <- catalogue_profit(y=y_true, contact_cost = CONTACT_COST, offer_cost = 0, value = basket_value,
                                 w=targeting_policy(
                                   tau_hat = tau_score, offer_cost = CONTACT_COST, customer_value = basket_value)
                                 )
  }
  
  return(as.data.table(res))
}

# Calculate performance per iteration
res <- lapply(model_library, function(ITER) lapply(ITER$pred, function(PRED){
  performance_CATE(tau_score= PRED,
                   y_true=ITER$true$y, w = ITER$true$w, prop_score = ITER$true$prop_score, tau_true = ITER$true$tau, value=VALUE)
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
res_KPI <- dcast(res[!grepl("profit",metric)&!grepl("targeted",metric),], metric~model, value.var = "perf_mean")

# CATE model targeted ratio
res_targeted <- dcast(res[grepl("targeted",metric),], metric~model, value.var = "perf_mean")

# CATE model performance for profit
res_profit <- dcast(res[grepl("profit",metric),], metric~model, value.var = "perf_mean")
res_profit <- as.data.frame(res_profit)
res_profit$metric <- gsub(as.character(res_profit$metric),pattern = "profit_",replacement = "")
res_profit[,-1] <- res_profit[,-1]/(0.5*N_CUSTOMER)

res_profit_net <- res_profit
# Calculate profit as net profit compared to balanced
for(model in c("ATE","CF","logit")){
  res_profit_net[,grepl(model, colnames(res_profit))] <- res_profit[,grepl(model, colnames(res_profit))] - 
                                                   unlist(res_profit[grepl(paste0(model,"_balanced"), colnames(res_profit))])
}

# Final tables
fwrite(res_targeted, "../CATE_targeted.txt")
fwrite(res_KPI, "../CATE_KPI.txt")
fwrite(res_profit_net, "../CATE_profit_net.txt")
fwrite(res_profit, "../CATE_profit.txt")

res_KPI
res_profit
res_profit_net

# t.test(perf_CATE$t_logit$balanced,perf_CATE$t_logit_DR$balanced)
# t.test(perf_CATE$t_logit_DR$balanced,perf_CATE$CF$balanced)
# 
# t.test(perf_CATE$t_logit$individual,perf_CATE$t_logit_DR$individual)
# t.test(perf_CATE$t_logit_DR$individual,perf_CATE$CF$individual)

