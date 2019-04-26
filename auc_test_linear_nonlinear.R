M = 50
auc_list <- matrix(NA,nrow=M,ncol=2)

N_VAR=20
N_CUSTOMER=1e5
#RATIO_SAMPLE=0.05

# Full set of customers
X_Customers <- sapply(1:N_VAR, function(x)rnorm(N_CUSTOMER, mean = 0, sd=1))
X <- data.frame(X_Customers)


for(i in 1:M){

expCtrl <- expControl(n_var = N_VAR, mode = "classification", beta_zero = -1.75,  # >0 indicates more than 50% purchasers
                      tau_zero =   0.75, # >0 indicates positive treatment effect)
                      DGP="nonlinear")

exp <- list()
exp$none <- do_experiment(X, expControl = expCtrl, prop_score = 0)


#### Examplary basic churn reponse model
# This can be any model that maps some business value to a treatment probability
# For churn, the churn probability is an obvious candidate, but could be rescaled
response_model <- glm(y~., cbind(X, y=exp$none$y), family = binomial(link="logit"))
churn_pred <- predict(response_model, X, type = "response")
auc_list[i,1] <-  ModelMetrics::auc(exp$none$y, churn_pred)

expCtrl <- expControl(n_var = N_VAR, mode = "classification", beta_zero = -1.75,  # >0 indicates more than 50% purchasers
                      tau_zero =   0.75, # >0 indicates positive treatment effect)
                      DGP="linear")

exp <- list()
exp$none <- do_experiment(X, expControl = expCtrl, prop_score = 0)


#### Examplary basic churn reponse model
# This can be any model that maps some business value to a treatment probability
# For churn, the churn probability is an obvious candidate, but could be rescaled
response_model <- glm(y~., cbind(X, y=exp$none$y), family = binomial(link="logit"))
churn_pred <- predict(response_model, X, type = "response")
auc_list[i,2] <-  ModelMetrics::auc(exp$none$y, churn_pred)

}

colnames(auc_list) <- c("nonlinear DGP", "linear DGP")
summary(auc_list)
