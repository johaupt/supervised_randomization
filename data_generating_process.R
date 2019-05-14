# #### Test  ####
# N_VAR=20
# N_CUSTOMER=10000
# EXPERIMENT_SIZE=N_CUSTOMER
# expCtrl <- expControl(n_var = N_VAR, mode = "classification", beta_zero = -3,  # >0 indicates more than 50% purchasers
#                       tau_zero =   0.425, # >0 indicates positive treatment effect)
#                       DGP="nonlinear")
# 
# # Location and scale of treatment effects
# round(quantile(probs=seq(0,1,0.05),
# do_experiment(
#   make_customers(EXPERIMENT_SIZE, 20),
#   expControl = expCtrl, prop_score = 1)$tau
# ),3)
# 
# mean(do_experiment(
#   make_customers(EXPERIMENT_SIZE, 20),
#   expControl = expCtrl, prop_score = 1)$tau)
# 
# # Response rate without treatment
# mean(probs=seq(0,1,0.1),
#          do_experiment(
#            make_customers(EXPERIMENT_SIZE, 20),
#            expControl = expCtrl, prop_score = 0)$y
# )
# 
# # Response rate with treatment
# mean(probs=seq(0,1,0.1),
#      do_experiment(
#        make_customers(EXPERIMENT_SIZE, 20),
#        expControl = expCtrl, prop_score = 1)$y
# )
# 

#### Experiment Functions ####
# Generate a population X
make_customers <- function(n_customer, n_var){
    X_Customers <- sapply(1:n_var, function(x)rnorm(n_customer, mean = 0, sd=1))
    return(data.frame(X_Customers))
}


# Define the experiment to be able to run it several times with 
# the same coefficients

expControl <- function(n_var, mode="regression", tau_zero=NULL, beta_zero=NULL,DGP="nonlinear",random_state=NULL){

  ###
  # n_var (int): Number of variables
  # model (str): 'regression' or binary 'classification'
  # tau_zero (float): Homogeneous/baseline treatment effect
  # beta_zero (float): Response constant
  ###
  if(!is.null(random_state)) set.seed(random_state)
  
  beta_zero = ifelse(is.null(beta_zero),0,beta_zero)
  beta = rnorm(n_var, mean = 0, sd = 0.5)
  beta_x2 = rnorm(n_var, mean = 0, sd = 1)
  beta_tau = rnorm(n_var, mean = 0, sd = 0.11)
  if(is.null(tau_zero)){
    tau_zero = rnorm(1, mean=0, sd=0.01)
  }
  
  return(list("beta_zero"=beta_zero, "beta"=beta, "beta_x2"=beta_x2, 
              "beta_tau"=beta_tau, "tau_zero" = tau_zero,
              "mode"=mode,"DGP"=DGP))
  
}

# Create data given the data generating process defined in 
# the experiment control object
do_experiment <- function(X, expControl=NULL, y=NULL, w=NULL, prop_score=NULL, X_out=TRUE, random_state=NULL){
  ###
  # X (array):
  #    Matrix of observations
  # g (vector of int): 
  #    Binary treatment assignment (0: control, 1:treatment) 
  # prop_score (single float or vector of floats): 
  #    Propensity score for random treatment. If g is not NULL it will override prop_score
  #
  ###
  if(is.null(expControl) & any(is.null(Y), is.null(W))){
    stop("Either expControl (simulation) or Y and W (usable obersvations) must be given.")
  }
  
  tau <- NULL
  class <- NULL
  logit <- function(x) 1/(1+exp(-x))
  n_obs = nrow(X)
  
  # X needs to be a matrix for calculation
  if(!"matrix" %in% class(X)){
    X <- as.matrix(X)
  }
  
  # Set random state
  if(!is.null(random_state)) set.seed(random_state)
  
  # Sample treatment group based on propensity
  g = rbinom(nrow(X), 1, prop_score)
  if(length(prop_score)==1){
      prop_score <- rep(prop_score, times=nrow(X))
  }

  
  # Simulate data if expControl is not NULL
  if(!is.null(expControl)){
  beta_zero = expControl$beta_zero
  beta = expControl$beta
  beta_x2 = expControl$beta_x2
  beta_tau = expControl$beta_tau
  tau_zero = expControl$tau_zero
  mode = expControl$mode
  DGP = expControl$DGP

  
  if(DGP %in% c('linear', 'nonlinear')){
  
    if(DGP=="linear"){


  if(mode %in% c('regression','classification')){
  
      
    if(mode == "regression"){
      # Functional form of treatment DGP
      tau = tau_zero + X%*%beta_tau + rnorm(n_obs, 0, 0.01)
      # Functional form of the response DGP
      y = beta_zero + X%*%beta + g*tau + rnorm(n_obs, 0, 0.5)
    }
    
    if(mode == "classification"){
      # Functional form of treatment DGP
      tau_logit = tau_zero + X%*%beta_tau + rnorm(n_obs, 0, 0.1)
      # Functional form of the response DGP
      p_logit = beta_zero + X%*%beta + rnorm(n_obs, 0, 1)
      
      y0 = logit(p_logit)
      y1 = logit(p_logit+tau_logit)
      tau = y1-y0
      
      p = logit(p_logit+g*tau_logit)
      y = rbinom(length(p),1,p)
      class = "linear"
    }
  }else{
    stop("Only mode 'classification' or 'regression' currently implemented")
  } 
    }
  
      if(DGP=="nonlinear"){
        
        if(mode %in% c('regression','classification')){
          
          if(mode == "regression"){
            tau = tau_zero + X%*%beta_tau + rnorm(n_obs, 0, 0.01)
            y = beta_zero + X%*%beta + g*tau + rnorm(n_obs, 0, 0.5)
          }
          
          if(mode == "classification"){

            tau_logit = tau_zero  + cos(X[,seq(2,20,2)])%*%beta_tau[seq(2,20,2)]     + sin(X[,seq(2,20,2)])%*%beta_tau[seq(2,20,2)] # + rnorm(n_obs, 0, 0.01)
            
            p_logit =   beta_zero + sin(X[,seq(1,19,2)])%*%beta[seq(1,19,2)]         + sin(X[,seq(2,20,2)])%*%beta_tau[seq(2,20,2)]
            #p_logit = beta_zero +  X[,seq(1,19,2)]%*%beta[seq(1,19,2)] + abs(X[,seq(1,19,2)]) %*% beta[seq(1,19,2)] + X[,seq(2,20,2)]%*%beta_tau[seq(2,20,2)]

            y0 = logit(p_logit)
            y1 = logit(p_logit+tau_logit)
            tau = y1-y0
            
            p = logit(p_logit+g*tau_logit)
            y = rbinom(length(p),1,p)

            class = "nonlinear"

          }
        }else{
          stop("Only mode 'classification' or 'regression' currently implemented")
        } 
      }
  } else{ stop("Choose between 'linear' or 'nonlinear'")}
  
  #### Sample from existing randomized controlled trial
  }else{
    usable <- W==g
    X <- X[usable,]
    y <- y[usable,]
    g <- g[usable,]
    prop_score <- prop_score[usable,]
  }
        
  
  if(X_out==FALSE){
    X = NULL
  }
  
  return(list("X"=X, "y"=y, "g"=g, "tau"=tau, "prop_score"=prop_score, "class"=class))
  
}

