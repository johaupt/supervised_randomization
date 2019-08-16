
#### Optimal threshold for targeting ####
opt_tau <- function(offer_cost, customer_value){
  threshold <- offer_cost/customer_value
  return(threshold)
}

#### Optimal binary targeting based on optimal threshold given cost setting ####
targeting_policy <- function(tau_hat, offer_cost, customer_value){
  threshold <- offer_cost/customer_value
  treatment <- as.numeric(tau_hat > threshold)
  return(treatment)
}

top_decile_targeting_policy <- function(tau_hat){
  tau_hat_order <- order(tau_hat, decreasing = TRUE)
  cutoff <- floor(0.1 * length(tau_hat))
  treatment <- ifelse(tau_hat_order <= cutoff, 1, 0)
  return(treatment)
}

#### Transformed outcome loss ####
transformed_outcome_loss <- function(treatment_effect, Y, W, p_treatment){
  ###
  # Biased estimate of the MSE (tau - tau_hat)^2 via the transformed outcome
  # See Athey & Imbens (2016) Recursive Partitioning for Heterogeneous Causal Effects
  #     Hitsch & Misra (2018): Heterogeneous Treatment Effects and Optimal Targeting Policy Evaluation
  ###
  Y_transformed <- Y * (W - p_treatment) / (p_treatment * (1-p_treatment))
  return( mean( (Y_transformed - treatment_effect)^2) )
}


#### Churn costs ####
churn_cost <- function(y, w, contact_cost, offer_cost, value){
  total <- sum((1-w)*(1-y)) * 0 +
    sum(w*y)         * -(contact_cost + customer_value) +
    sum(w*(1-y))     * -(contact_cost + offer_cost)+
    sum((1-w)*y)     * -customer_value
  
  return(total)
  
  #(mean(w)*               -cost_treatment_fix +
  #mean(w)*(1-mean(y))*   -cost_treatment_var +
  #mean(y)*               -cost_churn)*
  #n_customer
}

#### Catalogue Campaign Profit ####
catalogue_profit <- function(y, w, contact_cost, offer_cost=0, value, mode="total"){
  if(!mode %in% c("total","individual")){
    stop("Argument mode must be one of ['total', 'individual]")
  }
  
  profit <- (
    # Not treated/no purchase
    (1-w)*(1-y) * 0 +
    # Treated/purchase
    w*y         * (value - contact_cost) +
    # Treated/no purchase
    w*(1-y)     * (-contact_cost)+
    # Not treated/purchase
    (1-w)*y     * (value)
  )
  
  if(mode=="total"){
      profit <- sum(profit)
  }
  
  return(profit)
}

#### Expected Profit on Usable Observations
expected_profit <- function(policy_decision, profit, treatment_group, p_treatment){
  ### 
  # Calculate expected profit of a candidate targeting policy on randomized trial data.
  # See Hitsch & Misra (2018): Heterogeneous Treatment Effects and Optimal Targeting Policy Evaluation for details
  ###
  
  # Intuitively: When treatment group matches policy decision then add the profit corrected by the prob. that the 
  # random treatment group matches the policy decision (either e(x) or 1-e(x) )
  return( sum((1-treatment_group)*(1-policy_decision)/(1-p_treatment)*profit + (treatment_group)*(policy_decision)/(p_treatment)*profit) )
}

  
  
  