#### Optimal threshold for targeting ####
opt_tau <- function(offer_cost, customer_value){
  threshold <- offer_cost/customer_value
  return(threshold)
}


#### Churn costs ####
churn_cost <- function(y, g, contact_cost, offer_cost, value){
  total <- sum((1-g)*(1-y)) * 0 +
    sum(g*y)         * -(contact_cost + customer_value) +
    sum(g*(1-y))     * -(contact_cost + offer_cost)+
    sum((1-g)*y)     * -customer_value
  
  return(total)
  
  #(mean(g)*               -cost_treatment_fix +
  #mean(g)*(1-mean(y))*   -cost_treatment_var +
  #mean(y)*               -cost_churn)*
  #n_customer
}

#### Catalogue Campaign Profit ####

catalogue_profit <- function(y, g, contact_cost, offer_cost=0, value){
  total <- 
    # Not treated/no purchase
    sum((1-g)*(1-y)) * 0 +
    # Treated/purchase
    sum(g*y)         * (basket_value - contact_cost) +
    # Treated/no purchase
    sum(g*(1-y))     * (-contact_cost)+
    # Not treated/purchase
    sum((1-g)*y)     * (basket_value)
  
  return(total)
}