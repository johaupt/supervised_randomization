#### Churn costs per customer ####
churn_cost <- function(y, g, COST_TREATMENT_FIX, 
                       COST_TREATMENT_VAR, COST_CHURN){
  total <- sum((1-g)*(1-y)) * 0 +
    sum(g*y)         * -(COST_TREATMENT_FIX + COST_CHURN) +
    sum(g*(1-y))     * -(COST_TREATMENT_FIX + COST_TREATMENT_VAR)+
    sum((1-g)*y)     * -COST_CHURN
  
  return(total)
  
  #(mean(g)*               -cost_treatment_fix +
  #mean(g)*(1-mean(y))*   -cost_treatment_var +
  #mean(y)*               -cost_churn)*
  #n_customer
}

#### Catalogue Campaign Profit per customer ####
# Expected churn costs per customer
catalogue_cost <- function(y, g, COST_TREATMENT_FIX, 
                       COST_TREATMENT_VAR, COST_CHURN){
  total <- sum((1-g)*(1-y)) * 0 +
    sum(g*y)         * -(COST_TREATMENT_FIX + COST_CHURN) +
    sum(g*(1-y))     * -(COST_TREATMENT_FIX + COST_TREATMENT_VAR)+
    sum((1-g)*y)     * -COST_CHURN
  
  return(total)
  
  #(mean(g)*               -cost_treatment_fix +
  #mean(g)*(1-mean(y))*   -cost_treatment_var +
  #mean(y)*               -cost_churn)*
  #n_customer
}