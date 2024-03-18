# Indirectly biased CPWL

#Input: Vector of CPWL slopes (C), vector of CPWL breakponts (D), 
#predicted values form SE model and true values
 
CPWL_indirect_bias = function(C, D, preds, y){
  residuals = y - preds
  empirical_cdf = ecdf(residuals) 
  tmp = numeric(0) 
  evaluate_function = function(beta){  # Calculate value of the function to optimize for suggested beta
    for(i in 1:(length(C)-1)){
      tmp[i] = (C[i]-C[i+1])*empirical_cdf(D[i] - beta)
    }
    return(sum(tmp) + C[length(C)])
  }
  a = -1000  # Initialize to some point where we know function will be negative
  b = 1000 # Some point where we know function will be positive 
  iter = 1
  f_a = evaluate_function(a)
  while(iter < 10000){
    beta = (a+b)/2 # midpoint of a and b 
    f_beta = evaluate_function(beta)
    if(f_beta == 0|(b-a)/2 < 0.001){break} # break if we have found optimum or if a and b are very close
    if(f_a*f_beta > 0){a = beta         #if f_a and f_beta have the same sign
    f_a = f_beta}else{b = beta}   
    iter = iter + 1
  }
  if(iter == 10000){warning("Algorithm did not converge")}
  return(beta)
}

cpwl_indirect_bias_linear = CPWL_indirect_bias(cpwl_slopes_irr,D_cpwl_irr, SE_model$fitted.values, 
                                               lulea_train_recent$hour_ahead)

# Use linex predict function bcause it doesnt matter
irr_ind_pred_cpwl_linear = predict_irradiance_indirect_linex(SE_model, lulea_test, 
                                                             cpwl_indirect_bias_linear)

en_ind_pred_cpwl_linear = predict_energy(irr_ind_pred_cpwl_linear, lulea_test, 
                             lead(lulea_test$energy_produced), 
                             a_en, b_en)

#CPWL loss ratio 

sum(CPWL_loss(cpwl_slopes_en, D_cpwl_en, en_ind_pred_cpwl_linear$Predictions, lead(lulea_test$energy_produced))$loss, 
    na.rm = TRUE)/sum(CPWL_loss(cpwl_slopes_en, D_cpwl_en, 0, lead(lulea_test$energy_produced))$loss, 
                      na.rm = TRUE)

#CPWL loss ratio unbiased model

sum(CPWL_loss(cpwl_slopes_en, D_cpwl_en, en_unb_pred$Predictions, lead(lulea_test$energy_produced))$loss, 
    na.rm = TRUE)/sum(CPWL_loss(cpwl_slopes_en, D_cpwl_en, 0, lead(lulea_test$energy_produced))$loss, 
                      na.rm = TRUE)


# CPWL ind histogram 
library(ggplot2)

hist_cpwl_ind_linear = en_ind_pred_cpwl_linear$Residuals %>%
  as.data.frame() %>% setNames(c("X")) %>% 
  ggplot(aes(x = X)) + 
  geom_histogram(aes(X), bins = 50) +
  labs(x = "Error in Wh/m^2", title = "Indirectly 
optimized") +
  theme(axis.text = element_text(size = 18), axis.title = element_text(size = 21), 
        plot.title = element_text(size = 25))


# plot of cpwl linear

library(cowplot)

plot_grid(hist_unb_linear_linex, hist_cpwl_ind_linear, hist_cpwl_dir, ncol = 3)


# CPWL Tree based indirect 

cpwl_bias_tree = CPWL_indirect_bias(cpwl_slopes_irr, D_cpwl_irr, se_tree_in_sample, 
                   df_2021$hour_ahead)


# Predict irradiance 
irr_ind_pred_cpwl_tree = predict_irradiance_indirect_linex_tree(irr_unb_rf_pred,
                                                                cpwl_bias_tree)

en_ind_pred_cpwl_tree = predict_energy(irr_ind_pred_cpwl_tree, lulea_test, 
                                         lead(lulea_test$energy_produced), 
                                         a_en, b_en)

#CPWL loss ratio 

sum(CPWL_loss(cpwl_slopes_en, D_cpwl_en, en_ind_pred_cpwl_tree$Predictions, lead(lulea_test$energy_produced))$loss, 
    na.rm = TRUE)/sum(CPWL_loss(cpwl_slopes_en, D_cpwl_en, 0, lead(lulea_test$energy_produced))$loss, 
                      na.rm = TRUE)

#CPWL loss ratio unbiased model

sum(CPWL_loss(cpwl_slopes_en, D_cpwl_en, en_unb_rf_pred$Predictions, lead(lulea_test$energy_produced))$loss, 
    na.rm = TRUE)/sum(CPWL_loss(cpwl_slopes_en, D_cpwl_en, 0, lead(lulea_test$energy_produced))$loss, 
                      na.rm = TRUE)


# CPWL ind histogram 

hist_cpwl_ind_tree = en_ind_pred_cpwl_tree$Residuals %>%
  as.data.frame() %>% setNames(c("X")) %>% 
  ggplot(aes(x = X)) + 
  geom_histogram(aes(X), bins = 50) +
  labs(x = "Error in Wh/m^2", title = "Indirectly 
optimized") +
  theme(axis.text = element_text(size = 18), axis.title = element_text(size = 21), 
        plot.title = element_text(size = 25))

plot_grid(hist_unb_rf, hist_cpwl_ind_tree, ncol = 3)






