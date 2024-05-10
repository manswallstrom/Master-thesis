# Linex_linear gbg
library(tidyverse)
# Train directly biased model 

model_matrix = cbind(gbg_train$cos_zenith_ahead, 
                     gbg_train$normalized_irradiance)
colnames(model_matrix) = paste0("alpha_", c(0,1))
Y = gbg_train$hour_ahead

linex_direct_model = batch_gd(Y, model_matrix, a_irr, b_en, 0.00001, 1000, 0)

plot(linex_direct_model$Loss)

model_matrix_test = cbind(gbg_test$cos_zenith_ahead, gbg_test$normalized_irradiance)

# Hour ahead irradiance predictions 
irr_dir_pred = predict_irradiation(linex_direct_model, model_matrix_test, gbg_test$hour_ahead, 
                                         a_en, b_en)

# Hour ahead energy predictions
en_dir_pred = predict_energy(irr_dir_pred$Predictions, gbg_test, 
                                   lead(gbg_test$energy_produced), 
                                   a_en, b_en)

en_dir_pred$cost_ratio_linex 

# Indirect 

SE_model = lm(hour_ahead ~ 0 + cos_zenith_ahead + normalized_irradiance, 
              data = gbg_train)

beta_linex = fit_indirect_linex(SE_model, gbg_train$hour_ahead, 
                                a_irr, b_en)

irr_ind_pred = predict_irradiance_indirect_linex(SE_model, gbg_test, 
                                                       beta_linex)

en_ind_pred = predict_energy(irr_ind_pred, gbg_test, 
                                   lead(gbg_test$energy_produced), 
                                   a_en, b_en)

en_ind_pred$cost_ratio_linex 

#Unbiased 

irr_se_pred = predict(SE_model, gbg_test)

en_se_pred = predict_energy(irr_se_pred, gbg_test, 
                             lead(gbg_test$energy_produced), 
                             a = a_en, b = b_en, C = cpwl_slopes_en, D = D_cpwl_en, 
                             intercepts = intercepts_cpwl_en)

en_se_pred$cost_ratio_linex 


# Regression trees 

X_tree_train = cbind(gbg_train$normalized_irradiance, gbg_train$cos_zenith_ahead)
y_tree_train = gbg_train$hour_ahead

library(foreach)
library(doParallel)

cl = makeCluster(7)
registerDoParallel(cl)

start.time <- Sys.time()
tree_dir = train_random_forest(X_tree_train, y_tree_train, 100, 100, 2, linex_loss_rf,
                               a = a_irr, b = b_en)
end.time <- Sys.time()
time.taken <- end.time - start.time

stopCluster(cl)

X_tree_test = cbind(gbg_test$normalized_irradiance, gbg_test$cos_zenith_ahead)

irr_tree_dir_pred = predict_with_bagged_trees(X_tree_test, tree_dir)

en_tree_dir_pred = predict_energy(irr_tree_dir_pred, gbg_test, 
                                  lead(gbg_test$energy_produced), a = a_en, b = b_en)

en_tree_dir_pred$cost_ratio_linex

# SE regression tree model 

library(dplyr)

cl = makeCluster(7)

registerDoParallel(cl)

start.time <- Sys.time()
tree_se = train_random_forest(X_tree_train, y_tree_train, 250, 100, 2, se_loss_rf)
end.time <- Sys.time()
time.taken <- end.time - start.time

stopCluster(cl)

irr_tree_se_pred = predict_with_bagged_trees(X_tree_test, tree_se)

en_tree_se_pred = predict_energy(irr_tree_se_pred, gbg_test, 
                                lead(gbg_test$energy_produced), a = a_en, b = b_en)

en_tree_se_pred$cost_ratio_linex 

# Indirectly optimized tree model 

# Input training data predictions from SE tree model and y 
# we dont lead y because it is the "hour ahead" variable. 

se_tree_in_sample = predict_with_bagged_trees(X_tree_train, tree_se)

fit_indirect_linex_tree = function(preds, y, a, b){
  train_resids = preds - y 
  M = length(train_resids)
  beta_linex = -(1/a)*log((1/M)*sum(exp(a*train_resids)))
  return(beta_linex)
}

beta_linex_tree = fit_indirect_linex_tree(se_tree_in_sample, y_tree_train, a_irr, b_en)

# Input out of sample predictions from unbiased model, and the bias
# This function just adds the bias and sets negative predictions to 0 

predict_irradiance_indirect_linex_tree = function(preds, bias){
  preds = preds + bias
  preds[which(preds < 0)] = 0
  return(preds)
}

irr_tree_ind_pred = predict_irradiance_indirect_linex_tree(irr_tree_se_pred,beta_linex_tree)

en_tree_ind_pred = predict_energy(irr_tree_ind_pred, gbg_test, 
                                        lead(gbg_test$energy_produced), 
                                        a_en, b_en)

en_tree_ind_pred$cost_ratio_linex


## Plot of prediction errors for all models and methods 

as.data.frame(
  cbind(en_se_pred$Residuals, en_ind_pred$Residuals, en_dir_pred$Residuals, 
        en_tree_se_pred$Residuals, en_tree_ind_pred$Residuals, en_tree_dir_pred$Residuals)) %>%
  setNames(c("Linear_Least squares", "Linear_Indirect", "Linear_Direct", "Regression tree_Least squares",
             "Regression tree_Indirect", "Regression tree_Direct")) %>% 
  pivot_longer(cols = 1:6, names_to = c("Model", "Method"), names_sep = "_", 
               values_to = "Error") %>% 
  ggplot(aes(x = Error)) + 
  geom_histogram(bins = 50)  + 
  geom_vline(xintercept = 0) +
  facet_grid(rows = vars(Model), 
             cols = vars(factor(Method, levels = c("Least squares", "Indirect", "Direct")))) + 
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 20), 
        strip.text = element_text(size = 14))

library(xtable)
xtable(cbind(en_se_pred$cost_ratio_linex, en_ind_pred$cost_ratio_linex, en_dir_pred$cost_ratio_linex, 
             en_tree_se_pred$cost_ratio_linex, en_tree_ind_pred$cost_ratio_linex, 
             en_tree_dir_pred$cost_ratio_linex))


