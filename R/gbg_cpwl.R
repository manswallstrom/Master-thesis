# gbg cpwl

# Package  CVXR 

# Batch linear model

# CPWL directly biased forecast

# Note, only works for four-piece CPWL 
# Sum of CPWL losses 


CPWL_loss = function(c, d, pred, true){
  b1 = c[2]*d[1]  - c[1]*d[1]  
  b2 = 0
  b3 = 0
  b4 = c[3]*d[3]  - c[4]*d[3] 
  b_matr = cbind(rep(b1, length(true)), rep(b2, length(true)), rep(b3, length(true)), 
                 rep(b4, length(true)))
  resids = pred-true
  losses = resids%*%t(as.matrix(c)) + b_matr
  loss = apply(losses, 1, max)
  return(list(loss = loss, b = b_matr[1,]))
}

library(tidyverse)
library(CVXR)

alpha_hat = Variable(2)

model_matrix = cbind(gbg_train$cos_zenith_ahead, gbg_train$normalized_irradiance)
Y = gbg_train$hour_ahead

colnames(model_matrix) = paste0("alpha_", c(0,1))

b_matrix = cbind(rep(intercepts_cpwl_irr[1], length(Y)), rep(intercepts_cpwl_irr[2], length(Y)), 
                 rep(intercepts_cpwl_irr[3], length(Y)), rep(intercepts_cpwl_irr[4], length(Y)))

w = max_entries((model_matrix%*%alpha_hat - as.matrix(Y))%*%t(as.matrix(cpwl_slopes_irr)) + b_matrix,
                axis = 1)

objective = Minimize(sum(w))

loss_constraint_1 = (model_matrix%*%alpha_hat - as.matrix(Y))*cpwl_slopes_irr[1] + intercepts_cpwl_irr[1] <= w

loss_constraint_2 = (model_matrix%*%alpha_hat - as.matrix(Y))*cpwl_slopes_irr[2] + intercepts_cpwl_irr[2]  <= w

loss_constraint_3 = (model_matrix%*%alpha_hat - as.matrix(Y))*cpwl_slopes_irr[3] + intercepts_cpwl_irr[3]  <= w

loss_constraint_4 = (model_matrix%*%alpha_hat - as.matrix(Y))*cpwl_slopes_irr[4] + intercepts_cpwl_irr[4]  <= w

prob = Problem(objective)

cpwl_solve = solve(prob, solver = "SCS")

model_matrix_test = cbind(gbg_test$cos_zenith_ahead, gbg_test$normalized_irradiance)

#hour ahead irradiance predictions 
cpwl_irr_dir_pred = model_matrix_test%*%cpwl_solve$getValue(alpha_hat)

#hour ahead energy predictions
cpwl_en_dir_pred = predict_energy(cpwl_irr_dir_pred, gbg_test, 
                                  lead(gbg_test$energy_produced), C = cpwl_slopes_en, 
                                  D = D_cpwl_en, intercepts = intercepts_cpwl_en)

cpwl_en_dir_pred$cost_ratio_cpwl

SE_model = lm(hour_ahead ~ 0 + cos_zenith_ahead + normalized_irradiance, 
              data = gbg_train)

cpwl_indirect_bias_linear = CPWL_indirect_bias(cpwl_slopes_irr,D_cpwl_irr, SE_model$fitted.values, 
                                               gbg_train$hour_ahead)

# poorly named function predict_irradiance_indirect_linex

cpwl_irr_ind_pred = predict_irradiance_indirect_linex(SE_model, gbg_test, 
                                                      cpwl_indirect_bias_linear)

cpwl_en_ind_pred = predict_energy(cpwl_irr_ind_pred, gbg_test, 
                                         lead(gbg_test$energy_produced),C = cpwl_slopes_en, 
                                         D = D_cpwl_en, intercepts = intercepts_cpwl_en)

#CPWL loss ratio 

cpwl_en_ind_pred$cost_ratio_cpwl

#CPWL loss ratio unbiased model

en_se_pred$cost_ratio_cpwl
  
# CPWL tree models
  
# Regression trees 

X_tree_train = cbind(gbg_train$normalized_irradiance, gbg_train$cos_zenith_ahead)
y_tree_train = gbg_train$hour_ahead

View(gbg_train)

library(foreach)
library(doParallel)

cl = makeCluster(7)
registerDoParallel(cl)

start.time <- Sys.time()
cpwl_tree_dir = train_random_forest(X_tree_train, y_tree_train, 250, 100, 2, cpwl_loss_rf, a = a_irr, b_en, 
                                    C = cpwl_slopes_irr, D = D_cpwl_irr, intercepts = intercepts_cpwl_irr)
end.time <- Sys.time()
time.taken <- end.time - start.time

stopCluster(cl)

X_tree_test = cbind(gbg_test$normalized_irradiance, gbg_test$cos_zenith_ahead)

cpwl_tree_irr_dir_pred = predict_with_bagged_trees(X_tree_test, cpwl_tree_dir)

cpwl_tree_en_dir_pred = predict_energy(cpwl_tree_irr_dir_pred, gbg_test, 
                                      lead(gbg_test$energy_produced), C = cpwl_slopes_en, 
                                      D = D_cpwl_en, intercepts = intercepts_cpwl_en)

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
                                lead(gbg_test$energy_produced), a = a_en, b = b_en, 
                                C = cpwl_slopes_en, D = D_cpwl_en, intercepts = intercepts_cpwl_en)

en_tree_se_pred$cost_ratio_cpwl 

# Indirectly optimized tree model 

# Input training data predictions from SE tree model and y 
# we dont lead y because it is the "hour ahead" variable. 

se_tree_in_sample = predict_with_bagged_trees(X_tree_train, tree_se)

cpwl_indirect_bias_linear = CPWL_indirect_bias(cpwl_slopes_irr,D_cpwl_irr, se_tree_in_sample, 
                                               gbg_train$hour_ahead)

# Input out of sample predictions from unbiased model, and the bias
# This function just adds the bias and sets negative predictions to 0 

predict_irradiance_indirect_linex_tree = function(preds, bias){
  preds = preds + bias
  preds[which(preds < 0)] = 0
  return(preds)
}

cpwl_irr_tree_ind_pred = predict_irradiance_indirect_linex_tree(irr_tree_se_pred,
                                                           cpwl_indirect_bias_linear)

cpwl_en_tree_ind_pred = predict_energy(cpwl_irr_tree_ind_pred, gbg_test, 
                                       lead(gbg_test$energy_produced), 
                                       C = cpwl_slopes_en, D = D_cpwl_en, 
                                       intercepts = intercepts_cpwl_en)

cpwl_en_tree_ind_pred$cost_ratio_cpwl


## Plot of prediction errors for all models and methods 

as.data.frame(
  cbind(en_se_pred$Residuals, cpwl_en_ind_pred$Residuals, cpwl_en_dir_pred$Residuals, 
        en_tree_se_pred$Residuals, cpwl_en_tree_ind_pred$Residuals, cpwl_tree_en_dir_pred$Residuals)) %>%
  setNames(c("Linear_Least squares", "Linear_Indirect", "Linear_Direct", "Regression tree_Least squares",
             "Regression tree_Indirect", "Regression tree_Direct")) %>% 
  pivot_longer(cols = 1:6, names_to = c("Model", "Method"), names_sep = "_", 
               values_to = "Error") %>% 
  ggplot(aes(x = Error)) + 
  geom_histogram(bins = 50)  + 
  geom_vline(xintercept = D_cpwl_en) + 
  facet_grid(rows = vars(Model), 
             cols = vars(factor(Method, levels = c("Least squares", "Indirect", "Direct")))) + 
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 20), 
        strip.text = element_text(size = 14))

library(xtable)
xtable(cbind(en_se_pred$cost_ratio_cpwl, cpwl_en_ind_pred$cost_ratio_cpwl, cpwl_en_dir_pred$cost_ratio_cpwl, 
             en_tree_se_pred$cost_ratio_cpwl, cpwl_en_tree_ind_pred$cost_ratio_cpwl, cpwl_tree_en_dir_pred$cost_ratio_cpwl))

