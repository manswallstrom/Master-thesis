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

library(CVXR)

alpha_hat = Variable(2)

model_matrix = cbind(lulea_train_recent$cos_zenith_ahead, lulea_train_recent$normalized_irradiance)
Y = lulea_train_recent$hour_ahead

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

model_matrix_test = cbind(lulea_test$cos_zenith_ahead, lulea_test$normalized_irradiance)

#hour ahead irradiance predictions 
cpwl_preds = model_matrix_test%*%cpwl_solve$getValue(alpha_hat)

#hour ahead energy predictions
cpwl_dir_en_pred = predict_energy(cpwl_preds, lulea_test, 
                                  lead(lulea_test$energy_produced), 0, 0)

hist_cpwl_dir = cpwl_dir_en_pred$Residuals %>%
  as.data.frame() %>% setNames(c("X")) %>% 
  ggplot(aes(x = X)) + 
  geom_histogram(aes(X), bins = 50) +
  labs(x = "Error in Wh/m^2", title = "Directly 
optimized") +
  theme(axis.text = element_text(size = 18), axis.title = element_text(size = 21), 
        plot.title = element_text(size = 25))

sum(CPWL_loss(cpwl_slopes_en, D_cpwl_en, cpwl_dir_en_pred$Predictions, lead(lulea_test$energy_produced))$loss, 
    na.rm = TRUE)/
  sum(CPWL_loss(cpwl_slopes_en, D_cpwl_en, 0, lead(lulea_test$energy_produced))$loss, 
      na.rm = TRUE)






