# Linex linear models 

# Takes as input an MSE-trained lm object, vector of true values, and per kwh 
#linex cost parameters. Returns optimal bias. 

fit_indirect_linex = function(SE_model, y, a, b){
  train_resids = -SE_model$residuals
  M = length(train_resids)
  beta_linex = -(1/a)*log((1/M)*sum(exp(a*train_resids)))
  return(beta_linex)
}

#Predict with indirectly biased linex 

predict_irradiance_indirect_linex = function(SE_model, newdata, bias){
  preds = predict(SE_model, newdata = newdata) + bias
  preds[which(preds < 0)] = 0
  return(preds)
}

#Calculate linex loss

linex_loss = function(yhat, y, a, b){
  resids = yhat - y
  resids = na.omit(resids)
  mean(b*(exp(a*resids)-a*resids - 1))
}

# Function to calculate gradient of LinEx cost with respect to parameters theta

linex_grad = function(y, X, theta, a, b){
  n = length(y)
  p = length(theta)
  fitted_vals = X%*%theta
  resids = fitted_vals - y
  gradient = apply(X, 2, function(x){a*b*sum(x*(exp(a*resids)-1))})
  return(gradient)
}

#One-observation gradient for SGD  (without apply so that computational time is not wasted)
linex_grad_sgd = function(y, X, theta, a, b){
  fitted_value = sum(X*theta)
  resid = fitted_value - y
  gradient = a*b*X*(exp(a*resid)-1)
  return(gradient)
}

#Estimate parameters by gradient descent 

#This function takes a model matrix as input (where instead of 1:s, the first column is 
#The 1 hour-ahead cosine zenith angle). The y-input is the hour-ahead solar irradiance. 

# Standardize variables to effectivize the algorithm. 

standardize  = function(x){std = (x-mean(x, na.rm = TRUE))/sd(x, na.rm = TRUE)
return(std)}

batch_gd <- function(y, model_matrix, a, b, eta, epochs, momentum){
  #Standardize the irradiance and lagged irradiances . 
  #  for(j in 2:ncol(model_matrix)){
  # model_matrix[,j] = standardize(model_matrix[,j])
  #}
  
  #Extract the y-values for which we have predictors to make fitted values, remove NA:s from 
  #model matrix 
  y = y[which(apply(model_matrix, 1, function(x){!any(is.na(x))}))]
  model_matrix = na.omit(model_matrix)
  
  # Setup output matrix
  results <- matrix(0.0, ncol = ncol(model_matrix) + 2L, nrow = epochs)
  colnames(results) <- c("epoch", "Loss", colnames(model_matrix))
  # Run algorithm
  theta <- rep(0.0, ncol(model_matrix)) # Init theta to the 0 vector
  change_prev = 0 # Initialize previous change for momentum 
  for(j in 1:epochs){
    gradient = linex_grad(y, model_matrix, theta, a, b)
    change_new = eta*gradient + momentum*change_prev
    theta = theta - change_new
    change_prev = change_new
    # Store epoch, LinEx and output results
    preds = model_matrix%*%theta
    results[j, "epoch"] <- j
    results[j, "Loss"] <- linex_loss(preds, y, a, b)
    results[j, -(1:2)] <- theta
  }
  return(as.data.frame(results))
}

### Try Stochastic gradient descent instead 

sgd_linex <- function(y, model_matrix, a, b, eta, epochs, momentum){
  #Standardize the irradiance and lagged irradiances for better convergence. 
  for(j in 2:ncol(model_matrix)){
    model_matrix[,j] = standardize(model_matrix[,j])
  }
  
  #Extract the y-values for which we have predictors to make fitted values, remove NA:s from 
  #model matrix (Applying listwise deletion)
  y = y[which(apply(model_matrix, 1, function(x){!any(is.na(x))}))]
  model_matrix = na.omit(model_matrix)
  
  # Setup output matrix
  results <- matrix(0.0, ncol = ncol(model_matrix) + 2L, nrow = epochs)
  colnames(results) <- c("epoch", "Loss", colnames(model_matrix))
  theta <- rep(0.0, ncol(model_matrix)) # Init theta to the 0 vector
  change_prev = 0 #Initialize previous change for momentum
  # Run algorithm
  for(j in 1:epochs){
    rows = sample(length(y))
    y = y[rows]
    model_matrix = model_matrix[rows,]
    for(i in 1:length(y)){
      gradient = linex_grad_sgd(y[i], model_matrix[i,], theta, a, b)
      change_new = eta*gradient + momentum*change_prev
      theta = theta - change_new
      change_prev = change_new
      #?    theta = theta - eta*gradient
    }
    # Store epoch, LinEx and output results
    preds = model_matrix%*%theta
    results[j, "epoch"] <- j
    results[j, "Loss"] <- linex_loss(preds, y, a, b)
    results[j, -(1:2)] <- theta
  }
  return(as.data.frame(results))
}

# Train directly biased model 

model_matrix = cbind(lulea_train_recent$cos_zenith_ahead, 
                     lulea_train_recent$normalized_irradiance)
colnames(model_matrix) = paste0("alpha_", c(0,1))
Y = lulea_train_recent$hour_ahead


linex_direct_model = batch_gd(Y, model_matrix, a_irr, b_en, 0.00001, 1000, 0)

plot(linex_direct_model$Loss)

model_matrix_test = cbind(lulea_test$cos_zenith_ahead, lulea_test$normalized_irradiance)

# Hour ahead irradiance predictions 
linex_irr_dir_pred = predict_irradiation(linex_direct_model, model_matrix_test, lulea_test$hour_ahead, 
                                         a_en, b_en)

# Hour ahead energy predictions
linex_en_dir_pred = predict_energy(linex_irr_dir_pred$Predictions, lulea_test, 
                            lead(lulea_test$energy_produced), 
                            a_en, b_en)

linex_en_dir_pred$Cost_Ratio 

hist_dir_linear_linex = linex_en_dir_pred$Residuals %>%
  as.data.frame() %>% setNames(c("X")) %>% 
  ggplot(aes(x = X)) + 
  geom_histogram(aes(X), bins = 50) +
  labs(x = "Error in Wh/m^2", title = "Directly 
optimized") +
  theme(axis.text = element_text(size = 18), axis.title = element_text(size = 21), 
        plot.title = element_text(size = 25))


# Indirect 

SE_model = lm(hour_ahead ~ 0 + cos_zenith_ahead + normalized_irradiance, 
              data = lulea_train_recent)

beta_linex = fit_indirect_linex(SE_model, lulea_train_recent$hour_ahead, 
                                a_irr, b_en)

linex_irr_ind_pred = predict_irradiance_indirect_linex(SE_model, lulea_test, 
                                                 beta_linex)

linex_en_ind_pred = predict_energy(linex_irr_ind_pred, lulea_test, 
                             lead(lulea_test$energy_produced), 
                             a_en, b_en)

linex_en_ind_pred$Cost_Ratio


hist_ind_linear_linex = linex_en_ind_pred$Residuals %>%
  as.data.frame() %>% setNames(c("X")) %>% 
  ggplot(aes(x = X)) + 
  geom_histogram(aes(X), bins = 50) +
  labs(x = "Error in Wh/m^2", title = "Indirectly 
optimized") +
  theme(axis.text = element_text(size = 18), axis.title = element_text(size = 21), 
        plot.title = element_text(size = 25))

#Unbiased 

irr_unb_pred = predict(SE_model, lulea_test)

en_unb_pred = predict_energy(irr_unb_pred, lulea_test, 
                             lead(lulea_test$energy_produced), 
                             a_en, b_en)

en_unb_pred$Cost_Ratio

hist_unb_linear_linex = en_unb_pred$Residuals %>%
  as.data.frame() %>% setNames(c("X")) %>% 
  ggplot(aes(x = X)) + 
  geom_histogram(aes(X), bins = 50) +
  labs(x = "Error in Wh/m^2", title = "Least 
squares") +
  theme(axis.text = element_text(size = 18), axis.title = element_text(size = 21), 
        plot.title = element_text(size = 25))

# Plots for all models residuals 

library(cowplot)
plot_grid(hist_unb_linear_linex, hist_ind_linear_linex, hist_dir_linear_linex, ncol = 3)

# Histogram with all in same 
  
library(tidyverse)

cbind(en_unb_pred$Residuals, en_ind_pred$Residuals, en_dir_pred$Residuals) %>% 
  as.data.frame() %>% setNames(c("Unbiased","Directly optimized","Indirectly optimized")) %>% 
  pivot_longer(cols = c(1:3), names_to = "Model") %>% 
  ggplot(aes(x = value, color = Model))  + 
  geom_histogram(aes(value))

# Perhaps better to present the plots separately.. 




