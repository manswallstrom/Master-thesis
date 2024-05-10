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

#Estimate parameters by gradient descent 
#This function takes a model matrix as input (where instead of 1:s, the first column is 
#The 1 hour-ahead cosine zenith angle). The y-input is the hour-ahead solar irradiance. 

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

# Perhaps better to present the plots separately.. 




