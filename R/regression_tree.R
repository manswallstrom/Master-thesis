# Regression tree algorithm 

tree_split <- function(X, y, l, loss, a = NULL, b = NULL){
  checkmate::assert_matrix(X)
  checkmate::assert_numeric(y, len = nrow(X))
  checkmate::assert_int(l)
  if(nrow(X)<2*l){
    stop("nrow(X)<2*l")
  }
  
  losses <- matrix(Inf, nrow = nrow(X), ncol = ncol(X))
  split_point = numeric(2)
  for(j in 1:ncol(X)){
    for(k in 1:nrow(X)){
      s <- X[k,j] # Split point
      R1 <- which(X[,j] <= s)
      R2 <- which(X[,j] > s)
      
      if(length(R1)>=l & length(R2)>=l){ 
        #If one of the regions is smaller than the minimum node size, let loss remain inf so it is not selected as the split point
        
        losses[k, j] <- loss(y[R1], y[R2], a, b)
        
        if(losses[k,j] == min(losses)){
          split_point[1] = j
          split_point[2] = s #update optimal split point (we will end up with the lowest one)
        }
        
      } else{next}
    }
  }
  
  R1 <- which(X[,split_point[1]] <= split_point[2])         
  R2 <- which(X[,split_point[1]] > split_point[2]) 
  list(j = split_point[1],
       s = split_point[2],
       R1 = R1,
       R2 = R2,
       loss_at_split = min(losses))
}

grow_tree <- function(X, y, l, loss){
  checkmate::assert_matrix(X)
  checkmate::assert_numeric(y, len = nrow(X))
  checkmate::assert_int(l)
  init = tree_split(X, y, l, loss, a, b)
  S_m <- list(init$R1, init$R2)
  results <- data.frame(j = init$j,
                        s = init$s,
                        R1_i = 2, #In the initial row of my results matrix, representing the first split,
                        R2_i = 3,  #We go to row 2 to follow region 1 and row 3 to follow region 2 
                        gamma = NA)
  m = 1    # Keeps track of iteration of the while loop, index used to acces regions
  M = 2    # This index is also used to append and access regions
  while (m<=M) {   #This will stop being true eventually when there are no more possible splits
    
    #The following lines are used to determine whether to split.
    #This condition instead of the suggested length(S[[m]]) >= 2*l takes into 
    #account that sometimes even vectors where length(S[[m]]) >= 2*l cant be split 
    #if the variables are discrete since there may be several observations at the exact splitting value, 
    #making it possible for one region to be smaller than the smallest allowed number of leaves.
    
    truth_vector = logical(ncol(X))
    for(j in 1:ncol(X)){
      truth_vector[j] = 
        isTRUE(length(which(X[S_m[[m]],j] <= median(X[S_m[[m]],j]))) >=l & length(which(X[S_m[[m]],j] > median(X[S_m[[m]],j]))) >=l)
    }
    if(any(truth_vector)){
      
      new_split = tree_split(X[S_m[[m]], ], y[S_m[[m]]], l, loss, a, b)
      
      # Get indices on total data form
      S_m[[M+1]] = S_m[[m]][new_split$R1] 
      S_m[[M+2]] = S_m[[m]][new_split$R2]
      
      #Append split variables and values to result matrix
      results[m+1,1] = new_split$j
      results[m+1,2] = new_split$s
      
      #Increment R1_i and R2_i index 
      results[m+1,3] = M + 2
      results[m+1,4] = M + 3
      
      M = M + 2
      
    } else{
      results[m+1,5] = mean(y[S_m[[m]]])
    }
    m = m + 1
  }
  return(round(results,1))
}

pred = function(obs, tree){
  terminal_node = FALSE
  row = 1
  while(terminal_node == FALSE){
    
    if(obs[tree$j[row]] <= tree$s[row]){
      row = tree[row,"R1_i"]
    } else if(obs[tree$j[row]] > tree$s[row]){ 
      row = tree[row,"R2_i"]
    }
    terminal_node = isFALSE(is.na(tree$gamma[row]))
  }
  return(tree$gamma[row])
} 

predict_with_tree = function(new_data,tree){
  checkmate::assert_matrix(new_data)
  predictions <- apply(new_data, 1, function(x){pred(x, tree)})
  return(predictions)
}


# Linex loss in RF split format 
linex_loss_rf = function(R1, R2, a, b){
  N_region_1 = length(R1)
  N_region_2 = length(R2)
  within_region_min_1 = (log(N_region_1/sum(exp(-a*R1))))/a 
  within_region_min_2 = (log(N_region_2/sum(exp(-a*R2))))/a 
  errors = c(within_region_min_1 - R1, within_region_min_2 - R2)
  loss = sum(b*(exp(a*errors)-a*errors-1))
  names(loss) = "linex"
  return(loss)
}

# SE loss for RF split 

se_loss_rf = function(R1, R2, a = NULL, b = NULL){
  within_region_min_1 = mean(R1) 
  within_region_min_2 = mean(R2)
  errors = c(within_region_min_1 - R1, within_region_min_2 - R2)
  loss = sum(errors^2)
  names(loss) = "se"
  return(loss)
}

# CPWL loss in RF split format (computationally infeasible)

cpwl_loss_rf = function(R1, R2, C, D, B){
  
  b_matrix_1 = cbind(rep(B[1], length(R1)), rep(B[2], length(R1)), 
                   rep(B[3], length(R1)), rep(B[4], length(R1)))
  alpha_hat_1 = Variable(1)
  w_1 = max_entries((alpha_hat_1 - as.matrix(R1))%*%t(as.matrix(C)) + b_matrix_1,
                  axis = 1)
  
  objective_1 = Minimize(sum(w_1))
  r1_loss_constraint_1 = (alpha_hat_1 - as.matrix(R1))*C[1] + b_vector[1] <= w_1
  r1_loss_constraint_2 = (alpha_hat_1 - as.matrix(R1))*C[2] + b_vector[2]  <= w_1
  r1_loss_constraint_3 = (alpha_hat_1 - as.matrix(R1))*C[3] + b_vector[3]  <= w_1
  r1_loss_constraint_4 = (alpha_hat_1 - as.matrix(R1))*C[4] + b_vector[4]  <= w_1
  
  prob_1 = Problem(objective_1)
  cpwl_solve_1 = solve(prob_1, solver = "SCS")
  
  b_matrix_2 = cbind(rep(B[1], length(R2)), rep(B[2], length(R2)), 
                     rep(B[3], length(R2)), rep(B[4], length(R2)))
  alpha_hat_2 = Variable(1)
  w_2 = max_entries((alpha_hat_2 - as.matrix(R2))%*%t(as.matrix(C)) + b_matrix_2,
                    axis = 1)
  
  objective_2 = Minimize(sum(w_2))
  r2_loss_constraint_1 = (alpha_hat_2 - as.matrix(R2))*C[1] + b_vector[1] <= w_2
  r2_loss_constraint_2 = (alpha_hat_2 - as.matrix(R2))*C[2] + b_vector[2]  <= w_2
  r2_loss_constraint_3 = (alpha_hat_2 - as.matrix(R2))*C[3] + b_vector[3]  <= w_2
  r2_loss_constraint_4 = (alpha_hat_2 - as.matrix(R2))*C[4] + b_vector[4]  <= w_2
  
  prob_2 = Problem(objective_2)
  cpwl_solve_2 = solve(prob_2, solver = "SCS")  
  
  errors = c(cpwl_solve_1$getValue(alpha_hat_1) - R1, cpwl_solve_2$getValue(alpha_hat_2) - R2)
  losses = errors%*%t(as.matrix(C)) + B
  loss = sum(apply(losses, 1, max), na.rm = TRUE)
  names(loss) = "cpwl"
  return(loss)
}

# Implement RF algorithm 

#Bagging function
train_bagged_trees = function(X,y, l, B, loss){
  trees = vector(mode = "list", length = B)
  for(i in 1:B){
    
    index = sample(1:nrow(X), replace = TRUE)
    
    subset_X = matrix(NA, nrow = nrow(X), ncol = ncol(X)) 
    for(n in 1:length(index)){
      subset_X[n,] = X[index[n],]
    }
    subset_y = numeric(length(y))
    for(n in 1:length(index)){
      subset_y[n] = y[index[n]]
    }
    trees[[i]] = grow_tree(subset_X, subset_y, l, loss) 
  }
  return(trees)
}

#Grow trees with a subset of covariates. (Decorrelation of trees)

grow_tree_m = function(X, y, l, m, loss, region_min, a = NULL, b = NULL){
  checkmate::assert_matrix(X)
  checkmate::assert_numeric(y, len = nrow(X))
  checkmate::assert_int(l)
  m_covs  = m  # Dont want to get "m" as in the iteration tracker and "m" as in nr of covariates mixed up.# The old assignmetn template and my implemetations involve "m" as an interation tracker. 
  init = tree_split(X, y, l, loss, a, b)
  S_m <- list(init$R1, init$R2)
  results <- data.frame(j = init$j,
                        s = init$s,
                        R1_i = 2,
                        R2_i = 3,
                        gamma = NA)
  m = 1 
  M = 2
  while (m<=M) {
    
    covariate_indices = sample(1:ncol(X), m_covs, replace = FALSE)  
    truth_vector = logical(m_covs) 
    for(j in covariate_indices){
      truth_vector[j] = 
        isTRUE(length(which(X[S_m[[m]],j] <= median(X[S_m[[m]],j]))) >=l & length(which(X[S_m[[m]],j] > median(X[S_m[[m]],j]))) >=l)
    }
    if(any(truth_vector)){
      
      new_split = tree_split(as.matrix(X[S_m[[m]], covariate_indices]), 
                             y[S_m[[m]]], l, loss, a, b)
      
      # Get indices on total data form
      S_m[[M+1]] = S_m[[m]][new_split$R1]  
      S_m[[M+2]] = S_m[[m]][new_split$R2]
      
      #Append split variables and values to result matrix
      results[m+1,1] = covariate_indices[new_split$j] # if there were more covariates this would be necessary to select the correct variable
      results[m+1,2] = new_split$s
      
      # R1,R2 
      results[m+1,3] = M + 2
      results[m+1,4] = M + 3
      
      M = M + 2
    } else{results[m+1,5] = region_min(y[S_m[[m]]])}      
    m = m + 1
  }
  return(round(results,1))
}


#Train a random forest 
train_random_forest = function(X,y,l,B,m,loss, a = NULL, b = NULL){
  # Define the region estimate depending on the specified loss funciton. 
  if(as.character(substitute(loss)) == "linex_loss_rf"){region_min = function(R){
    return(log(length(R)/sum(exp(-a*R)))/a)}
  } 
  if(as.character(substitute(loss)) == "se_loss_rf"){region_min = function(R){return(mean(R))}}
  m_covs = m
  trees = vector(mode = "list", length = B)
  trees = foreach(i = 1:B, .export = c("grow_tree_m", "tree_split"), 
          .packages = c("checkmate", "uuml")) %dopar% {
    
    index = sample(1:nrow(X), replace = TRUE)
    
    subset_X = matrix(NA, nrow = nrow(X), ncol = ncol(X)) 
    for(n in 1:length(index)){
      subset_X[n,] = X[index[n],]
    }
    subset_y = numeric(length(y))
    for(n in 1:length(index)){
      subset_y[n] = y[index[n]]
    }
    grow_tree_m(subset_X, subset_y, l,m_covs, loss, region_min, a, b) 
  }
  return(trees)
}

#Predict with bagged trees 

predict_with_bagged_trees = function(new_data,trees){
  tmp = lapply(trees, function(x){predict_with_tree(new_data,x)})
  preds = numeric(nrow(new_data))
  for(i in 1:nrow(new_data)){
    preds[i] = mean(sapply(tmp, `[`, i))
  }
  return(preds)
}


X_rf_2021 = cbind(df_2021$normalized_irradiance, df_2021$cos_zenith_ahead)
y_rf_2021 = df_2021$hour_ahead

library(foreach)
library(doParallel)

cl = makeCluster(7)

registerDoParallel(cl)

start.time <- Sys.time()
rf_2021_linex = train_random_forest(X_rf_2021, y_rf_2021, 100, 200, 2, linex_loss_rf, a = a_irr, b_en)
end.time <- Sys.time()
time.taken <- end.time - start.time

stopCluster(cl)

# Increasing training data from 1 to two years increases the computation time by a factor of 4 
# 13 to 53 seconds for the sleected nmbr of cores and bootstrap iterations 
# Inreasing from 2 to 3 years increased by factor of 2 

library(dplyr)

X_rf_test = cbind(lulea_test$normalized_irradiance, lulea_test$cos_zenith_ahead)

nrow(X_rf_test)

length(irr_dir_rf_linex_pred)

irr_dir_rf_linex_pred = predict_with_bagged_trees(X_rf_test, rf_2021_linex)

en_dir_rf_linex_pred = predict_energy(irr_dir_rf_linex_pred, lulea_test, 
                                      lead(lulea_test$energy_produced), a = a_en, b = b_en)

en_dir_rf_linex_pred$Cost_Ratio 

hist_dir_rf_linex = en_dir_rf_linex_pred$Residuals %>%
  as.data.frame() %>% setNames(c("X")) %>% 
  ggplot(aes(x = X)) + 
  geom_histogram(aes(X), bins = 50) +
  labs(x = "Error in Wh/m^2", title = "Directly 
optimized") +
  theme(axis.text = element_text(size = 18), axis.title = element_text(size = 21), 
        plot.title = element_text(size = 25)) 



# SE regression tree model 

library(dplyr)

cl = makeCluster(7)

registerDoParallel(cl)

start.time <- Sys.time()
rf_2021_SE = train_random_forest(X_rf_2021, y_rf_2021, 100, 200, 2, se_loss_rf)
end.time <- Sys.time()
time.taken <- end.time - start.time

stopCluster(cl)

irr_unb_rf_pred = predict_with_bagged_trees(X_rf_test, rf_2021_SE)

en_unb_rf_pred = predict_energy(irr_unb_rf_pred, lulea_test, 
                                      lead(lulea_test$energy_produced), a = a_en, b = b_en)

en_unb_rf_pred$Cost_Ratio 

hist_unb_rf = en_unb_rf_pred$Residuals %>%
  as.data.frame() %>% setNames(c("X")) %>% 
  ggplot(aes(x = X)) + 
  geom_histogram(aes(X), bins = 50) +
  labs(x = "Error in Wh/m^2", title = "Least 
squares") +
  theme(axis.text = element_text(size = 18), axis.title = element_text(size = 21), 
        plot.title = element_text(size = 25))

# Indirectly optimized tree model 

# Input training data predictions from SE tree model and y 
 # we dont lead y because it is the "hour ahead" variable. 

se_tree_in_sample = predict_with_bagged_trees(X_rf_2021, rf_2021_SE)

fit_indirect_linex_tree = function(preds, y, a, b){
  train_resids = preds - y 
  M = length(train_resids)
  beta_linex = -(1/a)*log((1/M)*sum(exp(a*train_resids)))
  return(beta_linex)
}

beta_linex_tree = fit_indirect_linex_tree(preds, y_rf_2021, a_irr, b_en)

# Input out of sample predictions from unbiased model, and the bias
# This function just adds the bias and sets negative predictions to 0 

predict_irradiance_indirect_linex_tree = function(preds, bias){
  preds = preds + bias
  preds[which(preds < 0)] = 0
  return(preds)
}

linex_irr_ind_pred_tree = predict_irradiance_indirect_linex_tree(irr_unb_rf_pred,beta_linex_tree)

linex_en_ind_pred_tree = predict_energy(linex_irr_ind_pred_tree, lulea_test, 
                                   lead(lulea_test$energy_produced), 
                                   a_en, b_en)

linex_en_ind_pred_tree$Cost_Ratio


hist_ind_tree_linex = linex_en_ind_pred_tree$Residuals %>%
  as.data.frame() %>% setNames(c("X")) %>% 
  ggplot(aes(x = X)) + 
  geom_histogram(aes(X), bins = 50) +
  labs(x = "Error in Wh/m^2", title = "Indirectly 
optimized") +
  theme(axis.text = element_text(size = 18), axis.title = element_text(size = 21), 
        plot.title = element_text(size = 25))

library(cowplot)
plot_grid(hist_unb_rf, hist_ind_tree_linex, hist_dir_rf_linex, ncol = 3)










