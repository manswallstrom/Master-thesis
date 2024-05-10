# Regression tree algorithm 

tree_split <- function(X, y, l, loss, a = NULL, b = NULL, C = NULL, D = NULL, 
                       intercepts = NULL){
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
        
        losses[k, j] <- loss(y[R1], y[R2], a, b, C, D, intercepts)
        
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
                        inner_min = NA)
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
    terminal_node = isFALSE(is.na(tree$inner_min[row]))
  }
  return(tree$inner_min[row])
} 

predict_with_tree = function(new_data,tree){
  checkmate::assert_matrix(new_data)
  predictions <- apply(new_data, 1, function(x){pred(x, tree)})
  return(predictions)
}


# Loss functions in RF split format 

# Note, the way Ive coded it, all loss functions need to take the parameters of all possible loss
# functions as an argument. I just set the default to NULL. 
# (This is why linex loss takes CPWL parameters
# as one of its arguments etc..) 

# Linex loss in RF split format 
linex_loss_rf = function(R1, R2, a, b, C = NULL, D = NULL, intercepts = NULL){
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

se_loss_rf = function(R1, R2, a = NULL, b = NULL, C = NULL, D = NULL, intercepts = NULL){
  within_region_min_1 = mean(R1) 
  within_region_min_2 = mean(R2)
  errors = c(within_region_min_1 - R1, within_region_min_2 - R2)
  loss = sum(errors^2)
  names(loss) = "se"
  return(loss)
}

#CPWL loss 

cpwl_region_min = function(C, D, preds, y){
  residuals =  preds - y 
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

cpwl_loss_rf = function(R1, R2, a = NULL, b = NULL, C, D, intercepts){
  # Within-region loss minimizing value is an intercept-only regression fit to the  
  # loss function in question. 
  # This is equivalent to an intercept-only least squares regression (the mean)
  # plus the indirect CPWL bias for that model. 
  within_region_min_1 = cpwl_region_min(C, D, 0, R1)
  within_region_min_2 = cpwl_region_min(C, D, 0, R2)
  
  b_matrix = cbind(rep(intercepts[1], length(R1) + length(R2)), 
                   rep(intercepts[2], length(R1) + length(R2)), 
                   rep(intercepts[3], length(R1) + length(R2)), 
                   rep(intercepts[4], length(R1) + length(R2)))
  errors = c(within_region_min_1 - R1, within_region_min_2 - R2)
  losses = errors%*%t(as.matrix(C)) + b_matrix
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

grow_tree_m = function(X, y, l, m, loss, region_min, a = NULL, b = NULL, C = NULL, 
                       D = NULL, intercepts = NULL){
  checkmate::assert_matrix(X)
  checkmate::assert_numeric(y, len = nrow(X))
  checkmate::assert_int(l)
  m_covs  = m  # Dont want to get "m" as in the iteration tracker and "m" as in nr of covariates mixed up.
  init = tree_split(X, y, l, loss, a, b, C, D, intercepts)
  S_m <- list(init$R1, init$R2)
  results <- data.frame(j = init$j,
                        s = init$s,
                        R1_i = 2,
                        R2_i = 3,
                        inner_min = NA)
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
                             y[S_m[[m]]], l, loss, a, b, C, D, intercepts)
      
      # Get indices on total data form
      S_m[[M+1]] = S_m[[m]][new_split$R1]  
      S_m[[M+2]] = S_m[[m]][new_split$R2]
      
      #Append split variables and values to result matrix
      results[m+1,1] = covariate_indices[new_split$j] 
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
train_random_forest = function(X,y,l,B,m,loss, a = NULL, b = NULL, C = NULL, 
                               D = NULL, intercepts = NULL){
  # Define the region estimate depending on the specified loss function. 
  if(as.character(substitute(loss)) == "linex_loss_rf"){region_min = function(R){
    return(log(length(R)/sum(exp(-a*R)))/a)}
  } 
  if(as.character(substitute(loss)) == "se_loss_rf"){region_min = function(R){return(mean(R))}}
  if(as.character(substitute(loss)) == "cpwl_loss_rf"){region_min = function(R){
    return(cpwl_region_min(C, D, 0, R))}}
  m_covs = m
  trees = vector(mode = "list", length = B)
  trees = foreach(i = 1:B, .export = c("grow_tree_m", "tree_split", "cpwl_region_min"), 
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
    grow_tree_m(subset_X, subset_y, l,m_covs, loss, region_min, a, b, C, D, intercepts) 
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

fit_indirect_linex_tree = function(preds, y, a, b){
  train_resids = preds - y 
  M = length(train_resids)
  beta_linex = -(1/a)*log((1/M)*sum(exp(a*train_resids)))
  return(beta_linex)
}

# Input out of sample predictions from unbiased model, and the bias
# This function just adds the bias and sets negative predictions to 0 

predict_irradiance_indirect_linex_tree = function(preds, bias){
  preds = preds + bias
  preds[which(preds < 0)] = 0
  return(preds)
}
