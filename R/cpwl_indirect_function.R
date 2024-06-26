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

