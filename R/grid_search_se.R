# gbg regression tree least squares grid search

X_rf_3 = cbind(gbg_train_3$normalized_irradiance, gbg_train_3$cos_zenith_ahead)
y_rf_3 = gbg_train_3$hour_ahead

X_rf_val = cbind(df_2021_gbg$normalized_irradiance, df_2021_gbg$cos_zenith_ahead)
y_rf_val = df_2021_gbg$hour_ahead

# Grid search for se tree : Try minimum leaf sizes 25, 50, 100, 250, 500 

# Store trees in a list of lists. 

se_grid_search_trees = vector(mode = "list", length = 6)
names(se_grid_search_trees) = c("L_25", "L_50", "L_100", "L_250", "L_500", "L_1000")
se_grid_search_errors = vector(mode = "list", length = 6)
names(se_grid_search_errors) = c("L_25", "L_50", "L_100", "L_250", "L_500", "L_1000")
se_grid_search_rmse = vector(mode = "list", length = 6)
names(se_grid_search_rmse) = c("L_25", "L_50", "L_100", "L_250", "L_500", "L_1000")

library(foreach)
library(doParallel)
library(dplyr)

grid = c(25, 50, 100, 250, 500, 1000)

start.time <- Sys.time()
for(i in 1:length(grid)){
  
  cl = makeCluster(7)
  
  registerDoParallel(cl)
  
  rf_se = train_random_forest(X_rf_3, y_rf_3, grid[i], 100, 2, se_loss_rf, 
                                 a = a_irr, b = b_en)
  
  stopCluster(cl)
  
  irr_pred = predict_with_bagged_trees(X_rf_val, rf_se)
  en_pred = predict_energy(irr_pred, df_2021_gbg, lead(df_2021_gbg$energy_produced), 
                           a = a_en, b = b_en)
  
  se_grid_search_trees[[i]] = rf_se
  se_grid_search_errors[[i]] = en_pred$Residuals
  se_grid_search_rmse[[i]] = sqrt(mean(en_pred$Residuals^2, na.rm = TRUE))
  print(se_grid_search_rmse[[i]])
}

end.time = Sys.time()
end.time - start.time 


library(tidyverse)

as.data.frame(
  cbind(se_grid_search_errors[[1]], se_grid_search_errors[[2]], se_grid_search_errors[[3]], 
        se_grid_search_errors[[4]], se_grid_search_errors[[5]], se_grid_search_errors[[6]])) %>%
  setNames(c("L = 25", "L = 50", "L = 100", "L = 250", "L = 500", "L = 1000")) %>% 
  pivot_longer(cols = 1:6, names_to = "Leaf size", values_to = "Error") %>% 
  ggplot(aes(x = Error)) + 
  geom_histogram(bins = 50)  + 
  facet_wrap(.~factor(`Leaf size`, levels = c("L = 25", "L = 50", "L = 100", 
                                              "L = 250", "L = 500", "L = 1000"))) + 
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 20), 
        strip.text = element_text(size = 14))

library(xtable)
xtable(cbind(se_grid_search_rmse$L_25, se_grid_search_rmse$L_50, se_grid_search_rmse$L_100, 
             se_grid_search_rmse$L_250, se_grid_search_rmse$L_500, se_grid_search_rmse$L_1000), 
       digits = 4)


