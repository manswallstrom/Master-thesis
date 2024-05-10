# gbg regression tree LinEx grid search

X_rf_3 = cbind(gbg_train_3$normalized_irradiance, gbg_train_3$cos_zenith_ahead)
y_rf_3 = gbg_train_3$hour_ahead

View(gbg_train_3)
View(df_2021_gbg)

X_rf_val = cbind(df_2021_gbg$normalized_irradiance, df_2021_gbg$cos_zenith_ahead)
y_rf_val = df_2021_gbg$hour_ahead

# Grid search for linex tree : Try minimum leaf sizes 25, 50, 100, 250, 500 

# Store trees in a list of lists. 

linex_grid_search_trees = vector(mode = "list", length = 6)
names(linex_grid_search_trees) = c("L_25", "L_50", "L_100", "L_250", "L_500", "L_1000")
linex_grid_search_errors = vector(mode = "list", length = 6)
names(linex_grid_search_errors) = c("L_25", "L_50", "L_100", "L_250", "L_500", "L_1000")
linex_grid_search_CR = vector(mode = "list", length = 6)
names(linex_grid_search_CR) = c("L_25", "L_50", "L_100", "L_250", "L_500", "L_1000")

library(foreach)
library(doParallel)
library(dplyr)

grid = c(25, 50, 100, 250, 500, 1000)

start.time <- Sys.time()
for(i in 1:length(grid)){
  
  cl = makeCluster(7)
  
  registerDoParallel(cl)
  
  rf_linex = train_random_forest(X_rf_3, y_rf_3, grid[i], 100, 2, linex_loss_rf, 
                                a = a_irr, b = b_en)
  
  stopCluster(cl)
  
  irr_pred = predict_with_bagged_trees(X_rf_val, rf_linex)
  en_pred = predict_energy(irr_pred, df_2021_gbg, lead(df_2021_gbg$energy_produced), 
                           a = a_en, b = b_en)
  
  linex_grid_search_trees[[i]] = rf_linex
  linex_grid_search_errors[[i]] = en_pred$Residuals
  linex_grid_search_CR[[i]] = en_pred$cost_ratio_linex
  print(linex_grid_search_CR[[i]])
}

end.time = Sys.time()
end.time - start.time 


library(tidyverse)

as.data.frame(
  cbind(linex_grid_search_errors[[1]], linex_grid_search_errors[[2]], linex_grid_search_errors[[3]], 
        linex_grid_search_errors[[4]], linex_grid_search_errors[[5]], linex_grid_search_errors[[6]])) %>%
  setNames(c("L = 25", "L = 50", "L = 100", "L = 250", "L = 500", "L = 1000")) %>% 
  pivot_longer(cols = 1:6, names_to = "Leaf size", values_to = "Error") %>% 
  ggplot(aes(x = Error)) + 
  geom_histogram(bins = 50)  + 
  facet_wrap(.~factor(`Leaf size`, levels = c("L = 25", "L = 50", "L = 100", 
                                              "L = 250", "L = 500", "L = 1000"))) + 
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 20), 
        strip.text = element_text(size = 14))
        
library(xtable)
xtable(cbind(linex_grid_search_CR$L_25, linex_grid_search_CR$L_50, linex_grid_search_CR$L_100, 
                     linex_grid_search_CR$L_250, linex_grid_search_CR$L_500, linex_grid_search_CR$L_1000), 
       digits = 4)


