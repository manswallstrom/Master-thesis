# Loss function shape/scale parameters

# The mean of this vector gives the approximate "exchange rate" 
#between solar irradiance and energy production. 
rate = mean(gbg_train_recent$energy_produced/gbg_train_recent$Global.Irradians..svenska.stationer., 
            na.rm = TRUE)

# A certain prediction in solar irradiance represents approximately
#which prediction in energy production? 
# predicted irradiance*rate tends to represent the prediction in energy production. So generally, 
# if the loss of prediction error in energy is 69/kWh = 0.069/Wh
#, then the loss of overprediction in irradiance is what? 
#(hat(x) - x) = epsilon = (rate*hat(x_irr) - rate*x_irr) == > hat(x_irr) - x_irr = epsilon/rate
# So we have loss(epsilon) = C2*epsilon = 0.069*epsilon (suppose loss is linlin).
# When training the irradiance model, the loss of epsilon/rate  (C_2_irr*epsilon/rate) 
#should be equivalent to 0.069*epsilon. Hence, C_2_irr/rate = 0.069 => C_2_irr = 0.069*rate 
# so, any linear coefficient for irradiance loss is simply 
#the coefficient for the corresponding segment of energy loss times the exchange rate. 

library(ggplot2)
library(dplyr)

# Rated power of french panel: 425. per square meter: 425/(2.005*1.042) = about 203.4
# 45 percent of rated power per sq m is 203.4*0.45 = 93.53 

# One sharp panel is 1.642088 square meters. (Should multiply all power produced by that..)
# Or simply: rated power per unit of area (m^2) is 152.2. Meaninng the 45 percent is 
# 68.49. So, the breakpoint on the negative side where we start having a steeper slope should be 
# an underestimate of 68.49. (now we are talking in wH/m^2 )
# Slope of the steeper negative segment? In fatemi it is -300/pu while the gentler one i s
# -50/pu. Maintain this ratio. 
#300/50 = 6. Hence, multiply the gentler slope with 6 
# for C3, maintaining the ratio to the price of energy means we multiply by 6/5 

# 20 percent of energy price. first slope 6x steeper than 2nd, 3rd 6/5 times steeper than 2nd
cpwl_slopes_en = c(-6*0.0005*0.2, -0.2*0.0005, 0.2*0.0005*(6/5), 0.069)

# Maintain same ratio between breakpoints => multiply negative breakpoint by 0.15/0.4 for 
# the positive breakpoint
D_cpwl_en = c(-91.53, 0, 91.53*(0.15/0.4))

intercepts_cpwl_en = CPWL_loss(cpwl_slopes_en, D_cpwl_en, 1, 1)$b

#Also maintain ratio between d1 and d3 (in terms of distance from 0)


plot_errors = seq(-180,50,0.01)

plot_loss = c(intercepts_cpwl_en[1] + cpwl_slopes_en[1]*plot_errors[which(plot_errors <= D_cpwl_en[1])], 
              cpwl_slopes_en[2]*plot_errors[which(plot_errors <= D_cpwl_en[2] & plot_errors > D_cpwl_en[1])], 
              cpwl_slopes_en[3]*plot_errors[which(plot_errors <= D_cpwl_en[3] & plot_errors > D_cpwl_en[2])], 
              intercepts_cpwl_en[4] + cpwl_slopes_en[4]*plot_errors[which(plot_errors > D_cpwl_en[3])])

max(gbg_train_recent$energy_produced)

a_en = 0.075
b_en = 0.0035

cbind(plot_errors, plot_loss) %>%
  as.data.frame() %>% setNames(c("X", "Y")) %>% 
  ggplot(aes(x = X, y = Y)) + 
  geom_line(aes(color = "CPWL"), size = 0.7) + 
  stat_function(fun=function(x) b_en*(exp(a_en*x)-a_en*x - 1), 
                aes(color = "LinEx"), size = 0.7)  + 
  ylim(c(0,0.1)) + 
  xlim(c(-180, 50)) + 
  labs(y = "Loss (SEK per Wh/m^2)", x = "Prediction error (Wh/m^2)") + 
  scale_color_manual(name='Loss function',
                     breaks=c('CPWL', 'LinEx'),
                     values=c('CPWL'='blue', 'LinEx'='red')) +
  theme(legend.title=element_text(size=17),
        legend.text=element_text(size=17), axis.title = element_text(size = 20), 
        axis.text = element_text(size = 14))


#"Exchange rate" for irradiance to energy production 

a_irr = a_en*rate

cpwl_slopes_irr = cpwl_slopes_en*rate
D_cpwl_irr = D_cpwl_en/rate
intercepts_cpwl_irr = intercepts_cpwl_en
