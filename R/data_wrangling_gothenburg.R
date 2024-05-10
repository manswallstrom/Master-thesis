sol_gbg = read.csv("gbg_sol.csv", sep = ";")
View(sol_gbg)

#Calculate the monthly means from this. To use when transforming solar irradiance to 
#generated power. 
vind_gbg = read.csv("gbg_vind.csv", sep = ";")
View(vind_gbg)

temp_gbg = read.csv("gbg_temp.csv", sep = ";")
View(temp_gbg)


#Calculate solar zenith angle 

t = paste(sol_gbg$Datum, 
          sol_gbg$Tid..UTC.)

#The time-date vector but in UTC (delayed 1 hour)
t_UTC = c("1982-12-31 23:00:00", t[1:(length(t)-1)])

library(oce)

elevation_refr = sunAngle(t_UTC, longitude = 11.9796, latitude = 57.6879, 
                          useRefraction = TRUE)$altitude


#Cos function calculates cosine in radians. First transform degrees into radians. 

sol_gbg$cos_zenith = cos((90 - elevation_refr)*(pi/180))

# Set negative coisne angles to 0. 
#Because otherwise we get negative weighted irradiance which doesnt make sense

sol_gbg$Global.Irradians..svenska.stationer.[which(sol_gbg$cos_zenith < 0)] = NA

sol_gbg$cos_zenith[which(sol_gbg$cos_zenith < 0 )] = NA

### Set the irradiance where cos of zenith angle is < 0 
#to NA. Otherwise we get division by 0. 

sol_gbg$weighted_irradiance = sol_gbg$Global.Irradians..svenska.stationer./
  sol_gbg$cos_zenith

## Explore the properties of cos zenith weighted time series 

#Would need to visualize density of values to make anything of this pllot

# cosine zenith weighting appears to remove seasonality. There are some very large values whne 
#cosine weighting (because we divide by 0.0042 and stuff like that) 

###

library(dplyr)

#Create lagged normalized variables 

sol_gbg$normalized_irradiance = sol_gbg$weighted_irradiance*
  lead(sol_gbg$cos_zenith)  

sol_gbg$lag_1_normalized_irradiance = lag(sol_gbg$weighted_irradiance, n = 1)*
  lead(sol_gbg$cos_zenith, n = 1) 

sol_gbg$lag_2_normalized_irradiance = lag(sol_gbg$weighted_irradiance, n = 2)*
  lead(sol_gbg$cos_zenith, n = 1) 

sol_gbg$hour_ahead = lead(sol_gbg$Global.Irradians..svenska.stationer.)

sol_gbg$cos_zenith_ahead = lead(sol_gbg$cos_zenith)

### Estimate produced electricity using monthly means of wind and temperature 

#Calculate monthly mean of wind speed (1990 - 2020)
# We always want to use only modern values to convert irradiance into energy. Because 
# we are forecasting now, not 20 years ago. Hence, wind and temperature should be based on recent 
# values, (also the a parameter in linex) regardless of whether we are training the 
#irradiance model on only recent values or on the whole data set.  
# Usually 30 year periods are used for this purpose. Here I use the average between 1990 and 
# 2020. 

View(vind_gbg)

monthly_means = numeric(12)
names(monthly_means) = month.name
month_indices = c("01","02","03","04","05","06","07","08","09","10","11","12")

# 1990-2020  
vind_gbg_90_20 = vind_gbg[which(vind_gbg$Datum>1990 & vind_gbg$Datum < 2020), ]
split = strsplit(vind_gbg_90_20$Datum, "-")

for(i in 1:12){
  tmp = unlist(lapply(split, function(x){x[2] == month_indices[i]}))
  which_obs = which(tmp)
  monthly_means[i] = mean(vind_gbg_90_20[which_obs, "Vindhastighet"], na.rm = TRUE)
}

#Calculate monthly mean of temperature (1990-2020)

monthly_means_temp = numeric(12)
names(monthly_means_temp) = month.name

#From 1990-2020  
temp_gbg_90_20 = temp_gbg[which(temp_gbg$Från.Datum.Tid..UTC. > 1990 & 
                                      temp_gbg$Från.Datum.Tid..UTC. < 2020), ]

split_temp = strsplit(temp_gbg_90_20$Från.Datum.Tid..UTC., "-")

for(i in 1:12){
  tmp = unlist(lapply(split_temp, function(x){x[2] == month_indices[i]}))
  which_obs = which(tmp)
  monthly_means_temp[i] = mean(temp_gbg_90_20[which_obs, "Lufttemperatur"], na.rm = TRUE)
}

# Put monthly mean wind speed and temperature in data 

sol_gbg$ambient_temp = NA
sol_gbg$ambient_wind = NA
split_sol = strsplit(sol_gbg$Datum, "-")
for(i in 1:12){
  tmp = unlist(lapply(split_sol, function(x){x[2] == month_indices[i]}))
  which_obs = which(tmp)
  sol_gbg[which_obs, "ambient_temp"] = monthly_means_temp[i]
  sol_gbg[which_obs, "ambient_wind"] = monthly_means[i]
}


#Calculate module temperature and temperature efficiency for all observations 

# Note, I may need to change the coefficient which is set to 6.62 to something more 
#location specific. (But they used the same for all locations in the Fatemi paper)

sol_gbg$module_temperature = sol_gbg$ambient_temp + 
  sol_gbg$Global.Irradians..svenska.stationer.*(45 - 20)/
  (800 + 6.62*(sol_gbg$ambient_wind - 1)*(45 - 20))  
sol_gbg$temperature_efficiency = 1 - 0.00353*(sol_gbg$module_temperature - 25) 

sol_gbg$energy_produced = 0.93*0.2034*sol_gbg$temperature_efficiency*sol_gbg$Global.Irradians..svenska.stationer.



#adf test
# library(tseries)

#NA:s.. 
# adf.test(sol_gbg$weighted_irradiance)

#Which series should i test for unit root? We are regressing solar irradiance on 
#lagged versions of solar irradiance which are normalized by cosine weighting so 
# can I perform normal unit root tests to test for stationarity? 


#This function takes a gradient descent object from the above function and predicts irrdadiation.
#Returns preictions, residuals and linex cost. 
#Newx should be a model matrix of the same kind as in the training of the gd 

predict_irradiation = function(gd, newx, newy, a, b){
  parameters = gd[nrow(gd), -c(1:2)]
  # for(j in 2:ncol(newx)){
  #  newx[,j] = standardize(newx[,j])
  #  }
  preds = as.matrix(newx)%*%t(as.matrix(parameters))
  preds[which(preds < 0)] = 0  # Set negative predictions to 0 
  return(list(Predictions = preds, Residuals = preds - newy, 
              LinEx = linex_loss(preds, newy, a, b))) 
}

### Predictions of energy production 

#Predicted module temperature, efficiency (one ahead predictions )

predict_energy = function(irradiation_preds, weather, newy, a = NULL, b =  NULL, C = NULL, 
                            D = NULL, intercepts = NULL){
  #Predict module temperature 
  module_temperature_pred = weather[, "ambient_temp"] +    #NOTE, Hard coded so that it needs thsi specific column name.
    irradiation_preds*(45 - 20)/                              
    (800 + 6.62*(weather[, "ambient_wind"] - 1)*(45 - 20))
  #Predict temperature-efficienct coefficient
  temperature_efficiency_pred = 
    1 - 0.00353*(module_temperature_pred - 25)
  #Predict energy production
  energy_produced_pred = 
    0.93*0.2034*temperature_efficiency_pred*irradiation_preds
  
  output = list(Predictions = NULL, Residuals = NULL, 
                LinEx = NULL, CR_LinEx = NULL, CR_CPWL = NULL)
  output$Residuals = energy_produced_pred - newy
  if(!is.null(a) & !is.null(b)){output$LinEx = linex_loss(energy_produced_pred, newy, a, b)
  output$cost_ratio_linex = linex_loss(energy_produced_pred, newy, a, b)/
    linex_loss(0, newy, a, b)}
  if(!is.null(C) & !is.null(D) & !is.null(intercepts)){
    output$CPWL = CPWL_loss(C, D, energy_produced_pred, newy)$loss
    output$cost_ratio_cpwl = sum(CPWL_loss(C, D, energy_produced_pred, newy)$loss, na.rm = TRUE)/
    sum(CPWL_loss(C, D, 0, newy)$loss, na.rm = TRUE)}
  return(output)
}

# Create a subset with only brght hours (cosine zenith angle during the time forecasted 
#should be greater than 0.1)


sol_gbg_bright = sol_gbg[which(sol_gbg$cos_zenith_ahead > 0.15), ]
rownames(sol_gbg_bright) = 1:nrow(sol_gbg_bright)

View(sol_gbg_bright)

# Train on  - september 2021, test on rest 
#full data set 
gbg_train_full = sol_gbg_bright[1:(which(sol_gbg_bright$Datum == "2021-10-01")[1]-1), ]
#Only from 2010 onwards 
gbg_train_recent = gbg_train_full[which(gbg_train_full$Datum > 2010), ]
rownames(gbg_train_recent) = 1:nrow(gbg_train_recent)

gbg_test = sol_gbg_bright[which(sol_gbg_bright$Datum == "2021-10-01")[1]:
                            (which(sol_gbg_bright$Datum == "2023-10-01")[1] - 1), ]
rownames(gbg_test) = 1:nrow(gbg_test)

# Validation df is same as df_2021. 

# one-year df 

df_2021_gbg = gbg_train_recent[which(gbg_train_recent$Datum >= "2020-10-01"), ]

# Training set for the grid search 

gbg_train_3 = gbg_train_full[which(gbg_train_full$Datum < "2020-10-01" & 
                                    gbg_train_full$Datum >= "2017-10-01"), ]
rownames(gbg_train_3) = 1:nrow(gbg_train_3)

# Final trianing set 

gbg_train = gbg_train_full[which(gbg_train_full$Datum >= "2018-10-01"), ]
rownames(gbg_train) = 1:nrow(gbg_train)



