sol_lulea = read.csv("sol_lulea.csv", sep = ";")
View(sol_lulea)

#Calculate the monthly means from this. To use when transforming solar irradiance to 
#generated power. 
vind_lulea = read.csv("vind_lulea.csv", sep = ";")
View(vind_lulea)

temp_lulea = read.csv("temp_lulea.csv", sep = ";")
View(temp_lulea)


#Calculate solar zenith angle 

t = paste(sol_lulea$Datum, 
          sol_lulea$Tid..UTC.)

#The time-date vector but in UTC (delayed 1 hour)
t_UTC = c("1982-12-31 23:00:00", t[1:(length(t)-1)])

library(oce)

elevation_refr = sunAngle(t_UTC, longitude = 22.1113, latitude = 65.5436, 
                                    useRefraction = TRUE)$altitude


#Cos function calculates cosine in radians. First transform degrees into radians. 

sol_lulea$cos_zenith = cos((90 - elevation_refr)*(pi/180))

# Set negative coisne angles to 0. 
#Because otherwise we get negative weighted irradiance which doesnt make sense

sol_lulea$Global.Irradians..svenska.stationer.[which(sol_lulea$cos_zenith < 0)] = NA

sol_lulea$cos_zenith[which(sol_lulea$cos_zenith < 0 )] = NA

### Set the irradiance where cos of zenith angle is < 0 
#to NA. Otherwise we get division by 0. 

sol_lulea$weighted_irradiance = sol_lulea$Global.Irradians..svenska.stationer./
  sol_lulea$cos_zenith

## Explore the properties of cos zenith weighted time series 

#Would need to visualize density of values to make anything of this pllot
plot(sol_lulea$cos_zenith, sol_lulea$Global.Irradians..svenska.stationer.)


# cosine zenith weighting appears to remove seasonality. There are some very large values whne 
#cosine weighting (because we divide by 0.0042 and stuff like that) 
plot(ts(sol_lulea$Global.Irradians..svenska.stationer.))
plot(ts(sol_lulea$weighted_irradiance))

###
                                
library(dplyr)

#Create lagged normalized variables 

sol_lulea$normalized_irradiance = sol_lulea$weighted_irradiance*
  lead(sol_lulea$cos_zenith)  

sol_lulea$lag_1_normalized_irradiance = lag(sol_lulea$weighted_irradiance, n = 1)*
  lead(sol_lulea$cos_zenith, n = 1) 

sol_lulea$lag_2_normalized_irradiance = lag(sol_lulea$weighted_irradiance, n = 2)*
  lead(sol_lulea$cos_zenith, n = 1) 

sol_lulea$hour_ahead = lead(sol_lulea$Global.Irradians..svenska.stationer.)

sol_lulea$cos_zenith_ahead = lead(sol_lulea$cos_zenith)

### Estimate produced electricity using monthly means of wind and temperature 

#Calculate monthly mean of wind speed (1990 - 2020)
# We always want to use only modern values to convert irradiance into energy. Because 
# we are forecasting now, not 20 years ago. Hence, wind and temperature should be based on recent 
# values, (also the a parameter in linex) regardless of whether we are training the 
#irradiance model on only recent values or on the whole data set.  
# Usually 30 year periods are used for this purpose. Here I use the average between 1990 and 
# 2020. 

View(vind_lulea)

monthly_means = numeric(12)
names(monthly_means) = month.name
month_indices = c("01","02","03","04","05","06","07","08","09","10","11","12")

# 1990-2020  
vind_lulea_90_20 = vind_lulea[which(vind_lulea$Datum>1990 & vind_lulea$Datum < 2020), ]
split = strsplit(vind_lulea_90_20$Datum, "-")

for(i in 1:12){
  tmp = unlist(lapply(split, function(x){x[2] == month_indices[i]}))
  which_obs = which(tmp)
  monthly_means[i] = mean(vind_lulea_90_20[which_obs, "Vindhastighet"], na.rm = TRUE)
}

#Calculate monthly mean of temperature (1990-2020)

monthly_means_temp = numeric(12)
names(monthly_means_temp) = month.name

#From 1990-2020  
temp_lulea_90_20 = temp_lulea[which(temp_lulea$Från.Datum.Tid..UTC. > 1990 & 
                                      temp_lulea$Från.Datum.Tid..UTC. < 2020), ]

split_temp = strsplit(temp_lulea_90_20$Från.Datum.Tid..UTC., "-")

for(i in 1:12){
  tmp = unlist(lapply(split_temp, function(x){x[2] == month_indices[i]}))
  which_obs = which(tmp)
  monthly_means_temp[i] = mean(temp_lulea_90_20[which_obs, "Lufttemperatur"], na.rm = TRUE)
}

# Put monthly mean wind speed and temperature in data 

sol_lulea$ambient_temp = NA
sol_lulea$ambient_wind = NA
split_sol = strsplit(sol_lulea$Datum, "-")
for(i in 1:12){
  tmp = unlist(lapply(split_sol, function(x){x[2] == month_indices[i]}))
  which_obs = which(tmp)
  sol_lulea[which_obs, "ambient_temp"] = monthly_means_temp[i]
  sol_lulea[which_obs, "ambient_wind"] = monthly_means[i]
}


#Calculate module temperature and temperature efficiency for all observations 

# Note, I may need to change the coefficient which is set to 6.62 to something more 
#location specific. (But they used the same for all locations in the Fatemi paper)

sol_lulea$module_temperature = sol_lulea$ambient_temp + 
  sol_lulea$Global.Irradians..svenska.stationer.*(47.5 - 20)/
  (800 + 6.62*(sol_lulea$ambient_wind - 1)*(47.5 - 20))  
sol_lulea$temperature_efficiency = 1 - 0.0044*(sol_lulea$module_temperature - 25) 

sol_lulea$energy_produced = 0.93*0.152*sol_lulea$temperature_efficiency*sol_lulea$Global.Irradians..svenska.stationer.



#adf test
library(tseries)

#NA:s.. 
adf.test(sol_lulea$weighted_irradiance)

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

predict_energy = function(irradiation_preds, weather, newy, a, b){
  #Predict module temperature 
  module_temperature_pred = weather[, "ambient_temp"] +    #NOTE, Hard coded so that it needs thsi specific column name.
    irradiation_preds*(47.5 - 20)/                              
    (800 + 6.62*(weather[, "ambient_wind"] - 1)*(47.5 - 20))
  #Predict temperature-efficienct coefficient
  temperature_efficiency_pred = 
    1 - 0.0044*(module_temperature_pred - 25)
  #Predict energy production
  energy_produced_pred = 
    0.93*0.152*temperature_efficiency_pred*irradiation_preds
  
  resids = energy_produced_pred - newy
  LinEx = linex_loss(energy_produced_pred, newy, a, b)
  cost_ratio = linex_loss(energy_produced_pred, newy, a, b)/
    linex_loss(0, newy, a, b)
  
  return(list(Predictions = energy_produced_pred, Residuals = resids, 
              LinEx = LinEx, Cost_Ratio = cost_ratio))
}


# Create a subset with only brght hours (cosine zenith angle during the time forecasted 
#should be greater than 0.1)


sol_lulea_bright = sol_lulea[which(sol_lulea$cos_zenith_ahead > 0.15), ]
rownames(sol_lulea_bright) = 1:nrow(sol_lulea_bright)

View(sol_lulea_bright)

# Train on  - september 2021, test on rest 
#full data set 
lulea_train_full = sol_lulea_bright[1:(which(sol_lulea_bright$Datum == "2021-10-01")[1]-1), ]
#Only from 2010 onwards 
lulea_train_recent = lulea_train_full[which(lulea_train_full$Datum > 2010), ]
rownames(lulea_train_recent) = 1:nrow(lulea_train_recent)

lulea_test = sol_lulea_bright[which(sol_lulea_bright$Datum == "2021-10-01")[1]:nrow(sol_lulea_bright), ]
rownames(lulea_test) = 1:nrow(lulea_test)

# one-year df 

df_2021 = lulea_train_recent[which(lulea_train_recent$Datum >= "2020-10-01"), ]

