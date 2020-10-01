from WeatheringModel import runWeatheringModel, ModelInputs

#First create an input object. Use this object to change the input parameters.
inpts = ModelInputs()
inpts.lum = 0.9 #set the incident flux to 90% modern Earth's

#Call the model and pass in the inputs. Calling without inpts will run the 
#modern Earth
[pCO2_out, pH_out, Ts_out], status = runWeatheringModel(inputs=inpts)

#to run the modern Earth just uncomment the below line
#[pCO2_out, pH_out, Ts_out], status = runWeatheringModel()


#print the model outputs
print("pH=%0.2e"%(pH_out))
print("pCO2=%2.3e bar"%(pCO2_out))
print("Surface temperature=%0.2f K"%(Ts_out))
