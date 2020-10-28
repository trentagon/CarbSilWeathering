from weatheringmodel import runWeatheringModel, ModelInputs

#First create an input object. Use this object to change the input parameters.
inpts = ModelInputs() #default values are for modern Earth
inpts.lum = 0.9 #set the incident flux to 90% modern Earth's
guess = [3.301E-3,3.160E-3,8.242E-4]

#Call the model and pass in the inputs. 
times, Co, Ao, pCO2, pH, Ts = runWeatheringModel(inpts)


#print the model outputs
print("pH=%0.2e"%(pH[-1]))
print("pCO2=%2.3e bar"%(pCO2[-1]))
print("Surface temperature=%0.2f K"%(Ts[-1]))
