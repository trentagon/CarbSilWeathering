from weatheringmodel import runWeatheringModel, ModelInputs
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import time

#NOTE - running the model for 4E9 years (as done below) takes about ~15 minutes

#First create an input object. Use this object to change the input parameters.
inpts = ModelInputs()

#Configure the change in the model with time. The below equations are from the
#Krissansen-Totton et al. (2018) paper.

#biology with time
def f_bio_func(t, Bp):
    fb_EE=1/(1-Bp)
    return 1-1/(fb_EE + np.exp(-10*(t-0.6e9)/1e9))
inpts.getAtTime_f_bio = lambda t: f_bio_func(t, 0.6) #0.6 is the Cambrian time

#land fraction
def land_frac (t0, lfrac, growth_timing):
    LL=1/(1-lfrac)
    land_frac=np.max([0.0,1-1/(LL+np.exp(-10*(t0/1e9-growth_timing)))])
    return land_frac
f_land_func = lambda t: land_frac(t, 0.5, 2.3)
inpts.getAtTime_f_land = f_land_func

#internal heat
def internal_heat(t, n_out):
    return (1 - t/4.5E9)**(-n_out)
inpts.getAtTime_Q = lambda t: internal_heat(t, 0.5)

#incident flux
inpts.getAtTime_lum = lambda t: 1/(1+0.4*(t/4.6e9))


def printRuntime(dur):
    """
    Helper function to print the runtime
    """
    hours = dur//3600
    dur -= hours*3600
    minutes = dur//60
    dur -= minutes*60
    print("runtime %02d:%02d:%02.3f"%(hours, minutes, dur))

#run the model and track the total runtime
start = time.time()
times, Co, Ao, pCO2, pH, Ts = runWeatheringModel(inputs=inpts, endtime=4E9,
        steadystate=False)
end = time.time()
printRuntime(end - start) #print the time it took

#plot the results
fig, (ax0, ax1, ax2) = plt.subplots(3, 1, sharex=True, figsize=(5,6))
plt.subplots_adjust(hspace=0.05, left=0.16)

times = times/1E9

ax0.plot(times, pCO2)
ax0.set_ylabel(r"pCO$_{2}$ [bar]")
ax0.set_yscale('log')

ax1.plot(times, pH)
ax1.set_ylabel("pH")
ax1.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

ax2.plot(times, Ts)
ax2.set_ylabel("Sur. Temp. [K]")
ax2.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))

ax2.set_xlabel("Time Before Present [Ga]")

plt.show()
