# CarbSilWeathering

The CarbSilWeathering code will simulate the carbonate-silicate weathering
cycle on an Earth-like planet in the habitable zone (HZ). The model is
described in detail in:
[Krissansen-Totton & Catling (2017)](doi.org/10.1038/ncomms15423)
[Krissansen-Totton *et al.* (2018)](doi.org/10.1073/pnas.1721296115)
Lehmer et al. (in prep)

This file contains all of the functions necessary to simulate pCO2 in the 
atmosphere of a planet with liquid water at the surface. The model is written
using the equations of "Constraining the climate and ocean pH of the early
Earth with a geological carbon cycle model" by Krissansen-Totten et al. 
(2018). Unless otherwise noted, all equation references in this code refer 
to the equations of that paper, which is abbreviated as JKT in this code.

This code was written by Owen Lehmer, questions can be addressed to:
    owen.r.lehmer@nasa.gov

                           RUNNING THE MODEL
Steps: 
Simply call runWeatheringModel() to run the model. See the end of this file 
for an example call to runWeatheringModel(). Uncomment that example and run
this file to see the model in action. It's also recommended that you set DEBUG
to True for additional model output.

