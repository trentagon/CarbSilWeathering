# CarbSilWeathering

The CarbSilWeathering code will simulate the carbonate-silicate weathering
cycle on an Earth-like planet in the habitable zone (HZ). The model is
described in detail in:
* [Krissansen-Totton & Catling (2017)](https://doi.org/10.1038/ncomms15423)
* [Krissansen-Totton *et al.* (2018)](https://doi.org/10.1073/pnas.1721296115)
* Lehmer et al. (in prep)


The current version of the model, contained in the WeatheringModel directory, 
is based on the model description in Lehmer *et al.* (in prep).


### How to run the model
Import the runWeatheringModel function and the ModelInputs class from the
WeatheringModel.py file. The model can be run by simply calling
runWeatheringModel(), which will calculate the steady-state atmospheric carbon
dioxide level, surface temperature, and ocean pH. An example of how to use the
model and how to change the model inputs is shown in the example.py file.


#### Need help?
Please email Owen Lehmer at owen.r.lehmer@nasa.gov or Joshua Krissansen-Totton
at jkt@ucsc.edu for help running the model.

#### Funding
This model has been developed through funding from:
* NASA Exobiology Program grant NNX15AL23G
* NASA Astrobiology Institute’s Virtual Planetary Laboratory, grant NNA13AA93A
* NASA Earth and Space Science Fellowship program, grant NNX15AR63H
* Simons Collaboration on the Origin of Life Award 511570
* NASA's Virtual Planetary Laboratory, grant 80NSSC18K0829
* NASA Hubble Fellowship grant HF2-51437

