[![Build Status](https://travis-ci.org/pygeo/sense.svg?branch=master)](https://travis-ci.org/pygeo/sense)
[![Build Status](https://travis-ci.org/pygeo/sense.svg?branch=dev)](https://travis-ci.org/pygeo/sense)


# SenSE: Community SAR ScattEring model

SenSE is a generic community framework for radiative transfer (RT) modelling in the active microwave domain. It implements different existing models for scattering and emission for different surfaces in a coherent framework to simulate SAR backscattering coefficients as function of surface biogeophysical parameters. In the microwave domain the surface and canopy contribution of the total backscatter is usually estimated separately. Within the SenSE framework different model combination of surface and canopy models can be easily brought together and moreover analyzed. The analysis of the different model combination within one framework can be seen as the biggest advantage of the developed SenSE package. Currently implemented surface models are: Oh1992, Oh2004, Dubois95 and IEM and Water Cloud. Currently implemented canopy models are: SSRT and Water Cloud.

Within the MULTIPLY project the implemented rt-models will be used for the retrieval of biophysical parameters from the microwave remote sensing data. 
