[![Build Status](https://travis-ci.org/pygeo/sense.svg?branch=master)](https://travis-ci.org/pygeo/sense)
[![Build Status](https://travis-ci.org/pygeo/sense.svg?branch=dev)](https://travis-ci.org/pygeo/sense)


# SenSE: Community SAR ScattEring model

SenSE is a generic community framework for radiative transfer (RT) modelling in the active microwave domain. It implements different existing models for scattering and emission for different surfaces in a coherent framework to simulate SAR backscattering coefficients as function of surface biogeophysical parameters. In the microwave domain the surface and canopy contribution of the total backscatter is usually estimated separately. Within the SenSE framework different model combination of surface and canopy models can be easily brought together and moreover analyzed. The analysis of the different model combination within one framework can be seen as the biggest advantage of the developed SenSE package. Currently implemented surface models are: Oh1992, Oh2004, Dubois95 and IEM and Water Cloud. Currently implemented canopy models are: SSRT and Water Cloud.

Within the MULTIPLY project the implemented rt-models will be used for the retrieval of biophysical parameters from the microwave remote sensing data. 

# References
Attema, E.P.W. and Ulaby, F.T., 1978. Vegetation modeled as a water cloud. Radio science, 13(2), pp.357-364.
Baghdadi, N.; Holah, N.; Zribi, M. Calibration of the integral equation model for SAR data in C-band and HH and VV polarizations. Int. J. Remote Sens. 2006, 27, 805–816
Baghdadi, N., Zribi, M. and Bousbih, S., 2017 Calibration of the Water Cloud Model at C-Band for Winter Crop Fields and Grasslands. Remote Sensing, vol 9, 969
Bracaglia, M., et al. (1995). A fully polarimetric multiple scattering model for crops. Remote Sensing of Environment, 54(3), 170-179.
Chen, K.-S., Wu, T.-D., Tsang, L., Li, Q., Shi, J. and Fung, A.K., 2003, Emission of rough surfaces calculated by the integral equation method with comparison to three-dimensional moment method simulations. IEEE Trans. Geosci. Remote Sens. 41, 90–101
Choker, M., Baghdadi, N., Zribi, M. El Hajj ,M., Paloscia,S., Verhoest, N.E. C. , Lievens, H. and Mattia, F. 2017. Evaluation of the Oh, Dubois and IEM Backscatter Models Using a Large Dataset of SAR Data and Experimental Soil Measurements, Water, vol  9, 38
DeRoo, R.D., Du, Y., and Ulaby, F.T. 2001. A semi-empirical backscattering model at L-Band and C-Band for a soybean canopy with soil moisture inversion.  IEEE Transactions on Geoscience and Remote Sensing, 39(4), pp.864-872.
Dubois, P.C., Van Zyl, J. and Engman, T., (1995). Measuring soil moisture with imaging radars. IEEE Transactions on Geoscience and Remote Sensing, 33(4), pp.915-926.
Ferrazzoli, P. and Guerriero, L., (1995). Radar sensitivity to tree geometry and woody volume: A model analysis. IEEE Transactions on Geoscience and Remote Sensing, 33(2), pp.360-371.
Ferrazzoli, P. and Guerriero, L., (1996). Passive microwave remote sensing of forests: A model investigation. IEEE Transactions on Geoscience and Remote Sensing, 34(2), pp.433-443.
Jagdhuber, T.; Hajnsek, I. and Papathanassiou, K. P. , 2015. An iterative Generalized Hybrid Decomposition for Soil Moisture Retrieval Under Vegetation Cover Using Fully Polarimetric SAR IEEE Journal of Selected Topics in Applied Earth Observations and Remote Sensing, 8, 3911-3922
Loew, A. 2004. Coupled modelling of land surface microwave interactions using ENVISAT ASAR data Ludwig-Maximilian University Munich
Oh, Y., Sarabandi, K. and Ulaby, F.T., (1992). An empirical model and an inversion technique for radar scattering from bare soil surfaces. IEEE transactions on Geoscience and Remote Sensing, 30(2), pp.370-381.
Oh, Y., (2004). Quantitative retrieval of soil moisture content and surface roughness from multipolarized radar observations of bare soil surfaces. IEEE Transactions on Geoscience and Remote Sensing, 42(3), pp.596-601.
Paloscia, S., Pettinato S. and E Santi,“Combining L and X band SAR data for estimating biomass and soil moisture of agricultural fields”, European Journal of Remote Sensing 45: 99-109, 2012
Quast, R. and Wagner, W. (2016): "Analytical solution for first-order scattering in bistatic radiative transfer interaction problems of layered media," Appl. Opt. 55, 5379-5386.
Trouve, 1994.Adaptation of the MIMICS Backscattering Model to the Agricultural Context - Wheat and Canola at L and C Bands IEEE Transactions on Geoscience and Remote Sensing, vol. 32, no. 1, pp. 47-61, 
Ulaby, F.T., Sarabandi, K., McDonald, K., Whitt, M., Dobson, M.C. (1990): Michigan microwave canopy scattering model (MIMICS). Int. J. Remote Sensing. Vol. 11 (7). pp. 1223-1253.
Ulaby, F. T. and D. G. Long, 2014. Microwave Radar and Radiometric Remote Sensing, University of Michigan Press 
