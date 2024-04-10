*** Data_and_Model_of_Scale_Down_Simulator_for_Industrial_Syngas_Fermentation ***
Authors: L. Puiman, E. Almeida Benalcázar, C. Picioreanu, H.J. Noorman, C. Haringa
Downscaling industrial-scale syngas fermentation to simulate frequent and irregular gas shocks. University of Technology
Corresponding author: C. Haringa
Contact Information:
C.Haringa@tudelft.nl
Delft University of Technology - Faculty of Applied Sciences
Department of Biotechnology, 
Van der Maasweg 9, 
2629 HZ, Delft, 
The Netherlands

***General Introduction***
The data and MATLAB code are provided for the development of the scale-down simulator of industrial-scale syngas fermentation. 
This dataset is part of the publication: Downscaling industrial-scale syngas fermentation to simulate frequent and irregular gas shocks.
The data and script were developed in the Department of Biotechnology at Delft University of Technology, in September 2022.

***Methodological information***
The probability distributions of the residence time (seconds) and concentrations (mol/m3) during peaks and valleys in large-scale reactor are provided from the CFD simulations.
There are separate files for CO and H2, at the three biomass concentrations, they can be read using MATLAB (for example with version R2020b).

For example "prob-peaks-val-clCO-5-gl.mat" contains the probability distributions of the dissolved CO concentrations (mol/m3) and the residence times in the peaks and in the valleys from the simulation with 5 g/L. 
The dataset contains bins and actual probability values. 
The bins contain the ranges of concentration and the corresponding probability of that concentration is provided. 
The same is provided for the residence time (bins in seconds).

Every dataset contains the following vectors (with the compound specified: clCO or clH2):
-   peak_conc_bins_clCO,                [binned values of CO concentrations during the peak (mM)]
-   peak_time_bins_clCO, 
-   prob_peak_conc_clCO,                [probability-values of CO concentrations during the peak (-)]
-   prob_peak_time_clCO, 
-   prob_valley_conc_clCO, 
-   prob_valley_time_clCO, 
-   valley_conc_bins_clCO, 
-   valley_time_bins_clCO

This data is used as input for the file "Scale_down_simulator.m", which is used as conceptual model to scale-down industrial-scale syngas fermentation for a wide range of conditions,but may be used for other purposes.
The operating conditions, kinetics, geometry, etc, of the scale-down simulator might be adjusted to your system.

***Terms of use***
The data and script are free to use under the CC BY 4.0 license. 