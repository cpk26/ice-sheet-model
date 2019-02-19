# DEpth integrated Ice Sheet Model (DEISM)
This repository contains the code for the ice sheet model I wrote during my PhD thesis. The model is written in Matlab, implements the ice flow equations derived in [1], and is numerically based on [2]. The model uses automatic differentiation to perform inversions for basal drag. The model has been tested on the ISMIP-HOM experiments to verify its implementation. For full details of the model, see [3] or [4].

# Implementation
The model solves depth integrated ice flow equations presented in [1], often termed the 'hybrid equations'. The name arises from the fact that it combines the Shallow Ice Approximation (SIA) and Shallow Shelf Approximation (SSA) of ice sheet equations. This formulation allows you to solve for ice sheet velocities on a 2D-grid, and then reconstruct the vertical velocity profile after the fact. The equations are solved using finite difference methods on a Arakawa-C grid. 

The model can be used to invert for basal drag based on observations of surface velocities and topography. The inversion procedure is compatable with all three sliding laws implemented (Linear, Budd, and Schoof). A novel procedure for incorporating basal hydrology into the inversions for basal drag using the Budd and Schoof laws is described in [4]. In [4] and [5], the ice flow model was integrated with a subglacial hydrology model described in [6]. Please contact the author of [6] if you'd like a copy of the hydrology model. I am happy to share the code which integrates the two.

This model is written in Matlab 2016b, and uses the [ADiGator](https://sourceforge.net/projects/adigator/) automatic differentiation tool.

# Getting Started
For run files and help getting started contact me via Github, or alternatively, message [Dr. Neil Arnold](https://www.geog.cam.ac.uk/people/arnold/) at the Scott Polar Research institute (University of Cambridge). Due to the pressure of finishing my thesis, I never had the time to write a guide to using the model. However, the code is well documented, and I'd be happy to help you get started with it. 

# References

[1] Goldberg, Daniel N. (2011). A variationally derived, depth-integrated approximation to a higher-order glaciological flow model. Journal of Glaciology 57.201 : 157-170.

[2] Arthern, Robert J., Hindmarsh, Richard C.A., Williams, C. Rosie. (2015) Flow speed within the Antarctic ice sheet and its controls inferred from satellite observations. Journal of Geophysical Research, 120. 1171-1188. 10.1002/2014JF003239

[3] Koziol, C. P. (2018). Modelling the impact of surface melt on the hydrology and dynamics of the Greenland Ice Sheet (Doctoral thesis). University of Cambridge. https://doi.org/10.17863/CAM.20372 

[4]  Koziol, C., and Arnold, N. (2017). Incorporating modelled subglacial hydrology into inversions for basal drag. Cryosphere, 11 (6), 2783-2797. https://doi.org/10.5194/tc-11-2783-2017 

[5] Koziol, Conrad P., and Neil Stuart Arnold (2018). Modelling seasonal meltwater forcing of the velocity of the Greenland Ice Sheet. Cryosphere, 12 (3), 971-991. https://doi.org/10.5194/tc-12-971-2018

[6] Hewitt, I. J. (2013). Seasonal Changes in Ice Sheet Motion due to Meltwater Lubrication. Earth and Planetary Science Research Letters 371-372, 16-25 doi:10.1016/j.epsl.2013.04.022

If you find this work useful, please consider citing reference [4].

