Lag_mod is a C++ library implementing a hierarchical statistical model for 3D Laguerre tessellation

The model is introduced in paper "Fitting three-dimensional Laguerre tessellations by hierarchical marked point process models"
available at arXiv:

dependencies: 
  Voro++ (http://math.lbl.gov/voro++/ and https://github.com/chr1shr/voro) is used for computations of Laguerre tessellation
  Eigen
  
library files: Header.h data_creator.cpp dens_ener.cpp feasibility.cpp helping_fcs.cpp LAG_ORI_recompute.cpp LAG_recompute.cpp numeric.cpp orientations.cpp pipp.cpp prob_distributions.cpp radii.cpp statistics.cpp 

the capabilities of the library are demonstrated on several examples:

Source_Lag_tess.cpp - shows how to compute a Laguerre tessellation for a given list of generators using Voro++ 

Source_PP_sim.cpp - simulation of a 3D point pattern

pipp_estim.cpp

rad_sim.cpp

rad_estim.cpp






** the code still undergoes extensive changes **
