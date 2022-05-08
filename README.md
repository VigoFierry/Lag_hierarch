Lag_mod is a C++ library implementing various statistical models and related computations of 3D Laguerre tessellations

dependencies: 
  Voro++ (http://math.lbl.gov/voro++/ and https://github.com/chr1shr/voro) is used for computations of Laguerre tessellation;
  Eigen
  
library files: Header.h data_creator.cpp dens_ener.cpp feasibility.cpp helping_fcs.cpp LAG_ORI_recompute.cpp LAG_recompute.cpp numeric.cpp orientations.cpp pipp.cpp prob_distributions.cpp radii.cpp statistics.cpp subcells.cpp



The library consists of several parts (based on functionality; all of these are accompanied by examples):

1) a construction of hierarchical statistical model for 3D Laguerre tessellation as introduced in paper "Fitting three-dimensional Laguerre tessellations by hierarchical marked point process models" available at arXiv: https://arxiv.org/abs/2110.07485

related examples:

Source_Lag_tess.cpp - shows how to compute a Laguerre tessellation for a given list of generators using Voro++ 

Source_PP_sim.cpp - simulation of a 3D point pattern

Source_PP_estim.cpp - estimation of parameters of a given point process model

Source_rad_sim.cpp - simulation of radii conditioned on a point pattern

Source_rad_estim.cpp - estimation of parameters of a radii model






** the code still undergoes extensive changes **
