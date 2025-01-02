
% REQUIREMENTS: install chebfun: chebfun.org
% codes uses diffmat.m to generate differentiation matrices and chebpts.m
% to generate Chebyshev grids

To change model parameters, edit reverbTankParams.m and save.
 - Default set up is for a small spring from the Olson X-82 amplifier for faster run time.
 - Parameters for both springs modelled in the paper are available in reverbTankParams.m.
 - Remember to change the magnetic bead setup, not just the helical spring parameters.

To generate results, run modalAlg.m.
 - This can take a long time -- the matrix size for the Olson spring takes around 5 minutes,
   but the number of points specified in the paper for the (longer) springs in the Accutronics tank
   can take upwards of 15 minutes.
 - modalAlg.m calls tankModelEig.m to extract eigenvalues for the given setup.