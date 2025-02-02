# IVIM_Fitting_Methods
Different bi-exponential fit methods for IVIM signal

ALGORITHM 1 
- Mira original code 
- subtracts f0 and D mono-exponential from the total signal, and then fits f and D* to the remainder, at only b<200.

ALGORITHM 2 
- Salman - old code
- fits D first, then holds D as a constant to fit to f and D* (in log space)

ALGORITHM 3 
- Salman-Mira hybrid code 
- get f0 and D from the high b-values, then fit for f and D* with given f0 and D, at all bvalues.  

ALGORITHM 4 
- Bayesian fit
- minimum mean square error estimator as the mean of the posterior distribution



%% note to self, you need to check the startpoint order of parameters with the different matlab versions.