# IVIM_Fitting_Methods
Different bi-exponential fit methods for IVIM signal. 

See IVIM_Fit_Algorithm.pptx for overview of the quick simulation shown here. 


AlgorithmN.m:

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


ALGORITHM 5
- segmented fit with f and D held, D* fit in the second run

ALGORITHM 6
- segmented fit with f and D held, D* fit in the second run
- also forced D* upper bound (unlike algo5)

ALGORITHM 7
- segmented fit with f and D held, D* fit in the second run
- but also forced D* upper bound (unlike algo5)
- and not in log space (unlike algo6)


RunIVIM_Algos:
- overall function that gets a stack of dicoms and processes it into a volume of IVIM parameters, and also creates a brain mask.
- It saves both in the same folder as the dicoms
- It does not perform any smoothing or averaging beyond averaging across different diffusion directions
- Multiple the maps by the brain mask if you want them applied.
- To run: 
>> dicom_path = path/to/dicom
>> dicom = name_of_dicom
>> IVIM_maps = RunIVIM_Algos(dicom_path, dicom)

It should create a nx by ny by nz by 6 volume of images with 3D volumes for f, D, Dstar, SSE, rsq, and adj_rsq and it saves it in the same folder as the dicom.


view3DVolume: 
- view the resulting maps with slider on the bottom (slices), on the left (minimum contrast range) and on the right (maximum contrast range).
- input the 3D volume of interest (e.g. squeeze(IVIM_maps(:,:,:,1)) to look at the 3D volume if IVIM f parameter. 


AnomalousDiffusion_BiexponentialFitSimulation.ipynb
Python 3.11 (Anaconda inc) code that generates anisotropic anomalous diffusion signal using simulated ellipsoids and the anomalous diffusion bi-exponential equation.

Simulated_IVIM_Curves_20250202-1921 is a set of 1000 runs. Fit_Simulated_IVIM_CUrves_20250202-1921 is the corresponding fits of the generated signal using the four different algorithms. 
The - 2143 is the same, for both excel sheets, just with 2000 instead for a larger set. 

Tests.m
Was just the tests used to make sure algorithmN.m were outputting the same values

%% note to self, you need to check the startpoint order of parameters with the different matlab versions.
