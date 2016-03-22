Repository of the MATLAB file exchange contribution "Fast Gaussian Point Spread Function Fitting (MEX)"
http://www.mathworks.com/matlabcentral/fileexchange/52417-fast-gaussian-point-spread-function-fitting--mex-

The quick and accurate localization of many emitters showing up as 2D Gaussian like shapes in an image is an important tool in fluorescence microscopy. Mostly this is used in the context of "localization microscopy" which is able to yield super-resolved (resolution below the diffraction limit) images of biological samples. Popular tools for this task include RapidStorm and QuickPalm (ImageJ).
This software allows quick and accurate point spread function fitting using a MEX file interface for use directly in MATLAB programs. This is meant to facilitate the development of new and better customizable methods, as Matlab based fitting is usually much too slow for the amount of data that needs to be processed. The fitting code utilizes the ceres-solver library for optimization currently developed by Google (2016).
After a list of candidate positions is provided by the user the fitter returns the parameters [xpos; ypos; Amplitude; local background; standard deviation_x, standard_deviation_y, angle [degree], errorflag] for each candidate. Initial guesses can either be user-supplied or are estimated by the algorithm. Arbitrary parameters can be set fixed for the optimization. By default an isotropic Gaussian is fitted (only xpos, ypos, Amplitude, background; standard deviation).

The Gaussian PSF model can either be taken as point wise sampled (at pixel centers) or pixel integrated (usually the better fit for data recorded with a camera). Optionally a Poissonian noise based Maximum Likelihood refinement is performed after the initial least squared fit, improving accuracy of the fit at low light levels. Note that the input must be in photons rather than camera counts for this to be useful.

For a quickstart precompiled MEX/DLL Files for Windows 7 64-bit can be downloaded from:
  https://drive.google.com/file/d/0BzFq6JJEnf3fTVlDdkhjLWlVY0E/view?usp=sharing
You also need need the "Visual C++ Redistributable for Visual Studio 2012 Update 4" which you can download from www.microsoft.com.

Precompiled MEX/SO Files for Linux can be downloaded from:
  https://drive.google.com/file/d/0BzFq6JJEnf3fMUg0LWFtY0FpVGs/view?usp=sharing
Read the attached README2_Linux.txt to get it running.

If you want to compile the project yourself, follow the steps illustrated in the README.txt.