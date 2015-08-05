
 Matlab files:
psfFit.m       - Fit a gaussian point spread function model to the given image.
psfFit_Image.m - Fit multiple spot candidates in a given image with a gaussian point spread function model.
EXAMPLE.m      - Simple demonstration using psfFit_Image.

Compiled mex files are provided for 64-bit Windows built using Visual Studio Professional 2012.
For the precompiled mex files to work you may need the "Visual C++ Redistributable for Visual Studio 2012 Update 4" which you can download from www.microsoft.com. 
Also make sure the "ceres.dll" library lies in the same directory as the fitting functions for them to work.

The fitting code uses the ceres-solver library for optimization which is developed by Google.
  URL: http://ceres-solver.org


- Compiling source code by yourself -

(x) For windows download the VS project of ceres from Tal Ben-Nun's github repository, see
  URL: https://github.com/tbennun/ceres-windows
Put the "ceres-windows" folder next to the "FastPsfFitting" folder. 
(x) Download ceres (v 1.10.0 at time of compilation of this project) from 
  URL: http://ceres-solver.org" 
and put it in "ceres-windows/ceres-solver.
(x) Download the eigen library from
  URL: http://eigen.tuxfamily.org/ 
and put it in "ceres-windows/eigen"
(x) Download the google glog library from 
  URL: http://code.google.com/p/google-glog/  OR get the git repository -> URL: https://github.com/google/glog
and put it in "ceres-windows/glog

To build ceres with Visual Studio simply open the project files in ceres-windows and start the compilation of project "ceres" in Release mode.

Given a working MEX setup with Visual Studio, the FastPsfFitting functions can be compiled from MATLAB using the provided build scripts.
If you don't want to carry the shared library with you, you can compile with static linking (see the build scripts).
Make sure you compiled the project "ceres_static" beforehand so that the static library is available.

With Linux ceres can be compiled with CMake using the CMakeLists.txt provided by the ceres creators.