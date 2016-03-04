
 Matlab files:
psfFit.m       - Fit a gaussian point spread function model to the given image.
psfFit_Image.m - Fit multiple spot candidates in a given image with a gaussian point spread function model.
EXAMPLE.m      - Simple demonstration using psfFit_Image.


The fitting code uses the ceres-solver library for optimization which is developed by Google.
  URL: http://ceres-solver.org


- Pre-compiled files for Windows 7 64-bit -

Compiled mex files for 64-bit Windows built using Visual Studio Professional 2012 can be downloaded from:
     https://drive.google.com/file/d/0BzFq6JJEnf3fMFdUbGNNaDdyUDg/view?usp=sharing
For the precompiled mex files to work you may need the "Visual C++ Redistributable for Visual Studio 2012 Update 4" which you can download from www.microsoft.com. 
Make sure all files are in the directory of the Matlab functions before calling them.


- Compiling source code by yourself -

(x) Go to the FastPsfFitting directory
   For compiling ceres we use Tal Ben-Nun's github repository for windows (we use it for linux as well)
   Execute "git submodule init" and "git submodule update" to download it automatically
(x) Download the eigen library from
     URL: http://eigen.tuxfamily.org/ 
   and put it in "ceres-windows/eigen"

Windows:
To build ceres with Visual Studio simply open the project files in ceres-windows and 
start the compilation of project "ceres" in "Release" "x64" mode.

Linux (tested with Kubuntu):
With Linux ceres can be compiled with CMake as follows:
- Execute linux_install_forCeresBuild.sh to install neccessary libraries
- Switch to the "ceres-windows/ceres-solver/" directory
- Execute "mkdir build && cd build"
- Execute "cmake-gui ..", press "Configure" once, set "BUILD_SHARED_LIBS" to true, press "Generate"
- Close cmake-gui and execute "make -j4"
To build ceres with Visual Studio simply open the project files in ceres-windows and start the compilation of project "ceres" in Release mode.

Given a working MEX setup with Visual Studio, the FastPsfFitting functions can be compiled from MATLAB using the provided build_xxx.m-scripts.

