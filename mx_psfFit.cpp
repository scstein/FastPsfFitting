
// Author: Simon Christoph Stein
// E-Mail: scstein@phys.uni-goettingen.de
// Date: June 2015

// Optimization is done using the ceres-solver library.
// http://ceres-solver.org, New BSD license.
// Copyright 2015 Google Inc. All rights reserved.

// % Copyright (c) 2015, Simon Christoph Stein
// % All rights reserved.
// % 
// % Redistribution and use in source and binary forms, with or without
// % modification, are permitted provided that the following conditions are met:
// % 
// % 1. Redistributions of source code must retain the above copyright notice, this
// %    list of conditions and the following disclaimer.
// % 2. Redistributions in binary form must reproduce the above copyright notice,
// %    this list of conditions and the following disclaimer in the documentation
// %    and/or other materials provided with the distribution.
// % 
// % THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// % ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// % WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// % DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
// % ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// % (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// % LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// % ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// % (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// % SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// % 
// % The views and conclusions contained in the software and documentation are those
// % of the authors and should not be interpreted as representing official policies,
// % either expressed or implied, of the FreeBSD Project.


// Use OpenMP if available to parallelize computation
#ifdef _OPENMP
    #include <omp.h>
#endif

// mex
#include "mex.h"

// stdlib
#include <cmath>
#include <math.h>

// our headers
#include "mexUtil.h"

// -- Type definitions (must go before further inclusions!) --
typedef Array1D_t<double> Array1D;
typedef Array2D_t<double> Array2D;
typedef Array3D_t<double> Array3D;

#include "mx_psfFit_OptimizerFunctions.h"

using namespace std;

// --   psfFit  --
// % Short usage: params = psfFit( img );
// % Full usage: [ params, exitflag ] = psfFit( img, param_init, param_optimizeMask, useIntegratedGauss, useMLErefine)
// %  Fit a gaussian point spread function model to the given image.
// % 
// % Coordinate convention: Integer coordinates are centered in pixels. I.e.
// % the position xpos=3, ypos=5 is exactly the center of pixel img(5,3). Thus
// % pixel center coordinates range from 1 to size(img,2) for x and 1 to
// % size(img,1) for y coordinates.
// %
// % Use empty matrix [] for parameters you don't want to specify.
// %
// % Input:
// %   img - NxM image to fit to. 
// %   param_init - Initial values [xpos,ypos,A,BG,sigma]. If negative values
// %                are given, the fitter estimates a value for that
// %                parameter. | default: -1*ones(5,1) -> 'estimate all'
// %   param_optimizeMask - Must be true(1)/false(0) for every parameter [xpos,ypos,A,BG,sigma]. 
// %                Parameters with value 'false' are not fitted. | default: ones(5,1) -> 'optimize all'
// %   useIntegratedGauss - Wether to use pixel integrated gaussian or not | default: false
// %   useMLErefine - Use Poissonian noise based maximum likelihood estimation after
// %                  least squares fit. Make sure to input image intensities in photons 
// %                  for this to make sense. | default: false
// %
// % Output
// %   params - Final parameters [xpos,ypos,A,BG,sigma].
// %   exitflag - Return state of optimizer. Positive = 'good'. Negative = 'bad'.
// %         1 - CONVERGENCE
// %        -1 - NO_CONVERGENCE
// %        -2 - FAILURE
// %
// % Author: Simon Christoph Stein
// % Date: June 2015
// % E-Mail: scstein@phys.uni-goettingen.de
void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    // Replace the std stream with the 'matlab' stream to see cout outputs in the MATLAB console
    mexStream mout;
    std::streambuf *outbuf = std::cout.rdbuf(&mout);
    
    #ifdef DEBUG
        #ifdef _OPENMP
            printf("OpenMP number of threads: %i\n", omp_get_num_threads()); // Test if openMP is included correctly
        #endif
    #endif
    
        /// -- Input parsing -- ///      
    
    /* Check for proper number of arguments. */
    if(nrhs<1)
       mexErrMsgTxt("psfFit is hungry! Please feed me at least an image!");
    
     if(mxGetClassID(prhs[0]) != mxDOUBLE_CLASS)
         mexErrMsgTxt("psfFit: input image datatype must be double!");
            
    /* Map access to input data */
    Array2D img( prhs[0] );
    std::vector<double> param_init(5, -1); // -1 default: estimate value of each parameter
    std::vector<int> param_optimMask(5, 1); // 1 default: optimize every parameter
    bool usePixelIntegratedGauss = false;
    bool useMLErefine = false;
    
    if(nrhs>1)
    {
       Array1D p_init( prhs[1] );
       if(p_init.nElements > 0) 
       {
           if(p_init.nElements != 5)
               mexErrMsgTxt("Specify initial condition for every parameter. [xpos,ypos,A,BG,sigma].\n Use negative values for parameters you don't want to specify!");
           else
           {
            for(int iParam=0; iParam<5; ++iParam)
                param_init[iParam] = p_init[iParam];
           }
       }
    }
    
    if(nrhs>2)
    {
       Array1D p_optimMask( prhs[2] );
       if(p_optimMask.nElements > 0) 
       {
         if(p_optimMask.nElements != 5)
             mexErrMsgTxt("Specify for every parameter if it should be optimized. [xpos,ypos,A,BG,sigma].\n Use 0 to not optimize a parameter.");
         else
         {
            for(int iParam=0; iParam<5; ++iParam)
                param_optimMask[iParam] = p_optimMask[iParam];
         }        
       }       
    }    
    
    if(nrhs>3)
    {
       Array1D usePixelIntegration( prhs[3] );
       if(usePixelIntegration.nElements > 0) 
           usePixelIntegratedGauss = usePixelIntegration[0];
    }
    
    if(nrhs>4)
    {
       Array1D useMLE( prhs[4] );
       if(useMLE.nElements > 0) 
           useMLErefine = useMLE[0];
    }
    
    
        /// -- Here goes the fitting! -- ///
    // Built coordinate system to use
    std::vector<int> xCoords(img.nCols);
    std::vector<int> yCoords(img.nRows);    
    for(unsigned int iCol = 0; iCol<img.nCols; ++iCol)
          xCoords[iCol] = iCol+1; // We use Matlab coordinates
    for(unsigned int iRow = 0; iRow<img.nRows; ++iRow)
          yCoords[iRow] = iRow+1; // We use Matlab coordinates
    
    // Fit image, collect results
    std::vector<double> results = fitPSF(img, xCoords, yCoords, param_init, param_optimMask, usePixelIntegratedGauss, useMLErefine);
    
        ///  --  Output to MATLAB -- //
    const int nDims = 1; // Number of dimensions for output
    mwSize dim_out0[nDims] = { 5 };
    plhs[0] = mxCreateNumericArray( nDims, dim_out0 , mxDOUBLE_CLASS, mxREAL);
    
    // fitted parameter values
    Array1D fin_params ( plhs[0] ); // Map access to output
    fin_params[0] = results[0]; // x
    fin_params[1] = results[1]; // y
    fin_params[2] = results[2]; // A
    fin_params[3] = results[3]; // BG
    fin_params[4] = results[4]; // sigma
    
    // exit state of optimizer
    mwSize dim_out1[1] = { 1 };
    plhs[1] = mxCreateNumericArray( 1, dim_out1 , mxDOUBLE_CLASS, mxREAL);
    double& exitflag  = *mxGetPr(plhs[1]);
    exitflag = results[5]; // exitflag
    
    // Restore the std stream buffer
    std::cout.rdbuf(outbuf);
}


