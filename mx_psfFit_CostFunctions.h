// Author: Simon Christoph Stein
// E-Mail: scstein@phys.uni-goettingen.de
// Date: June 2015

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

#pragma once

using ceres::Jet;

// TODO: Replace finite differencing in MLE cost functions by automatic differentiation or analytic derivative.
//       This will be much faster and more accurate.


// Compute the complementary error function 1-erf(x)
// Templated for automatic differentiation with ceres
template<typename T>
        T erfc(T x)
{
    // Accurate to 1e-7
    // constants
    const T a1 =  T(0.254829592);
    const T a2 =  T(-0.284496736);
    const T a3 =  T(1.421413741);
    const T a4 =  T(-1.453152027);
    const T a5 =  T(1.061405429);
    const T p  =  T(0.3275911);
    
    // Save the sign of x
    int sign =  (x<T(0))?-1:1;
    x = abs(x);
    
    // Abramowitz and Stegun (equations 7.1.25–28), see Wikipedia "error function"
    T t = T(1.0)/(T(1.0) + p*x);
    T y = T(1.0) - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x);
    
    return T(1) - T(sign)*y;
};



//////////////////////////    Cost Functions    //////////////////////////

struct SampledGaussResidual {
    SampledGaussResidual(Array2D& I, std::vector<int>& x, std::vector<int>& y)
    : I_(I), x_(x), y_(y)
    {}
    
    // Parameters x,y,A,BG,sigma
    template <typename T>
            bool operator()(const T* const xpos, const T* const ypos, const T* const A, const T* const BG, const T* const sigma, T* residual) const
    {
        int cnt=0;
        T* xVals = new T[I_.nCols];
        T* yVals = new T[I_.nRows];
        
        // Here we make use of the fact that the gaussian is seperable
        for(unsigned int iCol = 0; iCol<I_.nCols; ++iCol)
            xVals[iCol] = exp( -(xpos[0]-T(x_[iCol]))*(xpos[0]-T(x_[iCol]))/(T(2)*sigma[0]*sigma[0]) );
        for(unsigned int iRow = 0; iRow<I_.nRows; ++iRow)
            yVals[iRow] = exp( -(ypos[0]-T(y_[iRow]))*(ypos[0]-T(y_[iRow]))/(T(2)*sigma[0]*sigma[0]) );
        
        for(unsigned int iCol = 0; iCol<I_.nCols; ++iCol)
            for(unsigned int iRow = 0; iRow<I_.nRows; ++iRow)
            {
                residual[cnt] = I_(iRow,iCol) - A[0]*xVals[iCol]*yVals[iRow] - BG[0];
                ++cnt;
            }
        
        delete[] xVals;
        delete[] yVals;
        
        return true;
    }
    
    // Observations of a sample
    Array2D& I_;
    const std::vector<int>& x_; // x coordinates
    const std::vector<int>& y_; // y coordinates
};


struct IntegratedGaussResidual {
    IntegratedGaussResidual(Array2D& I, std::vector<int>& x, std::vector<int>& y)
    : I_(I), x_(x), y_(y)
    {}
    
    // Parameters x,y,A,BG,sigma
    template <typename T>
            bool operator()(const T* const xpos, const T* const ypos, const T* const A, const T* const BG, const T* const sigma, T* residual) const
    {
        T* xVals = new T[I_.nCols];
        T* yVals = new T[I_.nRows];
        
        int cnt=0;
        
        // Here we make use of the fact that the gaussian is separable
        for(unsigned int iCol = 0; iCol<I_.nCols; ++iCol)
            xVals[iCol] = ( erfc<T>(-(T(x_[iCol])-xpos[0]+T(0.5))/(T(sqrt(2))*sigma[0])) - erfc<T>(-(T(x_[iCol])-xpos[0]-T(0.5))/(T(sqrt(2))*sigma[0])) );
        for(unsigned int iRow = 0; iRow<I_.nRows; ++iRow)
            yVals[iRow] = ( erfc<T>(-(T(y_[iRow])-ypos[0]+T(0.5))/(T(sqrt(2))*sigma[0])) - erfc<T>(-(T(y_[iRow])-ypos[0]-T(0.5))/(T(sqrt(2))*sigma[0])) );
        
        for(unsigned int iCol = 0; iCol<I_.nCols; ++iCol)
            for(unsigned int iRow = 0; iRow<I_.nRows; ++iRow)
            {
                residual[cnt] = I_(iRow,iCol) - A[0]*M_PI*sigma[0]*sigma[0]/T(2) * xVals[iCol]*yVals[iRow] - BG[0];
                ++cnt;
            }
        
        delete[] xVals;
        delete[] yVals;
        
        return true;
    }
    
    // Observations of a sample
    Array2D& I_;
    const std::vector<int>& x_; // x coordinates
    const std::vector<int>& y_; // y coordinates
};



class SampledGauss_MLE_Cost : public ceres::FirstOrderFunction {
    
public:
    SampledGauss_MLE_Cost(Array2D& I, std::vector<int>& x, std::vector<int>& y, std::vector<int>& param_optimMask, double xpos_init, double ypos_init, double A_init, double BG_init, double sigma_init )
    : I_(I), x_(x), y_(y), param_optimMask_(param_optimMask), nr_optim_params_(0), xpos_(xpos_init, 0), ypos_(ypos_init, 1), A_(A_init, 2), BG_(BG_init, 3), sigma_(sigma_init, 4)
    {
        for(int i =0; i<param_optimMask_.size(); ++i)
            nr_optim_params_ += (param_optimMask_[i]>0);
        
        // Set infinitesimals of constant parameters to zero
        if(param_optimMask_[0]>0) xpos_.v[0] = 0;
        if(param_optimMask_[1]>0) ypos_.v[1] = 0;
        if(param_optimMask_[2]>0) A_.v[2] = 0;
        if(param_optimMask_[3]>0) BG_.v[3] = 0;
        if(param_optimMask_[4]>0) sigma_.v[4] = 0;
    }
    
    virtual bool Evaluate(const double* parameters,
            double* cost,
            double* gradient) const {
        
        { // grab values from parameter vector
            int cnt = 0;
            if(param_optimMask_[0]>0){ xpos_.a  = parameters[cnt]; ++cnt;};
            if(param_optimMask_[1]>0){ ypos_.a  = parameters[cnt]; ++cnt;};
            if(param_optimMask_[2]>0){ A_.a     = parameters[cnt]; ++cnt;};
            if(param_optimMask_[3]>0){ BG_.a    = parameters[cnt]; ++cnt;};
            if(param_optimMask_[4]>0){ sigma_.a = parameters[cnt];};
        }

        if(gradient == NULL) // if only function value is requested
        {
            cost[0] =  cost_function<double>(xpos_.a, ypos_.a, A_.a, BG_.a, sigma_.a);
        }
        else // if gradient is requested
        {
            Jet<double, 5> J_result = cost_function(xpos_, ypos_, A_, BG_, sigma_);
            cost[0] = J_result.a;
            
            // Assign derivatives to result
            int cnt = 0;
            if(param_optimMask_[0]>0){ gradient[cnt] = J_result.v[0]; ++cnt;};
            if(param_optimMask_[1]>0){ gradient[cnt] = J_result.v[1]; ++cnt;};
            if(param_optimMask_[2]>0){ gradient[cnt] = J_result.v[2]; ++cnt;};
            if(param_optimMask_[3]>0){ gradient[cnt] = J_result.v[3]; ++cnt;};
            if(param_optimMask_[4]>0){ gradient[cnt] = J_result.v[4]; ++cnt;};
        };
        
        return true;
    }
    
    virtual int NumParameters() const { return nr_optim_params_; }
    
    template <typename T>
            T cost_function( const T& xpos, const T& ypos, const T& A, const T& BG, const T& sigma) const
    {
        T* xVals = new T[I_.nCols];
        T* yVals = new T[I_.nRows];
        
        // Here we make use of the fact that the gaussian is seperable
        for(unsigned int iCol = 0; iCol<I_.nCols; ++iCol)
            xVals[iCol] = exp( -(xpos-T(x_[iCol]))*(xpos-T(x_[iCol]))/(T(2)*sigma*sigma) );
        for(unsigned int iRow = 0; iRow<I_.nRows; ++iRow)
            yVals[iRow] = exp( -(ypos-T(y_[iRow]))*(ypos-T(y_[iRow]))/(T(2)*sigma*sigma) );
        
        T cost(0.0);
        T model_val(0.0);
        for(unsigned int iCol = 0; iCol<I_.nCols; ++iCol)
            for(unsigned int iRow = 0; iRow<I_.nRows; ++iRow)
            {
                model_val = A*xVals[iCol]*yVals[iRow] + BG;
                cost = I_(iRow,iCol) - model_val;
            }
        
        delete[] xVals;
        delete[] yVals;
        
        return cost;
    }
    
private:
    // Observations of a sample
    Array2D& I_;
    const std::vector<int>& x_; // x coordinates
    const std::vector<int>& y_; // y coordinates
    
    const std::vector<int>& param_optimMask_; // mask which parameters are optimized
    int nr_optim_params_;
    
    // To keep different parameters constant, we copy all initial values.
    // We use automatic differentiation here, so all variables are type Jet<double, nr_parameters>
    // Note: Need mutable because we need to change these inside the const. Evaluate function.
    mutable Jet<double,5> xpos_, ypos_, A_, BG_, sigma_;
};


class IntegratedGauss_MLE_Cost : public ceres::FirstOrderFunction {
    
public:
    IntegratedGauss_MLE_Cost(Array2D& I, std::vector<int>& x, std::vector<int>& y, std::vector<int>& param_optimMask, double xpos_init, double ypos_init, double A_init, double BG_init, double sigma_init )
    : I_(I), x_(x), y_(y), param_optimMask_(param_optimMask), nr_optim_params_(0), xpos_(xpos_init, 0), ypos_(ypos_init, 1), A_(A_init, 2), BG_(BG_init, 3), sigma_(sigma_init, 4)
    {
        for(int i =0; i<param_optimMask_.size(); ++i)
            nr_optim_params_ += (param_optimMask_[i]>0);
        
        // Set infinitesimals of constant parameters to zero
        if(param_optimMask_[0]>0) xpos_.v[0] = 0;
        if(param_optimMask_[1]>0) ypos_.v[1] = 0;
        if(param_optimMask_[2]>0) A_.v[2] = 0;
        if(param_optimMask_[3]>0) BG_.v[3] = 0;
        if(param_optimMask_[4]>0) sigma_.v[4] = 0;
    }
    
    virtual bool Evaluate(const double* parameters,
            double* cost,
            double* gradient) const {
        
        { // grab values from parameter vector
            int cnt = 0;
            if(param_optimMask_[0]>0){ xpos_.a  = parameters[cnt]; ++cnt;};
            if(param_optimMask_[1]>0){ ypos_.a  = parameters[cnt]; ++cnt;};
            if(param_optimMask_[2]>0){ A_.a     = parameters[cnt]; ++cnt;};
            if(param_optimMask_[3]>0){ BG_.a    = parameters[cnt]; ++cnt;};
            if(param_optimMask_[4]>0){ sigma_.a = parameters[cnt];};
        }

        if(gradient == NULL) // if only function value is requested
        {
            cost[0] =  cost_function<double>(xpos_.a, ypos_.a, A_.a, BG_.a, sigma_.a);
        }
        else // if gradient is requested
        {
            Jet<double, 5> J_result = cost_function(xpos_, ypos_, A_, BG_, sigma_);
            cost[0] = J_result.a;
            
            // Assign derivatives to result
            int cnt = 0;
            if(param_optimMask_[0]>0){ gradient[cnt] = J_result.v[0]; ++cnt;};
            if(param_optimMask_[1]>0){ gradient[cnt] = J_result.v[1]; ++cnt;};
            if(param_optimMask_[2]>0){ gradient[cnt] = J_result.v[2]; ++cnt;};
            if(param_optimMask_[3]>0){ gradient[cnt] = J_result.v[3]; ++cnt;};
            if(param_optimMask_[4]>0){ gradient[cnt] = J_result.v[4]; ++cnt;};
        };
        return true;
    }
    
    virtual int NumParameters() const { return nr_optim_params_; }
    
    template <typename T>
            T cost_function( const T& xpos, const T& ypos, const T& A, const T& BG, const T& sigma) const
    {
        T* xVals = new T[I_.nCols];
        T* yVals = new T[I_.nRows];
        
        // Here we make use of the fact that the gaussian is separable
        for(unsigned int iCol = 0; iCol<I_.nCols; ++iCol)
            xVals[iCol] = ( erfc<T>(-(T(x_[iCol])-xpos+T(0.5))/(T(sqrt(2))*sigma)) - erfc<T>(-(T(x_[iCol])-xpos-T(0.5))/(T(sqrt(2))*sigma)) );
        for(unsigned int iRow = 0; iRow<I_.nRows; ++iRow)
            yVals[iRow] = ( erfc<T>(-(T(y_[iRow])-ypos+T(0.5))/(T(sqrt(2))*sigma)) - erfc<T>(-(T(y_[iRow])-ypos-T(0.5))/(T(sqrt(2))*sigma)) );
        
        T cost(0.0);
        T model_val (0.0);
        for(unsigned int iCol = 0; iCol<I_.nCols; ++iCol)
            for(unsigned int iRow = 0; iRow<I_.nRows; ++iRow)
            {
                model_val = A*M_PI*sigma*sigma/double(2) * xVals[iCol]*yVals[iRow] + BG   + 1e-14;
                cost -= I_(iRow,iCol) * log(model_val) - model_val;
                // Note -= because we evaluate the negative log likelihood
            }
        
        delete[] xVals;
        delete[] yVals;
        
        return cost;
    }
    
private:
    // Observations of a sample
    Array2D& I_;
    const std::vector<int>& x_; // x coordinates
    const std::vector<int>& y_; // y coordinates
    
    const std::vector<int>& param_optimMask_; // mask which parameters are optimized
    int nr_optim_params_;
    
    // To keep different parameters constant, we copy all initial values.
    // We use automatic differentiation here, so all variables are type Jet<double, nr_parameters>
    // Note: Need mutable because we need to change these inside the const. Evaluate function.
    mutable Jet<double,5> xpos_, ypos_, A_, BG_, sigma_;  
};


// An example for setting up an optimization problem with per-pixel residuals
// //  Optimization using the per pixel functions (slow but easy implementation)
//   for(unsigned int iCol = 0; iCol<img.nCols; ++iCol)
//       for(unsigned int iRow = 0; iRow<img.nRows; ++iRow)
//       {
// //           Set up the residual per pixel. This uses
// //           auto-differentiation to obtain the derivative (jacobian).
//           CostFunction* cost_function =
//               new AutoDiffCostFunction<SampledGaussResidual_PerPix, 1, 1,1,1,1,1>(  new SampledGaussResidual_PerPix(iCol,iRow, img(iRow,iCol))  );
//           problem.AddResidualBlock(cost_function, NULL, &xpos, &ypos, &A, &BG, &sigma);
//       }

// A templated cost functor that implements the residual at every pixel.
struct SampledGaussResidual_PerPix {
    SampledGaussResidual_PerPix(double x, double y, double I)
    : x_(x), y_(y), I_(I) {}
    
    // Parameters x,y,A,BG,sigma
    template <typename T>
            bool operator()(const T* const xpos, const T* const ypos, const T* const A, const T* const BG, const T* const sigma, T* residual) const
    {
        residual[0] = T(I_) - A[0] *exp( -((xpos[0]-T(x_))*(xpos[0]-T(x_))+ (ypos[0]-T(y_))*(ypos[0]-T(y_)) )/(T(2)*sigma[0]*sigma[0]) ) - BG[0];
        return true;
    }
    
    // Observations of a sample
    const double x_;
    const double y_;
    const double I_;
};


// A templated cost functor that implements the residual at every pixel.
struct IntegratedGaussResidual_PerPix {
    IntegratedGaussResidual_PerPix(double x, double y, double I)
    : x_(x), y_(y), I_(I) {}
    
    // Parameters x,y,A,BG,sigma
    template <typename T>
            bool operator()(const T* const xpos, const T* const ypos, const T* const A, const T* const BG, const T* const sigma, T* residual) const
    {
        residual[0] = T(I_) - A[0]*M_PI*sigma[0]*sigma[0]/T(2) * ( erfc<T>(-(T(x_)-xpos[0]+T(0.5))/(T(sqrt(2))*sigma[0])) - erfc<T>(-(T(x_)-xpos[0]-T(0.5))/(T(sqrt(2))*sigma[0]))    ) \
                * ( erfc<T>(-(T(y_)-ypos[0]+T(0.5))/(T(sqrt(2))*sigma[0])) - erfc<T>(-(T(y_)-ypos[0]-T(0.5))/(T(sqrt(2))*sigma[0]))   ) - BG[0] ;
        return true;
    }
    
    // Observations of a sample
    const double x_;
    const double y_;
    const double I_;
};