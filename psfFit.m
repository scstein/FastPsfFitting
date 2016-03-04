function [ params, exitflag ] = psfFit( img, varargin )
% Short usage: params = psfFit( img );
% Full usage: [ params, exitflag ] = psfFit( img, param_init, param_optimizeMask, useIntegratedGauss, useMLErefine)
%  Fit a gaussian point spread function model to the given image.
% 
% Coordinate convention: Integer coordinates are centered in pixels. I.e.
% the position xpos=3, ypos=5 is exactly the center of pixel img(5,3). Thus
% pixel center coordinates range from 1 to size(img,2) for x and 1 to
% size(img,1) for y coordinates.
%
% Use empty matrix [] for parameters you don't want to specify.
%
% Input:
%   img - NxM image to fit to. (internally converted to double)
%   param_init - Initial values [xpos,ypos,A,BG,sigma]. If negative values
%                are given, the fitter estimates a value for that
%                parameter. | default: -1*ones(5,1) -> 'estimate all'
%   param_optimizeMask - Must be true(1)/false(0) for every parameter [xpos,ypos,A,BG,sigma_x,sigma_y,angle].
%                Parameters with value 'false' are not fitted. | default: [1,1,1,1,1,0,0] -> 'optimize x,y,A,BG,sigma' (2D case)
%   useIntegratedGauss - Wether to use pixel integrated gaussian or not | default: false
%   useMLErefine - Use Poissonian noise based maximum likelihood estimation after
%                  least squares fit. Make sure to input image intensities in photons 
%                  for this to make sense. | default: false
%
% Output
%   params - Final parameters [xpos,ypos,A,BG,sigma_x,sigma_y,angle].
%   exitflag - Return state of optimizer. Positive = 'good'. Negative = 'bad'.
%         1 - CONVERGENCE
%        -1 - NO_CONVERGENCE
%        -2 - FAILURE
%
% Optimization is done using the ceres-solver library.
% http://ceres-solver.org, New BSD license.
% Copyright 2015 Google Inc. All rights reserved.
%
% Author: Simon Christoph Stein, extended by Jan Thiart
% Date: June 2015
% E-Mail: scstein@phys.uni-goettingen.de
% 


% -- Licensing --

% Copyright (c) 2015, Simon Christoph Stein
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
% 
% 1. Redistributions of source code must retain the above copyright notice, this
%    list of conditions and the following disclaimer.
% 2. Redistributions in binary form must reproduce the above copyright notice,
%    this list of conditions and the following disclaimer in the documentation
%    and/or other materials provided with the distribution.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
% ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
% WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
% ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
% (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
% LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
% ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
% (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
% 
% The views and conclusions contained in the software and documentation are those
% of the authors and should not be interpreted as representing official policies,
% either expressed or implied, of the FreeBSD Project.


% License ceres-solver:
%
% Copyright 2015 Google Inc. All rights reserved.
%
% Redistribution and use in source and binary forms, with or without modification, 
% are permitted provided that the following conditions are met:
% 
%     1. Redistributions of source code must retain the above copyright notice, this list 
% of conditions and the following disclaimer.
%     2. Redistributions in binary form must reproduce the above copyright notice, this list
% of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
%     3. Neither the name of Google Inc., nor the names of its contributors may be used to
% endorse or promote products derived from this software without specific prior written permission.
%
% This software is provided by the copyright holders and contributors “AS IS” and any express or 
% implied warranties, including, but not limited to, the implied warranties of merchantability and 
% fitness for a particular purpose are disclaimed. In no event shall Google Inc. be liable for any 
% direct, indirect, incidental, special, exemplary, or consequential damages (including, but not limited
% to, procurement of substitute goods or services; loss of use, data, or profits; or business interruption)
% however caused and on any theory of liability, whether in contract, strict liability, or tort (including
% negligence or otherwise) arising in any way out of the use of this software, even if advised of the
% possibility of such damage.

% Make sure logicals are passed as correct datatype
if numel(varargin) >= 2 && ~isempty(varargin{2});  varargin{2} = logical(varargin{2});  end
if numel(varargin) >= 3 && ~isempty(varargin{3});  varargin{3} = logical(varargin{3});  end
if numel(varargin) >= 4 && ~isempty(varargin{4});  varargin{4} = logical(varargin{4});  end

% Convert img to double if neccessary
% sigma_x = varargin{1}(5);
% sigma_y = varargin{1}(6);
% angle = varargin{1}(7);
% q1 = cos(angle).^2./(2*sigma_x.^2)+sin(angle).^2./(2*sigma_y.^2);
% q2 = sin(angle).^2./(2*sigma_x.^2)+cos(angle).^2./(2*sigma_y.^2);
% q3 = -sin(2*angle)./(4*sigma_x.^2)+sin(2*angle)./(4*sigma_y.^2);
% varargin{1}(5:7) = [q1;q2;q3];

% Set default optimization x,y,A,BG,sigma
if( numel(varargin)<2 || isempty(varargin{2}) )
    varargin{2} = logical([1,1,1,1,1,0,0]);
end

[ params, exitflag ] = mx_psfFit( double(img), varargin{:});

% Compute sigma_x, sigma_y, angle from the output values
param_optimizeMask = varargin{2};

fitSigma_y = param_optimizeMask(6);
fitAngle = param_optimizeMask(7);

if fitAngle
    q1 = params(5); q2 = params(6); q3 = params(7);
    angle = 0.5*atan(2*q3./(q2-q1)); %angle towards largest eigenvector axis in radian
    sigma_x = 1./sqrt(q1+q2-2*q3./sin(2*angle));
    sigma_y = 1./sqrt(q1+q2+2*q3./sin(2*angle));
    
    % if the angle is 0 or very close to 0, calculating sigma will fail
    if angle == 0 || sum(isnan([sigma_x,sigma_y]),2)>0;
        sigma_x = 1./sqrt(2*q1);
        sigma_y = 1./sqrt(2*q2);
        angle = 0;
    end
    
    params(5:7) = [sigma_x;sigma_y;angle];
else
    if fitSigma_y
        params(5:6) = sqrt(0.5./params(5:6));
    else
        params(5) = sqrt(0.5./params(5));
    end
end

end

