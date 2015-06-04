% *************************************************************************
%   Title: phiweights.m
%   Author: Jason Kaye
%           minor modification by Yuan-Xun Bao
% 
%   Date: June 2013
%   Description: Compute one-dimensional kernel weights at Cartesian grid
%   points for kernels with up to 6-point support.
%   
%   Inputs:
%   r - Point at which to compute weights.
%   kernelID - can be 'stnd3pt', 'stnd4pt', 'bspline4pt', etc. 
%              or 1s derivative of a kernel, 'bspline4pt_d', etc.
%   
%   K - Parameter for flexible families, 
%       default: K = [],
%       new 6pt: K = 59/60 - sqrt(29)/20
%       standard 6pt: K = 0
%
%   Outputs:
%   w - Column vector of weights at Cartesian grid point.
% *************************************************************************

function w=phiweights(r,kernelID,K)

w = [0;0;0;0;0;0];

switch kernelID
    case 'bspline4pt'
          w(2) = bspline4pt(-1-r);
          w(3) = bspline4pt(-r);
          w(4) = bspline4pt(1-r);
          w(5) = bspline4pt(2-r);
    case 'bspline4pt_d'
         w(2) = bspline4pt_d(-1-r);
         w(3) = bspline4pt_d(-r);
         w(4) = bspline4pt_d(1-r);
         w(5) = bspline4pt_d(2-r);
    case 'flex6pt'
         w(1) = flex6pt(-2-r,K);
         w(2) = flex6pt(-1-r,K);
         w(3) = flex6pt(-r,K);
         w(4) = flex6pt(1-r,K);
         w(5) = flex6pt(2-r,K);
         w(6) = flex6pt(3-r,K);
    case 'flex6pt_d'
         w(1) = flex6pt_d(-2-r,K);
         w(2) = flex6pt_d(-1-r,K);
         w(3) = flex6pt_d(-r,K);
         w(4) = flex6pt_d(1-r,K);
         w(5) = flex6pt_d(2-r,K);
         w(6) = flex6pt_d(3-r,K);
    case 'stnd3pt'
        if r <= 1/2
            w(2) = stnd3pt(-1-r);
            w(3) = stnd3pt(-r);
            w(4) = stnd3pt(1-r);
        elseif r>1/2 && r<=1
            w(3) = stnd3pt(-r);
            w(4) = stnd3pt(1-r);
            w(5) = stnd3pt(2-r);
        end  
    case 'stnd4pt'
        w(2) = stnd4pt(-1-r);
        w(3) = stnd4pt(-r);
        w(4) = stnd4pt(1-r);
        w(5) = stnd4pt(2-r);
    case 'stnd6pt'
        K = 0;
        w(1) = flex6pt(-2-r,K);
        w(2) = flex6pt(-1-r,K);
        w(3) = flex6pt(-r,K);
        w(4) = flex6pt(1-r,K);
        w(5) = flex6pt(2-r,K);
        w(6) = flex6pt(3-r,K);        
    case 'new6pt'
        K = 59/60-sqrt(29)/20;
        w(1) = flex6pt(-2-r,K);
        w(2) = flex6pt(-1-r,K);
        w(3) = flex6pt(-r,K);
        w(4) = flex6pt(1-r,K);
        w(5) = flex6pt(2-r,K);
        w(6) = flex6pt(3-r,K);
    otherwise
        error('Not a valid kernelID.');
end

end
