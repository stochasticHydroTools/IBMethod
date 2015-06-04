% ************************************************************************
%   bspline4pt_d.m
%   Author: Yuan-Xun Bao
%   June 2014   
%
%   1st derivative of the 4pt B-spline kernel
%   can be used for interpolating derivative of a discrete delta function
%   in the divergence-free IB method
% *************************************************************************

function phi=bspline4pt_d(r)

s=sign(r);
phi =...
    (abs(r)<=1).*(-2*r+3/2*r.^2.*s)+...
    (1<abs(r) & abs(r)<=2).*(-2*s+2*r-1/2*r.^2.*s);

end