% ----------------------------------------------------------------------- %
% Author: Jason Kaye
% Date: June 2013
% Description: 6-point B-Spline Kernel for IBM
%
% Conditions:
% ---------------------------------------------------------------------
% | Support | Odd/Even | 1st Mom. | 2nd Mom. | 3rd Mom. | Sum Squares |
% ---------------------------------------------------------------------
% |  6-pt   |          |    x     |    x     |    x     |             |
% ---------------------------------------------------------------------
% |---------|
% | 4th Mom.|
% |---------|
% |    x    |
% |---------|
%
% Note: 2nd and 4th moment constants have not yet been analytically
% computed.
% ----------------------------------------------------------------------- %

function phi=bspline6pt(r)

r = abs(r);
phi =...
    (r<1).*(11/20-1/2*r.^2+1/4*r.^4-1/12*r.^5) +...
    (1<=r&r<2).*(17/40+5/8*r-7/4*r.^2+5/4*r.^3-3/8*r.^4+1/24*r.^5) +...
    (2<=r&r<3).*(81/40-27/8*r+9/4*r.^2-3/4*r.^3+1/8*r.^4-1/120*r.^5);

end