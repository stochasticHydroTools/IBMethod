% ----------------------------------------------------------------------- %
% Author: Jason Kaye
% Date: June 2013
% Description: Standard 3-point Kernel for IBM
%
% Conditions:
% ---------------------------------------------------------------------
% | Support | Odd/Even | 1st Mom. | 2nd Mom. | 3rd Mom. | Sum Squares |
% ---------------------------------------------------------------------
% |  3-pt   |          |     x    |          |          |   C = 1/2   |
% ---------------------------------------------------------------------
%
% ----------------------------------------------------------------------- %

function phi=stnd3pt(r)

r = abs(r);
phi =...
    (r<1/2).*((1/3)*(1+sqrt(1-3*r.^2))) +...
    (1/2<=r & r<3/2).*((1/6)*(5-3*r-sqrt(-3*(1-r).^2+1)));

end