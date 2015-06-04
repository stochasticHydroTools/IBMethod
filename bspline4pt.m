% ----------------------------------------------------------------------- %
% Author: Jason Kaye
% Date: June 2013
% Description: 4-point B-Spline Kernel for IBM
%
% Conditions:
% ---------------------------------------------------------------------
% | Support | Odd/Even | 1st Mom. | 2nd Mom. | 3rd Mom. | Sum Squares |
% ---------------------------------------------------------------------
% |  4-pt   |          |    x     | K = 1/3  |    x     |             |
% ---------------------------------------------------------------------
%
% ----------------------------------------------------------------------- %

function phi=bspline4pt(r)

r = abs(r);
phi =...
    (r<=1).*(2/3-r.^2+(1/2)*r.^3)+...
    (1<r & r<2).*(4/3-2*r+r.^2-(1/6)*r.^3);

end
