% ----------------------------------------------------------------------- %
% Author: Yuanxun Bill Bao
% Date: Feb 2015
% Description: A smoothed 3-point function (X.Yang et al. 2009)

function phi = stnd3ptYang(r)

r = abs(r);

phi = (r<1).*         (17/48 + sqrt(3)*pi/108 +r/4-r.^2/4 +(1-2*r)/16.*sqrt(-12*r.^2+12*r+1) - ...
                       sqrt(3)/12*asin(sqrt(3)/2*(2*r-1))) + ...
      (1<=r & r<2) .* (55/48 - sqrt(3)*pi/108 -13*r/12+r.^2/4 + (2*r-3)/48.*sqrt(-12*r.^2+36*r-23) + ...
                       sqrt(3)/36*asin(sqrt(3)/2*(2*r-3)));

