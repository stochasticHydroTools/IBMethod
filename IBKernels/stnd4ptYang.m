% ----------------------------------------------------------------------- %
% Author: Yuanxun Bill Bao
% Date: Feb 2015
% Description: A smoothed 3-point function (X.Yang et al. 2009)

function phi = stnd4ptYang(r)

r = abs(r);

phi = (r<0.5) .* (3/8 + pi/32 - r.^2/4) + ...
      (0.5<=r & r<1.5) .* (1/4+(1-r)/8.*sqrt(-2+8*r-4*r.^2)-1/8*asin(sqrt(2)*(r-1))) + ...
      (1.5<=r & r<2.5) .* (17/16-pi/64-3*r/4+r.^2/8+(r-2)/16.*sqrt(-14+16*r-4*r.^2)+1/16*asin(sqrt(2)*(r-2)));
  
