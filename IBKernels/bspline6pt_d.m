function phi=bspline6pt_d(r)

s=sign(r);
phi = ... 
      (abs(r)<=1) .* (-r+r.^3 - 5/12*r.^4.*s )+ ...
      (abs(r)>1 & abs(r)<=2) .* (5/8*s-7/2*r+15/4*r.^2.*s-3/2*r.^3+5/24*r.^4.*s) + ...
      (abs(r)>2 & abs(r)<=3) .* (-27/8*s+9/2*r-9/4*r.^2.*s+1/2*r.^3-5/120*r.^4.*s);
  
end
    
