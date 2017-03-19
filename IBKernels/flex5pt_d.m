function val = flex5pt_d(x,KK)

% KK = (38 - sqrt(69))/60; 
% KK = (38 + sqrt(69))/60; 

disc = @(r) sqrt(3123 - 6840*KK + 3600*KK.^2 - 12440*r.^2 + 25680*KK*r.^2 - 12600*KK.^2*r.^2 + 8080*r.^4 - 8400*KK*r.^4 - 1400*r.^6);

% phi = @(r) (136 - 40*KK - 40*r.^2 + sqrt(2)*disr(r))/280; 
dphi = @(r) ( -40*2*r + sqrt(2)*(1/2/disc(r))*(-24880*r + 25680*KK*2*r - 12600*KK.^2*2*r + 8080*4*r.^3 - 8400*KK*4*r.^3 - 1400*6*r.^5) )/280;  
% function is even, 
s = sign(x);
x = abs(x); 
if abs(x) < 0.5
    r = x; 
    val = dphi(r)*s; 
elseif abs(x) < 1.5
    r = x - 1; 
    val = (-4*dphi(r) - 4 + 3*KK - 2*r + 3*r.^2)/6.*s; 
elseif abs(x) < 2.5
    r = x - 2; 
    val = ( 2*dphi(r) + 1 - 3*KK + 4*r - 3*r.^2)/12.*s; 
else 3
    val = 0; 
end 