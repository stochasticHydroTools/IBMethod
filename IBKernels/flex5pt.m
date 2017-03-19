% ----------------------------------------------------------------------- %
% Author: Alex Kaiser & Bill Bao
% Date: Dec 2016
% Description: Flexible 5-pt Kernel for IBM
%
% Conditions:
% ---------------------------------------------------------------------
% | Support | Odd/Even | 1st Mom. | 2nd Mom. | 3rd Mom. | Sum Squares |
% ---------------------------------------------------------------------
% |  5-pt   |    -     |    x     |    x     |    x     |      x      |
% ---------------------------------------------------------------------
%
% Notes:
% (1) Letting K = 0 gives the so-called 'Old 6-point Kernel'.
% (2) Letting K = (38-sqrt(69))/60 gives the so-called 'New
% 5-point Kernel'. This choice of K ensures that the kernel has continuous
% first, second, and third derivatives. It is also the smallest value of K
% such that the kernel is everywhere non-negative.
% (3) This code is slow. However, the closed form of the kernel is really
% quite horrible. If speed becomes an issue, however, that will be the way
% to go, so it will be implemented eventually.
% ----------------------------------------------------------------------- %


function val = flex5pt(x,KK)

% KK = (38 - sqrt(69))/60; 
% KK = (38 + sqrt(69))/60; 

phi = @(r) (136 - 40*KK - 40*r.^2 + sqrt(2)*sqrt(3123 - 6840*KK + 3600*KK.^2 - 12440*r.^2 + 25680*KK*r.^2 - 12600*KK.^2*r.^2 + 8080*r.^4 - 8400*KK*r.^4 - 1400*r.^6))/280; 
  
% function is even, 
x = abs(x); 

if abs(x) < 0.5
    r = x; 
    val = phi(r); 
elseif abs(x) < 1.5
    r = x - 1; 
    val = (4 - 4*phi(r) - KK - 4*r + 3*KK*r - r.^2 + r.^3)/6; 
elseif abs(x) < 2.5
    r = x - 2; 
    val = (-2 + 2*phi(r) + 2*KK + r - 3*KK*r + 2*r.^2 - r.^3)/12; 
else 
    val = 0; 
end 



% Combined version 
% 
% KK = (38 - sqrt(69))/60; 
% 
% % function is even, 
% r = abs(x); 
% 
% if r < 0.5
%     val = (136 - 40*KK - 40*r.^2 + sqrt(2)*sqrt(3123 - 6840*KK + 3600*KK.^2 - 12440*r.^2 + 25680*KK*r.^2 - 12600*KK.^2*r.^2 + 8080*r.^4 - 8400*KK*r.^4 - 1400*r.^6))/280; 
% elseif r < 1.5
%     val = (324 - 10*r - 240*r.^2 + 70*r.^3 + 30*KK*(-8 + 7*r) - ... 
%           sqrt(2)*sqrt(-2637 + 960*r + 15040*r.^2 - 4320*r.^3 - 12920*r.^4 + 8400*r.^5 - 1400*r.^6 - ...  
%           1800*KK^2*(5 - 14*r + 7*r.^2) - 120*KK*(-87 + 148*r + 206*r.^2 - 280*r.^3 + 70*r.^4)))/420; 
% elseif r < 2.5
%     val = (-4 + 2*KK + (136 - 40*KK + sqrt(2)*sqrt(3123 - 6840*KK + 3600*KK^2 - 12440*(-2 + r).^2 + 25680*KK*(-2 + r).^2 - 12600*KK^2*(-2 + r).^2 + 8080*(-2 + r).^4 - ...
%        8400*KK*(-2 + r).^4 - 1400*(-2 + r).^6) - 40*(-2 + r).^2)/140 - 3*KK*(-2 + r) + 2*(-2 + r).^2 - (-2 + r).^3 + r)/12 ; 
% else 
%     val = 0; 
% end 



% With K put into the function and things simplified 
% Something funny happening in third derivative on this one 

% function is even, 
% r = abs(x); 
% 
% if r < 0.5
%     val = (332 + 2*sqrt(69) - 120*r^2 + 3*sqrt(76*(8 + sqrt(69)) - 27*(109 + 12*sqrt(69))*r^2 + 40*(138 + 7*sqrt(69))*r^4 - 2800*r^6))/840 ; 
% elseif r < 1.5
%     val = (344 + 8*sqrt(69) + (246 - 7*sqrt(69))*r - 480*r^2 + 140*r^3 - ...
%           2*sqrt(385 + 32*sqrt(69) + (606 - 472*sqrt(69))*r + 3*(-3941 + 452*sqrt(69))*r^2 - 160*(-212 + 7*sqrt(69))*r^3 + 40*(-912 + 7*sqrt(69))*r^4 + ...
%           16800*r^5 - 2800*r^6))/840; 
% elseif r < 2.5
%     val = (2340 - 18*sqrt(69) + (-2766 + 7*sqrt(69))*r + 1080*r^2 - 140*r^3 + sqrt(4*(-25511 + 815*sqrt(69)) + (372732 - 7664*sqrt(69))*r +  ... 
%            (-542463 + 6396*sqrt(69))*r^2 - 320*(-1262 + 7*sqrt(69))*r^3 + 40*(-4062 + 7*sqrt(69))*r^4 + 33600*r^5 - 2800*r^6))/1680 ; 
% else 
%     val = 0; 
% end
% 



