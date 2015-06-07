%**************************************************************************
%   deltaCKplot.m  
%   Author: Charles S. Peskin
%           minor modification by Yuan-Xun Bao   
%
%   This code plots (calling deltaCK.m) the new 6-point immersed boundary 
%   kernel in [1] and its three continuous derivatives.
%   
%   Input: 
%   r: real number in [0,1]
%   K: second moment constant[1], 
%      for the new 6-pt kernel, K = 59/60 - sqrt(29)/20;
%      for the standard 6-pt kernel, K = 0    
%
%   Notation:
%   p   = phi
%   dp  = derivative of phi
%   d2p = second derivative of phi
%   d3p = third derivative of phi
%
%   m3 = "minus 3", that is evaluate at r-3
%   m2 = "minus 2", that is evaluate at r-2
%   m1 = "minus 1", that is evaluate at r-1
%
%   p1 = "plus 1", that is evaluate at r+1
%   p2 = "plus 2", that is evaluate at r+2
%   
%   [1] "A gaussian-like immersed boundary kernel with three
%   continuous derivatives and improved translational invariance" 
%   http://arxiv.org/abs/1505.07529
%
%**************************************************************************


clear all
clf
r=0:.0001:1;
K= 59/60 - sqrt(29)/20;

[pm3 pm2 pm1 p pp1 pp2 dpm3 dpm2 dpm1 dp dpp1 dpp2 d2pm3 d2pm2 d2pm1 d2p d2pp1 d2pp2 d3pm3 d3pm2 d3pm1 d3p d3pp1 d3pp2 flag]=deltaCK(r,K);
    if(~flag)
      figure(1)
      plot(r-3,  pm3,r-2,  pm2,r-1,  pm1,r,  p,r+1,  pp1,r+2,  pp2)
      title('new 6-pt IB kernel')
      
      figure(2)
      hold on;
      plot(r-3, dpm3,r-2, dpm2,r-1, dpm1,r, dp,r+1, dpp1,r+2, dpp2)
      plot(r-3,d2pm3,r-2,d2pm2,r-1,d2pm1,r,d2p,r+1,d2pp1,r+2,d2pp2)
      plot(r-3,d3pm3,r-2,d3pm2,r-1,d3pm1,r,d3p,r+1,d3pp1,r+2,d3pp2)
      title('derivatives of the new 6-pt IB kernel')
      hold off
    end

%check that defining properties are satisfied:

sum_odd  = pm3 + pm1 + pp1;
sum_even = pm2 + p   + pp2;
moment1 = (r-3).*pm3 + (r-2).*pm2 + (r-1).*pm1 + r.*p + (r+1).*pp1 + (r+2).*pp2;
moment2 = ((r-3).^2).*pm3 + ((r-2).^2).*pm2 + ((r-1).^2).*pm1 + (r.^2).*p + ((r+1).^2).*pp1 + ((r+2).^2).*pp2;
moment3 = ((r-3).^3).*pm3 + ((r-2).^3).*pm2 + ((r-1).^3).*pm1 + (r.^3).*p + ((r+1).^3).*pp1 + ((r+2).^3).*pp2;
sum_of_squares = pm3.^2 + pm2.^2 + pm1.^2 + p.^2 + pp1.^2 + pp2.^2;

err_sum_odd = max(abs(sum_odd - (1/2)))
err_sum_even= max(abs(sum_even - (1/2)))
err_moment1 = max(abs(moment1))
err_moment2 = max(abs(moment2-K))
err_moment3 = max(abs(moment3))
err_sum_of_squares = max(sum_of_squares)-min(sum_of_squares)

