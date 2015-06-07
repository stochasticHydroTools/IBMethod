%**************************************************************************
%   deltaCK.m  
%   Author: Charles S. Peskin
%           minor modification by Yuan-Xun Bao   
%
%   This code produces (called by deltaCKplot.m) the new 6-point 
%   immersed boundary kernel in [1] and its three continuous 
%   derivatives.
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

function [pm3,pm2,pm1,p,pp1,pp2,dpm3,dpm2,dpm1,dp,...
          dpp1,dpp2,d2pm3,d2pm2,d2pm1,d2p,d2pp1,d2pp2,d3pm3,d3pm2,...
          d3pm1,d3p,d3pp1,d3pp2,flag] = deltaCK(r,K)

alpha=28;

beta=(9/4)-(3/2)*(K+r.^2)+((22/3)-7*K)*r-(7/3)*r.^3;
dbeta=((22/3)-7*K)-3*r-7*r.^2;
d2beta=-3-14*r;
d3beta=  -14;

gamma=(1/4)*( ((161/36)-(59/6)*K+5*K^2)*(1/2)*r.^2 + (-(109/24)+5*K)*(1/3)*r.^4 + (5/18)*r.^6 );
dgamma=(1/4)*(((161/36)-(59/6)*K+5*K^2)*r + (-(109/24)+5*K)*(4/3)*r.^3 + (5/3)*r.^5);
d2gamma=(1/4)*( ((161/36)-(59/6)*K+5*K^2) + (-(109/24)+5*K)*4*r.^2 + (25/3)*r.^4);
d3gamma= (-(109/24)+5*K)*2*r +(25/3)*r.^3;

discr=beta.^2-4*alpha*gamma;
flag=0;
if(any(discr<0))
   flag=1
end

pm3=(-beta+sign((3/2)-K)*sqrt(discr))/(2*alpha);
pm2= -3*pm3 - (1/16) + (1/8)*(K+r.^2) + (1/12)*(3*K-1)*r + (1/12)*r.^3;
pm1=  2*pm3 + (1/4)                   +  (1/6)*(4-3*K)*r -  (1/6)*r.^3;
p  =  2*pm3 + (5/8)  - (1/4)*(K+r.^2);
pp1= -3*pm3 + (1/4)                   -  (1/6)*(4-3*K)*r +  (1/6)*r.^3;
pp2=    pm3 - (1/16) + (1/8)*(K+r.^2) - (1/12)*(3*K-1)*r - (1/12)*r.^3;

dpm3= -(dbeta.*pm3+dgamma)./(2*alpha*pm3+beta);
dpm2= -3*dpm3 + (1/12)*(3*K-1) + (1/4)*r + (1/4)*r.^2;
dpm1=  2*dpm3 +  (1/6)*(4-3*K)           - (1/2)*r.^2;
dp  =  2*dpm3                  - (1/2)*r;
dpp1= -3*dpm3 -  (1/6)*(4-3*K)           + (1/2)*r.^2;
dpp2=    dpm3 - (1/12)*(3*K-1) + (1/4)*r - (1/4)*r.^2;

d2pm3= - (2*alpha*dpm3.^2 + 2*dbeta.*dpm3 + d2beta.*pm3 + d2gamma)./(2*alpha*pm3+beta);
d2pm2= -3*d2pm3 + (1/4) + (1/2)*r;
d2pm1=  2*d2pm3         -       r;
d2p  =  2*d2pm3 - (1/2);
d2pp1= -3*d2pm3         +       r;
d2pp2=    d2pm3 + (1/4) - (1/2)*r;

d3pm3= - (6*alpha*dpm3.*d2pm3 + d3beta.*pm3 + 3*d2beta.*dpm3 +3*dbeta.*d2pm3 + d3gamma) ./ (2*alpha*pm3 + beta);
d3pm2= -3*d3pm3 + (1/2);
d3pm1=  2*d3pm3 - 1;
d3p  =  2*d3pm3;
d3pp1= -3*d3pm3 + 1;
d3pp2=    d3pm3 - (1/2);


