% ************************************************************************
%   flex6pt_d.m
%   Author: Yuan-Xun Bao
%   June 2014   
%
%   1st derivative of the new 6pt kernel
%   can be used for interpolating derivative of a discrete delta function
%   in the divergence-free IB method
%
%   [1] "A gaussian-like immersed boundary kernel with three
%   continuous derivatives and improved translational invariance" 
%   http://arxiv.org/abs/1505.07529
% *************************************************************************


function dphi=flex6pt_d(r,K)

R=r;

% compute r in [0,1]
r = (-3<r & r<=-2).*(r+3) + (-2<r & r<=-1).*(r+2) + ...
    (-1<r & r<=0) .*(r+1) + (0<r & r<=1)  .*(r+0) + ...
    (1<r & r<=2)  .*(r-1) + (2<r & r<=3)  .*(r-2);
    

alpha=28;

beta=(9/4)-(3/2)*(K+r.^2)+((22/3)-7*K)*r-(7/3)*r.^3;
dbeta=((22/3)-7*K)-3*r-7*r.^2;

gamma=(1/4)*( ((161/36)-(59/6)*K+5*K^2)*(1/2)*r.^2 + (-(109/24)+5*K)*(1/3)*r.^4 + (5/18)*r.^6 );
dgamma=(1/4)*(((161/36)-(59/6)*K+5*K^2)*r + (-(109/24)+5*K)*(4/3)*r.^3 + (5/3)*r.^5);

discr=beta.^2-4*alpha*gamma;
% flag=0;
% if(any(discr<0))
%    flag=1
% end

pm3=(-beta+sign((3/2)-K)*sqrt(discr))/(2*alpha);

dpm3= -(dbeta.*pm3+dgamma)./(2*alpha*pm3+beta);
dpm2= -3*dpm3 + (1/12)*(3*K-1) + (1/4)*r + (1/4)*r.^2;
dpm1=  2*dpm3 +  (1/6)*(4-3*K)           - (1/2)*r.^2;
dp  =  2*dpm3                  - (1/2)*r;
dpp1= -3*dpm3 -  (1/6)*(4-3*K)           + (1/2)*r.^2;
dpp2=    dpm3 - (1/12)*(3*K-1) + (1/4)*r - (1/4)*r.^2;

dphi = (-3<R & R<=-2).* dpm3 + ...
       (-2<R & R<=-1).* dpm2 + ...
       (-1<R & R<=0) .* dpm1 + ...
       (0<R & R<=1)  .* dp   + ...
       (1<R & R<=2)  .* dpp1 + ...
       (2<R & R<=3)  .* dpp2;
       
       
       
       
