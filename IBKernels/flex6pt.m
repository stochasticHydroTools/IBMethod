% ----------------------------------------------------------------------- %
%   Author: Jason Kaye
%           modified by Yuan-Xun Bao
%   Date: June 2013
%   Description: Flexible 6-pt Kernel for IBM
%
%   Conditions:
%   ---------------------------------------------------------------------
%   | Support | Odd/Even | 1st Mom. | 2nd Mom. | 3rd Mom. | Sum Squares |
%   ---------------------------------------------------------------------
%   |  6-pt   |    x     |    x     |    x     |    x     |      x      |
%   ---------------------------------------------------------------------
%
%   Notes:
%   (1) Letting K = 0 gives the so-called 'standard 6-point Kernel'.
%   (2) Letting K = 59/60-sqrt(29)/20 gives the so-called 'New 6-point 
%   Kernel'. This choice of K ensures that the kernel has three continuous
%   derivatives. It is also the smallest value of K such that the kernel 
%   is everywhere non-negative.
%
%   [1] "A gaussian-like immersed boundary kernel with three
%   continuous derivatives and improved translational invariance" 
%   http://arxiv.org/abs/1505.07529
% ----------------------------------------------------------------------- %

function phi=flex6pt(r,K)

R = r;

% compute r in [0,1]
r = (-3<r & r<=-2).*(r+3) + (-2<r & r<=-1).*(r+2) + ...
    (-1<r & r<=0) .*(r+1) + (0<r & r<=1)  .*(r+0) + ...
    (1<r & r<=2)  .*(r-1) + (2<r & r<=3)  .*(r-2);

alpha=28;
beta=(9/4)-(3/2)*(K+r.^2)+((22/3)-7*K)*r-(7/3)*r.^3;
gamma=(1/4)*( ((161/36)-(59/6)*K+5*K^2)*(1/2)*r.^2 + (-(109/24)+5*K)*(1/3)*r.^4 + (5/18)*r.^6 );

discr=beta.^2-4*alpha*gamma;
% flag=0;
% if(any(discr<0))
%    flag=1
% end

pm3=(-beta+sign((3/2)-K)*sqrt(discr))/(2*alpha);
pm2= -3*pm3 - (1/16) + (1/8)*(K+r.^2) + (1/12)*(3*K-1)*r + (1/12)*r.^3;
pm1=  2*pm3 + (1/4)                   +  (1/6)*(4-3*K)*r -  (1/6)*r.^3;
p  =  2*pm3 + (5/8)  - (1/4)*(K+r.^2);
pp1= -3*pm3 + (1/4)                   -  (1/6)*(4-3*K)*r +  (1/6)*r.^3;
pp2=    pm3 - (1/16) + (1/8)*(K+r.^2) - (1/12)*(3*K-1)*r - (1/12)*r.^3;

phi = (-3<R & R<=-2).* pm3 + ...
      (-2<R & R<=-1).* pm2 + ...
      (-1<R & R<=0) .* pm1 + ...
      (0<R & R<=1)  .* p   + ...
      (1<R & R<=2)  .* pp1 + ...
      (2<R & R<=3)  .* pp2;

end