% ----------------------------------------------------------------------- %
% Author: Jason Kaye, Yuanxun Bill Bao
% Date: June 2013
% Description: Spread force from grid markers to Cartesian grid
%
% Inputs:
% F - Vector force on markers.
% X - Marker locations [x,y].
% N - Domain size [Nx,Ny].
% h - Cartesian grid spacing.
% kernelID - Can be 'stnd3pt', 'stnd4pt', 'bspline6pt', etc.
% K - Parameter for flexible kernels.
% mexFlag - 1 use mex, 0 use phiweights.m
%
% Reference: This routine is only a slight modification of that
% 'spreadMAC2D' written by Yuanxun Bill Bao on May 20, 2013, altered
% primarily to add flexibility in the choice of kernel.
% ----------------------------------------------------------------------- %

function f = spreadMAC2Dvector(F,X,N,h,kernelID,K,mexFlag)

Nx  = N(1); Ny  = N(2);
Nb  = size(X,1);
f = zeros(Nx,Ny,2);
c = 1/(h*h);

for k = 1 : Nb
    % g1 grid
    s  = [X(k,1)/h, X(k,2)/h-1/2];
    p  = floor(s);
    r  = s-p;
    j1 = mod((p(1)-2):(p(1)+3),Nx)+1;
    j2 = mod((p(2)-2):(p(2)+3),Ny)+1;
    
    if mexFlag == 0
        w  = phiweights(r(1),kernelID,K)*phiweights(r(2),kernelID,K).';
    else
        w  = KernelGrid2D(r(1), kernelID  , r(2), kernelID, K);
    end
    f(j1,j2,1) = f(j1,j2,1) + c*(F(k,1)*w);
    
    % g2 grid
    s  = [X(k,1)/h-1/2, X(k,2)/h];
    p  = floor(s);
    r  = s-p;
    j1 = mod((p(1)-2):(p(1)+3),Nx)+1;
    j2 = mod((p(2)-2):(p(2)+3),Ny)+1;
    
    if mexFlag == 0
        w  = phiweights(r(1),kernelID,K)*phiweights(r(2),kernelID,K).';
    else    
        w  = KernelGrid2D(r(1), kernelID  , r(2), kernelID, K);
    end
    f(j1,j2,2) = f(j1,j2,2) + c*(F(k,2)*w);    
end