% ----------------------------------------------------------------------- %
% Author: Jason Kaye, Yuanxun Bill Bao
% Date: June 2013
% Description: Interpolate velocity from Cartesian grid to markers
%
% Inputs:
% u - Vector fluid velocity.
% X - Marker locations [x,y].
% N - Domain size [Nx,Ny].
% h - Cartesian grid spacing.
% kernelID - Can be 'stnd3pt', 'stnd4pt', 'bspline6pt', etc.
% K - Parameter for flexible kernels.
%
% Reference: This routine is only a slight modification of that
% 'interpMAC2D' written by Yuanxun Bill Bao on May 20, 2013, altered
% primarily to add flexibility in the choice of kernel.
% ----------------------------------------------------------------------- %

function U = interpMAC2Dvector(u,X,N,h,kernelID,K)

Nx = N(1); Ny = N(2);
Nb = size(X,1);  
U  = zeros(Nb,2);

for k = 1 : Nb
    % g1 grid
    s  = [X(k,1)/h, X(k,2)/h-1/2];      
    p  = floor(s);                  
    r  = s-p;
    j1 = mod((p(1)-2):(p(1)+3),Nx)+1;   
    j2 = mod((p(2)-2):(p(2)+3),Ny)+1;
    %w  = phiweights(r(1),kernelID,K)*phiweights(r(2),kernelID,K).';
    w  = KernelGrid2D(r(1), kernelID  , r(2), kernelID, K);
    
    U(k,1) = sum(sum(w.*u(j1,j2,1)));
    
    % g2 grid
    s  = [X(k,1)/h-1/2, X(k,2)/h];
    p  = floor(s);
    r  = s-p;
    j1 = mod((p(1)-2):(p(1)+3),Nx)+1;
    j2 = mod((p(2)-2):(p(2)+3),Ny)+1;
    % w  = phiweights(r(1),kernelID,K)*phiweights(r(2),kernelID,K).';
    w  = KernelGrid2D(r(1), kernelID  , r(2), kernelID, K);
    U(k,2) = sum(sum(w.*u(j1,j2,2)));
  
end
