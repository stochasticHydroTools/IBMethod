% ----------------------------------------------------------------------- %
% Title: PoissonSolver2D
% Author: Jason Kaye
% Date: June 2013
%
% Description: Solves Poisson's equation with periodic boundary conditions
% in 2D using Fourier basis
%
% Inputs:
% f - Right hand side to Poisson's equation -\Delta u = f, defined on 2D
% grid
% N - Number of grid points for each coordinate, given in the form [Nx,Ny]
% h - Grid point spacing
% L_hat - Discrete Laplacian operator.
%
% Outputs:
% u - Solution defined on same grid as f
%
% Note: Code for discrete Laplacian given by
% [m2,m1] = meshgrid(0:Ny-1, 0:Nx-1);
% L_hat   = -4/(h*h) * (sin(pi/Nx*m1).^2 + sin(pi/Ny*m2).^2);
% L_hat(1,1) = 1; % Arbitrary choice to avoid divide by zero

% ----------------------------------------------------------------------- %

function u = PoissonSolver2D(f,L_hat)

% Fourier transform of f
f_hat = fft2(f);

% Solution in Fourier space
u_hat = f_hat./L_hat;

% Inverse Fourier transform
u = -real(ifft2(u_hat));

end
