% ----------------------------------------------------------------------- %
% NavierStokes2D_FFT
% Yuanxun Bill Bao
% July 1, 2014
%
% Description: FFT based 2D Navier Stokes solver
% 
% Inputs:
% N      -- grid size N=[Nx,Ny], Nx=Ny for now
% mu     -- viscosity
% rho    -- density
% S      -- advection term
% f      -- forcing
% L_hat1 -- Laplacian in Fourier space
% L_hat2 -- Laplacian in Fourier space with zeroth mode fixed to be 1
% ----------------------------------------------------------------------- %


function unew = NavierStokes2D_FFT(u,N,h,dt,mu,rho,S,f,L_hat1,L_hat2)

Nx = N(1); Ny = N(2);

c1 = dt/rho; c2 = mu*dt / (2*rho);

% poisson solve for pressure p(n+1/2)
w1=(-rho/dt)*u(:,:,1)+rho*S(:,:,1)-mu/2*Laplacian2D(u(:,:,1),N,h)-f(:,:,1);
w2=(-rho/dt)*u(:,:,2)+rho*S(:,:,2)-mu/2*Laplacian2D(u(:,:,2),N,h)-f(:,:,2);
Divw = (w1([2:Nx,1],:)-w1)/h + (w2(:,[2:Ny,1])-w2)/h;
pp = PoissonSolver2D(Divw,L_hat2);


% solve for velocity u(n+1)
r1 = u(:,:,1) - dt*S(:,:,1) + c2*Laplacian2D(u(:,:,1),N,h) + ...
     c1*f(:,:,1) - c1*(pp-pp([Nx,1:Nx-1],:))/h;
r2 = u(:,:,2) - dt*S(:,:,2) + c2*Laplacian2D(u(:,:,2),N,h) + ...
     c1*f(:,:,2) - c1*(pp-pp(:,[Ny,1:Ny-1]))/h;
u1_hat = fft2(r1)./(1-c2*L_hat1);
u2_hat = fft2(r2)./(1-c2*L_hat1);

unew(:,:,1) = real(ifft2(u1_hat));
unew(:,:,2) = real(ifft2(u2_hat));




