% ----------------------------------------------------------------------- %
% NavierStokes2D_FFT_RK
% Yuanxun Bill Bao
% June, 2015
%
% Description: FFT based 2D Navier Stokes solver for Runge-Kutta (2 step)
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


function [unew,uu] = NavierStokes2D_FFT_RK(u,N,h,dt,mu,rho,f,L_hat1,L_hat2)

Nx=N(1); Ny=N(2);
xp=[2:Nx,1];
xm=[Nx,1:Nx-1];
yp=[2:Ny,1];
ym=[Ny,1:Ny-1];

%% solve for pressure
ADV=advection2D(u,N,h);
w1=(-2*rho/dt)*u(:,:,1)+rho*ADV(:,:,1)-f(:,:,1);
w2=(-2*rho/dt)*u(:,:,2)+rho*ADV(:,:,2)-f(:,:,2);
Divw=(w1(xp,:)-w1)/h + (w2(:,yp)-w2)/h;
pp=PoissonSolver2D(Divw,L_hat2);

%% solve for uu=u(n+1/2)
c1=dt/(2*rho); c2=mu*dt/(2*rho);
r1=u(:,:,1)-(dt/2)*ADV(:,:,1)+c1*f(:,:,1)-c1*(pp-pp(xm,:))/h;
r2=u(:,:,2)-(dt/2)*ADV(:,:,2)+c1*f(:,:,2)-c1*(pp-pp(:,ym))/h;
uu1_hat=fft2(r1)./(1-c2*L_hat1);
uu2_hat=fft2(r2)./(1-c2*L_hat1);
uu(:,:,1)=real(ifft2(uu1_hat));
uu(:,:,2)=real(ifft2(uu2_hat));

%% solve for pp=p(n+1/2)
ADV=advection2D(uu,N,h);
w1=(-rho/dt)*u(:,:,1)+rho*ADV(:,:,1)-mu/2*Laplacian2D(u(:,:,1),N,h)-f(:,:,1);
w2=(-rho/dt)*u(:,:,2)+rho*ADV(:,:,2)-mu/2*Laplacian2D(u(:,:,2),N,h)-f(:,:,2);
Divw=(w1(xp,:)-w1)/h + (w2(:,yp)-w2)/h;
pp=PoissonSolver2D(Divw,L_hat2);

%% solve for unew=u(n+1/2)
c3=dt/rho;
r1=u(:,:,1)-dt*ADV(:,:,1)+c2*Laplacian2D(u(:,:,1),N,h)+c3*f(:,:,1)-c3*(pp-pp(xm,:))/h;
r2=u(:,:,2)-dt*ADV(:,:,2)+c2*Laplacian2D(u(:,:,2),N,h)+c3*f(:,:,2)-c3*(pp-pp(:,ym))/h;
u1_hat=fft2(r1)./(1-c2*L_hat1);
u2_hat=fft2(r2)./(1-c2*L_hat1);
unew(:,:,1)=real(ifft2(u1_hat));
unew(:,:,2)=real(ifft2(u2_hat));
