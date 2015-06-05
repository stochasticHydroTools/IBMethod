% ----------------------------------------------------------------------- %
% NS2D_IBMAC.m
% Yuanxun Bill Bao
% Nov 09, 2013
% 
% Description: solve the fluid-structure interaction of an elastic
% membrane by solving a 2D NS flow. The standard immersed 
% boundary on a staggered grid is used (IB-MAC). 
% AB2 for time stepping
% RK2 to get another initial value
%
% Inputs:
% L      -- domain length
% N      -- grid size N=[Nx,Ny], Nx=Ny for now
% nu     -- viscosity
% Kernel -- [KernelID, DerivativeID, K]
% X0     -- initial configuration of the memberane
% u0     -- intial velocity
% tend   -- final time
% dt     -- time step
% ----------------------------------------------------------------------- %

function NS2D_IBMAC(X,u)

global L mu rho Kernel showplot mexFlag
global Nx Ny Ns ds h dt Nt Nf

N=[Nx,Ny];
x=(0:Nx-1)*h;
y=(0:Ny-1)*h;
[yy,xx]=meshgrid(y,x);
xx2=xx+h/2; yy2=yy+h/2;
sk=1:1:Nx; % skip

% set Kernel
kID=Kernel{1}; K=Kernel{3};

% 2D Laplacian in Fourier space
[m2,m1] = meshgrid(fftshift(-Ny/2:Ny/2-1), fftshift(-Nx/2:Nx/2-1));
L_hat1  = -4/(h*h) * (sin(pi/N(1)*m1).^2 + sin(pi/N(2)*m2).^2);
L_hat2  = L_hat1; L_hat2(1,1)=1;


%% use 2nd-order RK scheme to get u(1*dt), X(1*dt)
uold=u; Xold=X; ADVold = advection2D(u,N,h);
% compute X(n+1/2)
XX=X+(dt/2)*interpMAC2Dvector(u,X,N,h,kID,K,mexFlag);
% force-spreading f(n+1/2)
FF=Force(XX,Ns,ds)*ds;
ff=spreadMAC2Dvector(FF,XX,N,h,kID,K,mexFlag);
% fluid solver, Navier Stokes, RK (2 step)
[u,uu]=NavierStokes2D_FFT_RK(u,N,h,dt,mu,rho,ff,L_hat1,L_hat2);
ADV=advection2D(u,N,h);
% compute X(n+1)
X=X+dt*interpMAC2Dvector(uu,XX,N,h,kID,K,mexFlag);

%% time stepping (AB2 for advection)
for n=2:Nt
    
    % compute X(n+1)
    Uold=interpMAC2Dvector(uold,Xold,N,h,kID,K,mexFlag);
    U   =interpMAC2Dvector(u,X,N,h,kID,K,mexFlag);
    UU  = 3/2*U - 1/2*Uold;
    Xold = X;
    X=X+dt*UU;
    
    % compute f(n+1/2)
    XX=(X+Xold)/2;
    FF=Force(XX,Ns,ds)*ds;
    ff=spreadMAC2Dvector(FF,XX,N,h,kID,K,mexFlag);
    
    % fluid solve
    S=3/2*ADV-1/2*ADVold;
    u = NavierStokes2D_FFT(u,N,h,dt,mu,rho,S,ff,L_hat1,L_hat2);
    
    ADVold = ADV;
    ADV=advection2D(u,N,h);

    % plotting
    figure(1)
    if mod(n,Nf) == 0
        disp([num2str(n/Nt*100),'%'])
        
        if strcmp(showplot, 'on')
                    
            uu = (u([2:Nx,1],:,1) + u(:,:,1))/2;
            vv = (u(:,[2:Nx,1],2) + u(:,:,2))/2;
            quiver(xx2(sk,sk),yy2(sk,sk),uu(sk,sk),vv(sk,sk),1); hold on;
            plot(mod(X(:,1),L),mod(X(:,2),L),'k-','linewidth',2);    
            xlabel(['t = ',num2str(n*dt),', N = ',num2str(Nx)])
            title(['NS, IB-MAC, ',kID])
            axis equal
            axis([0,L,0,L])
            %axis([.2,.5,.5,.8]);
            %axis([.25,.35,.6,.7])
            drawnow
            hold off

        end
    end    
end


