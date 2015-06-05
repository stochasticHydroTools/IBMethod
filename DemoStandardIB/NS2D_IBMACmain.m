clear; restoredefaultpath;
addpath('../IBKernels/')

% profile on;
    
L=1; 
nu=0.1;

%Kernel = {'flex6pt','flex6pt_d', (59/60)*(1-sqrt(261/3481))};
%Kernel = {'stnd4pt','stnd4pt_d', []};
%Kernel = {'stnd3pt','stnd3pt_d', []};
Kernel={'bspline4pt', 'bspline4pt_d',[]};
%Kernel = {'new4ptC3','new4ptC3_d', []};
%Kernel = {'new4ptC2','new4ptC2_d', []};

showplot = 'on';
mexFlag = 1;

% Eulerian grid
Nx=64; Ny=Nx;
N=[Nx,Ny];
h=L/Nx;

% Lagrangian grid
alpha = 5/28; beta = 7/20;
%alpha=1/4; beta=1/4;
MpC = 3; % number markers per cell
Ns = round(MpC*(2*pi*alpha)/L*Nx);
ds=2*pi/Ns;
s =(0:Ns-1)*ds;
X = [alpha*cos(s'), beta*sin(s')]*L+L/2;
%alpha = 1/5; beta = 1/5;
%X = [alpha*cos(s').*(1+0.05*cos(4*s')), beta*(1+0.05*cos(4*s')).*sin(s')]*L+L/2;

% initial velocity
u=zeros(Nx,Ny,2);

% time step
tend = 1;
dt   = h/15;
Nt   = floor(tend/dt);
dt   = tend/Nt;
tt   = 0:dt:tend;
Nf   = Nt;

tic
figure(1); clf;
%[uIBMAC,XIBMAC,areaIBMAC]=NS2D_IBMAC(L,N,nu,Kernel,X,u,tend,dt,Nf,MpC,showplot);
%save(['NS2D-IBMAC-circle-spacing-',num2str(MpC),'-',Kernel{1},'.mat'],'uIBMAC','XIBMAC','areaIBMAC','tt');
NS2D_IBMAC(L,N,nu,Kernel,X,u,tend,dt,Nf,MpC,showplot,mexFlag);
toc
% profile viewer
% p = profile('info');
% profsave(p,'profile_results')
% profile off;
