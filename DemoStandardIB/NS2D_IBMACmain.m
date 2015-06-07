% *************************************************************************
%   NS2D_IBMACmain.m
%   Author: Yuanxun Bill Bao
%   July 2014
%   
%   L:        domain [0,L]x[0,L]
%   mu:       viscosity
%   rho:      density
%   Kernel:   choice of IB kernel
%   showplot: plotting 'on' or 'off'
%   mexFlag:  '1' use mex excutables, '0' use matlab routines (slower)
%   Nx,Ny:    Eulerian grid size
%   h:        grid spacing
%   Ns:       Lagrangian grid size
%   ds:       Lagrangian grid spacing
%   X:        initial configuration
%   u:        initial velocity
%   tend:     final time
%   dt:       timestep
%   Nt:       number of timesteps
% *************************************************************************

clear; restoredefaultpath;
addpath('../IBKernels/')

% profile on

global L mu rho Kernel showplot mexFlag
global Nx Ny Ns ds h dt Nt Nf

L=1; 
mu=0.01;
rho=1;

Kernel = {'flex6pt','', 59/60-sqrt(29)/20};
%Kernel = {'stnd4pt','', []};
%Kernel = {'stnd3pt','', []};
%Kernel={'bspline4pt','',[]};

showplot = 'on';
mexFlag = 1;

% Eulerian grid size
Nx=64; Ny=Nx;
h=L/Nx;

% Lagrangian grid
Ns = 4*Nx;
ds=2*pi/Ns;
s =(0:Ns-1)*ds;

% initial configuration: ellipse
alpha = 5/28; beta = 7/20;
X = [alpha*cos(s'), beta*sin(s')]*L+L/2;

% initial velocity
u=zeros(Nx,Ny,2);

% time step
tend = 2;
dt   = h/2;
Nt   = floor(tend/dt);
dt   = tend/Nt;
Nf   = 2; % plot every Nf timesteps

tic
NS2D_IBMAC(X,u);
toc

% profile viewer
% p = profile('info');
% profsave(p,'profile_results')
% profile off;
