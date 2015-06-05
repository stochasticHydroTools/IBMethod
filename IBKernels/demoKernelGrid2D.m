% *************************************************************************
%   Title: demoKernelGrid2D.m
%   Author: Yuan-Xun Bao
%      
%   Description: a short demo testing MEX implementation of 2D IB weights
%
% *************************************************************************

addpath('../../../IBStaggered/KernelsTesting/Kernels');

K=[];
% check standard 3pt kernel
kID1='stnd3pt';
kID2='stnd3pt';
r = rand(1,2);  % r in [0,1]
wx = phiweights(r(1),kID1,K);
wy = phiweights(r(2),kID2,K);
w = wx*wy.';
y = KernelGrid2D(r(1),kID1,r(2),kID2,K);
disp('check stnd3pt:')
disp(['max diff = ', num2str(max(max(abs(w-y))))]);
disp(' ');

% check standard 4pt kernel
kID1='stnd4pt';
kID2='stnd4pt';
r = rand(1,2);  % r in [0,1]
wx = phiweights(r(1),kID1,K);
wy = phiweights(r(2),kID2,K);
w = wx*wy.';
y = KernelGrid2D(r(1),kID1,r(2),kID2,K);
disp('check stnd4pt:')
disp(['max diff = ', num2str(max(max(abs(w-y))))]);
disp(' ');

% check standard 6pt kernel
K = 0;
kID1='flex6pt';
kID2='flex6pt';
r = rand(1,2);  % r in [0,1]
wx = phiweights(r(1),kID1,K);
wy = phiweights(r(2),kID2,K);
w = wx*wy.';
y =  KernelGrid2D(r(1),kID1,r(2),kID2,K);
disp('check stnd6pt:')
disp(['max diff = ', num2str(max(max(abs(w-y))))]);
disp(' ');


% check bspline 4pt and its derivative
K = 0;
kID1='bspline4pt';
kID2='bspline4pt_d';
r = rand(1,2);  % r in [0,1]
wx = phiweights(r(1),kID1,K);
wy = phiweights(r(2),kID2,K);
w = wx*wy.';
y = KernelGrid2D(r(1),kID1,r(2),kID2,K);
disp('check bspline4pt and its derivative:')
disp(['max diff = ', num2str(max(max(abs(w-y))))]);
disp(' ');

% check standard 6pt kernel
K = 59/60-sqrt(29)/20;
kID1='flex6pt';
kID2='flex6pt';
r = rand(1,2);  % r in [0,1]
wx = phiweights(r(1),kID1,K);
wy = phiweights(r(2),kID2,K);
w = wx*wy.';
y = KernelGrid2D(r(1),kID1,r(2),kID2,K);
disp('check new 6pt and its derivative:')
disp(['max diff = ', num2str(max(max(abs(w-y))))]);
disp(' ');
