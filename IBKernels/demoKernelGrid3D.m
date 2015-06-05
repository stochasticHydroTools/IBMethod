% *************************************************************************
%   Title: demoKernelGrid3D.m
%   Author: Yuan-Xun Bao
%      
%   Description: a short demo testing MEX implementation of 3D IB weights
%
% *************************************************************************


K=[];
% check standard 3pt kernel
kID1='stnd3pt';
kID2='stnd3pt';
kID3='stnd3pt';
r = rand(1,3);  % r in [0,1]
wx = repmat(phiweights(r(1),kID1,K),[1 6 6]);
wy = repmat(phiweights(r(2),kID2,K).',[6 1 6]);
wz = repmat(reshape(phiweights(r(3),kID3,K),1,1,[]),[6 6 1]);
w = wx.*wy.*wz;
y =  KernelGrid3D(r(1),kID1,r(2),kID2,r(3),kID3,K);
disp('check stnd3pt:')
disp(['max diff = ', num2str(max(max(max(abs(w-y)))))]);
disp(' ');

% check standard 4pt kernel
kID1='stnd4pt';
kID2='stnd4pt';
kID3='stnd4pt';
r = rand(1,3);  % r in [0,1]
wx = repmat(phiweights(r(1),kID1,K),[1 6 6]);
wy = repmat(phiweights(r(2),kID2,K).',[6 1 6]);
wz = repmat(reshape(phiweights(r(3),kID3,K),1,1,[]),[6 6 1]);
w = wx.*wy.*wz;
y =  KernelGrid3D(r(1),kID1,r(2),kID2,r(3),kID3,K);
disp('check stnd4pt:')
disp(['max diff = ', num2str(max(max(max(abs(w-y)))))]);
disp(' ');


% check 4pt bspline and its derivative
kID1='bspline4pt';
kID2='bspline4pt_d';
kID3='bspline4pt';
r = rand(1,3);  % r in [0,1]
wx = repmat(phiweights(r(1),kID1,K),[1 6 6]);
wy = repmat(phiweights(r(2),kID2,K).',[6 1 6]);
wz = repmat(reshape(phiweights(r(3),kID3,K),1,1,[]),[6 6 1]);
w = wx.*wy.*wz;
y =  KernelGrid3D(r(1),kID1,r(2),kID2,r(3),kID3,K);
disp('check bspline4t and bspline4pt_d:')
disp(['max diff = ', num2str(max(max(max(abs(w-y)))))]);
disp(' ');


% check standard 6pt kernel
K = 0;
kID1='flex6pt';
kID2='flex6pt';
kID3='flex6pt';
r = rand(1,3);  % r in [0,1]
wx = repmat(phiweights(r(1),kID1,K),[1 6 6]);
wy = repmat(phiweights(r(2),kID2,K).',[6 1 6]);
wz = repmat(reshape(phiweights(r(3),kID3,K),1,1,[]),[6 6 1]);
w = wx.*wy.*wz;
y =  KernelGrid3D(r(1),kID1,r(2),kID2,r(3),kID3,K);
disp('check stnd6pt:')
disp(['max diff = ', num2str(max(max(max(abs(w-y)))))]);
disp(' ');

% check new 6pt kernel and its derivative
K = 59/60 - sqrt(29)/20;
kID1='flex6pt_d';
kID2='flex6pt';
kID3='flex6pt';
r = rand(1,3);  % r in [0,1]
wx = repmat(phiweights(r(1),kID1,K),[1 6 6]);
wy = repmat(phiweights(r(2),kID2,K).',[6 1 6]);
wz = repmat(reshape(phiweights(r(3),kID3,K),1,1,[]),[6 6 1]);
w = wx.*wy.*wz;
y =  KernelGrid3D(r(1),kID1,r(2),kID2,r(3),kID3,K);
disp('check new 6pt and its derivative:')
disp(['max diff = ', num2str(max(max(max(abs(w-y)))))]);
disp(' ');
