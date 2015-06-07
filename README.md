# IBMethod
The IBMethod repository contains codes developed for the Immersed Boundary Method.

IBKernels
* contains codes for computing IB kernels and routines for computing weights used in force spreading and velocity interpolation.
* Matlab codes for IB kernels are stnd3pt.m, stnd4pt.m, bspline4pt.m, flex6pt.m. First derivatives of the 4pt B-spline and the new 6pt kernels are also available in bspline4pt_d.m and flex6pt_d.m. The C version of these kernels are available in Kernels.c
* phiweights.m computes weights of spreading/interpolation in 1D. To get weights in 2D/3D, one needs to compute a tensor product. The C version is available in KernelGrid2D.c and KernelGrid3D.c. No tensor product is required as the C codes already take care of it. 
* *mexa64, *mexamaci64 are MEX executables of KernelGrid2D.c and KernelGrid3D.c on Linux and OSX platforms. One can use them in a Matlab program to speed up computing IB weights. 

DemoStandardIB
