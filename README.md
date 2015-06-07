# IBMethod
The IBMethod repository contains codes developed for the Immersed Boundary Method.

IBKernels
* contains codes for computing IB kernels and routines for computing weights used in force spreading and velocity interpolation.
* deltaCK.m deltaCKplot.m are Matlab codes demonstrating the new Gaussian-like 6pt IB kernel
* Matlab codes for IB kernels are stnd3pt.m, stnd4pt.m, bspline4pt.m, flex6pt.m. First derivatives of the 4pt B-spline and the new 6pt kernels are also available in bspline4pt_d.m and flex6pt_d.m. The C version of these kernels are available in Kernels.c
* phiweights.m computes weights of spreading/interpolation in 1D. To get weights in 2D/3D, one needs to compute a tensor product. The C version is available in KernelGrid2D.c and KernelGrid3D.c. No tensor product is required as the C codes already take care of it. 
* *mexa64, *mexamaci64 are MEX executables of KernelGrid2D.c and KernelGrid3D.c on Linux and OSX platforms. One can use them in a Matlab program to speed up computing IB weights. 

DemoStandardIB
* contains codes for a demo of a 2D elastic membrane immersed in a viscous incompressible fluid.
* NS2D_IBMACmain.m is the main program. Users can set fluid parameters, grid resolutions, choice of IB kernel and other parameters in this code. 
* NS2D_IBMAC.m is called by NS2D_IBMACmain.m to carry out the main steps of the simulation.
* spreadMAC2Dvector.m and interpMAC2Dvector.m: force spreading and velocity interpolation to a MAC grid. Users have the option to use MEX excutables to speed up calcuation of weights. 
* NavierStokes2D_FFT.m: a FFT-based 2D Navier-Stokes fluid solver, used in the main loop of time stepping
* NavierStokes2D_FFT_RK.m: a FFT-based 2D Navier-Stokes fluid solver, used only in the beginning to get another initial conditiong via RK2
