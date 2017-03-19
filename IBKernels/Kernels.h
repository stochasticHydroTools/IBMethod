/******************************************************
*	Kernels.h
* 
*
*	Created by yuanxun bao on 6/15/14.
******************************************************/

#ifndef _Kernels_h
#define _Kernels_h

double bspline4pt( double );

double bspline4pt_d( double );

double bspline6pt( double );

double bspline6pt_d( double );

double flex6pt( double, double );

double flex6pt_d( double, double );

double flex5pt( double, double );

double flex5pt_d( double, double );

double stnd4pt( double );

double stnd4pt_d( double );

double stnd3pt( double );

void phiweights(double[], char[], double);

#endif
