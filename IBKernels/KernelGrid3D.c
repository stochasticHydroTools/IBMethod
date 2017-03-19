/************************************************************************
*  KernelGrid3D.c
*  
*  MEX function to produce the 3D kernel weights w = phi(x)*phi(y)*phi(z)
*  
*  Created by yuanxun bao on 6/13/14.
*   
*************************************************************************/

/* #define char16_t UINT16_T */
#include <mex.h>
#include <string.h>
#include <math.h>
#include "Kernels.h"
#include <stdio.h>

#define w_SIZE 6
#define NDIMS 3

void outer_prod(double v1[], double v2[], double v3[], double out[], mwSize n){
  
    mwSize i,j,k;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            for (k = 0; k < n; k++)
            {
                out[i+j*n+k*(n*n)] = v1[i] * v2[j] * v3[k]; /* matlab's weird ordering */
            }
        }
    }
}


void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[]){
  
    double r1, r2, r3, K;                    /* input scalars */
    char   *kID_1, *kID_2, *kID_3;           /* input strings */
    double *outMatrix;                       /* output matirx */
    const mwSize dims[] = {w_SIZE, w_SIZE, w_SIZE};
  
    r1      = mxGetScalar(prhs[0]);
    kID_1   = mxArrayToString(prhs[1]);
    r2      = mxGetScalar(prhs[2]);
    kID_2   = mxArrayToString(prhs[3]);
    r3	    = mxGetScalar(prhs[4]);
    kID_3   = mxArrayToString(prhs[5]);
    K       = mxGetScalar(prhs[6]);
    plhs[0] = mxCreateNumericArray(NDIMS,dims,mxDOUBLE_CLASS, mxREAL);
    outMatrix = mxGetPr(plhs[0]);
  
    double v1[] = {-2-r1, -1-r1, -r1, 1-r1, 2-r1, 3-r1};
    double v2[] = {-2-r2, -1-r2, -r2, 1-r2, 2-r2, 3-r2};
    double v3[] = {-2-r3, -1-r3, -r3, 1-r3, 2-r3, 3-r3};

    phiweights(v1, kID_1, K);
    phiweights(v2, kID_2, K);
    phiweights(v3, kID_3, K);
  
    mxFree(kID_1); mxFree(kID_2); mxFree(kID_3);   /* free variables */
  
    outer_prod(v1,v2,v3, outMatrix, (mwSize)w_SIZE);
  
}
