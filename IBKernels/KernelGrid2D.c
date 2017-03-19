/********************************************************************
*  KernelGrid2D.c
*  
*  MEX function to produce the 2D kernel weights w = phi(x)*phi(y)
*  
*  Created by yuanxun bao on 6/13/14.
*   
*********************************************************************/

/* #define char16_t UINT16_T */
#include <mex.h>
#include <string.h>
#include <math.h>
#include "Kernels.h"
#include <stdio.h>

#define w_SIZE 6


void outer_prod(double v1[], double v2[], double out[], mwSize n){

    mwSize i,j;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            out[i+j*n] = v1[i] * v2[j]; /* matlab's weird ordering */
        }
    }
}


void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[]){
  
    double r1, r2, K;                        /* input scalars */
    char   *kID_1, *kID_2;                /* input strings */
    double *outMatrix;                    /* output matirx */
  
    r1      = mxGetScalar(prhs[0]);
    kID_1   = mxArrayToString(prhs[1]);
    r2      = mxGetScalar(prhs[2]); 
    kID_2   = mxArrayToString(prhs[3]);
    K       = mxGetScalar(prhs[4]);
    plhs[0] = mxCreateDoubleMatrix((mwSize)w_SIZE, (mwSize)w_SIZE, mxREAL);

    outMatrix = mxGetPr(plhs[0]);
  
    double v1[] = {-2-r1, -1-r1, -r1, 1-r1, 2-r1, 3-r1};
    double v2[] = {-2-r2, -1-r2, -r2, 1-r2, 2-r2, 3-r2};
  
    phiweights(v1, kID_1, K);
    phiweights(v2, kID_2, K);
  
    mxFree(kID_1); mxFree(kID_2);   /* free variables */
  
    outer_prod(v1,v2, outMatrix, (mwSize)w_SIZE);
}
