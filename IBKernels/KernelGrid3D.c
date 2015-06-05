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



void phiweights(double r[], char *kernel, double K)
{
  
    mwSize i;
  
    if (strcmp(kernel, "bspline4pt") == 0)
    {	
        for (i=1; i<5; i++)
        {
            r[i] = bspline4pt(r[i]);
        }
        r[0] = 0.; r[5] = 0.;   /* 4-pt kernel, 1st, last are zero */
    }
    else if (strcmp(kernel, "bspline4pt_d") == 0)
    {	
        for (i=1; i<5; i++)
        {
            r[i] = bspline4pt_d(r[i]);
        }
        r[0] = 0.; r[5] = 0.;   /* 4-pt kernel, 1st, last are zero */
    }
    else if (strcmp(kernel, "stnd4pt") == 0)
    {
        for (i=1; i<5; i++)
        {
            r[i] = stnd4pt(r[i]);
        }
        r[0] = 0.; r[5] = 0.;
    }
    else if (strcmp(kernel, "stnd3pt") == 0)
    {
        double R = -r[2];  /* for odd number kernel, we need to check r<1/2 or r>1/2 */
        if (R<=0.5)
        {
            for(i=1; i<4; i++)
            {
                r[i] = stnd3pt(r[i]);
            }
            r[0]=r[4]=r[5]=0.;
        }
        else
        {
            for (i=2; i<5; i++)
            {
                r[i] = stnd3pt(r[i]); 
            }
            r[0]=r[1]=r[5]=0.;
        }
       
    }
    else if (strcmp(kernel, "flex6pt") == 0)
    {	
        for (i=0; i<6; i++)
        {
            r[i] = flex6pt(r[i],K);
        }	
    }
    else if (strcmp(kernel, "flex6pt_d") == 0)
    {
        for (i=0; i<6; i++)
        {
            r[i] = flex6pt_d(r[i],K);
        }
    }
}

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
