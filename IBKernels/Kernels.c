/**********************************************************
*   Kernels.c
* 
*   Discrete delta kernels/derivatives for the 
*   Immersed Boundary Method 
*
*   Created by yuanxun bao on 6/15/14.
*
**********************************************************/

#include <stdio.h>
#include <math.h>

/* the new C3 6-pt kernel */
double flex6pt(double r, double K){
	
    double  alpha, beta, gamma, discr, R, R2, R3; 
	int sgn;
	double phi;
	
    R = r - ceil(r) + 1;  /* R between [0,1] */
    R2 = R * R; R3 = R2*R;
    alpha = 28.;
    beta  = 9./4 - 1.5 * (K + R2) + (22./3-7*K)*R - 7./3*R3;
    gamma = 0.25 * ( 0.5*(161./36 - 59./6*K + 5*K*K)*R2 + 1./3*(-109./24 + 5*K)*R2*R2 + 5./18*R3*R3  );
    discr = beta*beta - 4 * alpha * gamma;
    
    sgn = ((1.5-K)>0) ? 1 : -1;   /* sign(3/2 - K) */
  
    if (-3 < r && r <= -2)
    {
        phi =  1./(2*alpha) * ( -beta + sgn * sqrt(discr) );
    }
    else if (-2 < r && r <=-1 )
    {
		phi = -3./(2*alpha) * ( -beta + sgn * sqrt(discr) ) -
		       1./16 + 1./8*( K+(r+2)*(r+2) ) + 1./12*(3*K-1)*(r+2) + 1./12*(r+2)*(r+2)*(r+2); 
    }
    else if (-1 < r && r <= 0 )
    {	
		phi = 2./(2*alpha) * ( -beta + sgn * sqrt(discr) ) + 
		      1./4 + 1./6*(4-3*K)*(r+1) - 1./6*(r+1)*(r+1)*(r+1);	      
    }
    else if ( 0 < r && r <= 1 )
    {	
		phi = 2./(2*alpha) * ( -beta + sgn * sqrt(discr) ) +
		      5./8 - 1./4 * ( K+r*r );
    }
    else if ( 1 < r && r <= 2 )
    {	
		phi = -3./(2*alpha) * ( -beta + sgn * sqrt(discr) ) + 
		      1./4 - 1./6*(4-3*K)*(r-1) + 1./6*(r-1)*(r-1)*(r-1);
    }
    else if ( 2 < r && r <= 3 )
    {	
		phi = 1./(2*alpha) * ( -beta + sgn * sqrt(discr) ) - 
		      1./16 + 1./8*(K+(r-2)*(r-2)) - 1./12*(3*K-1)*(r-2) - 1./12*(r-2)*(r-2)*(r-2); 	
    }
    else
    {
	    phi = 0.;
    }
	
	return phi;
  	
}


/* derivative of the new C3 6-pt kernel */
double flex6pt_d(double r, double K){
	
    double alpha, beta, dbeta, gamma, dgamma, discr; 
    double pm3, dpm3; 
    double R, R2, R3; 
    int sgn = ((1.5-K)>0) ? 1 : -1;
    double phi;
	
	
    R = r - ceil(r) + 1;    /* shift to between [0,1] */
    R2 = R * R; R3 = R2*R;
    alpha = 28.;
  
    beta  = 9./4 - 1.5 * (K + R2) + (22./3-7*K)*R - 7./3*R3;
    dbeta =                 -3.*R + (22./3-7*K)   - 7*R2   ;
  
    gamma = 0.25 * ( 0.5*(161./36 - 59./6*K + 5*K*K)*R2 + 1./3*(-109./24 + 5*K)*R2*R2 + 5./18*R3*R3 );
    dgamma= 0.25 * (     (161./36 - 59./6*K + 5*K*K)*R  + 4./3*(-109./24 + 5*K)*R3    + 5./3 *R3*R2 );
  
    discr = beta*beta - 4 * alpha * gamma;
  
    pm3  =  (-beta + sgn*sqrt(discr)) / (2*alpha);
    dpm3 = -(dbeta*pm3+dgamma) / (2*alpha*pm3+beta);
  
    if ( -3 < r && r <= -2 )
    {	
        phi = dpm3;
    }
    else if( -2 < r && r <= -1 )
    {  
	    phi = -3*dpm3 + 1./12*(3*K-1) + 1./4*R + 1./4*R2;
    }
    else if( -1 < r && r <=  0 )
    {
		phi =  2*dpm3 +  1./6*(4-3*K)          - 1./2*R2;	
    }
    else if(  0 < r && r <=  1 )
    {	
		phi =  2*dpm3                 - 1./2*R; 	
    }
    else if(  1 < r && r <=  2 )
    {
		phi = -3*dpm3 - 1./6 *(4-3*K)          + 1./2*R2;
    }
    else if(  2 < r && r <=  3 )
    {	
		phi =    dpm3 - 1./12*(3*K-1) + 1./4*R - 1./4*R2;	
    }
    else
    {	
		phi = 0.;
    } 
  
    return phi;
	
}

/* the 4-pt B-spline kernel */
double bspline4pt(double r){
  
    r = fabs(r);
    double r2 = r*r;
    double r3 = r2*r;
    double phi;
    
    if (r <= 1.)
    {
       phi = 2./3. - r2 + .5*r3;
    }
    else if (r >1. && r<2.)
    {
       phi = 4./3 - 2*r + r2 - r3/6.;
    }
    else
    {
       phi = 0.;
    }
    
    return phi;
}

/* derivative of the 4-pt B-spline kernel */
double bspline4pt_d(double r){
  
    double phi;
    double r2 = r*r;

    if (r<=-1 && r>-2)
    {
        phi = 2 + 2*r + .5*r2;
    }
    else if (r<=0 && r>-1)
    {
        phi = -2*r - 1.5*r2;
    }
    else if (r<=1 && r>0)
    {
        phi = -2*r + 1.5*r2;
    }
    else if (r<=2 && r>1)
    {
        phi = -2 + 2*r - .5*r2;
    }
    else
    {
        phi = 0.;
    }
  
    return phi;
}

/* the standard 4-pt kernel */
double stnd4pt(double r){

    r = fabs(r);
    double r2 = r*r;
    double phi;   
 
    if (r < 1)
    {
        phi = (3 - 2*r + sqrt(1+4*r-4*r2)) / 8.; 
    }
    else if (r>=1 && r<2)
    {
        phi = (5 - 2*r - sqrt(-7+12*r-4*r2)) / 8.;
    }
    else
    {
        phi = 0.;
    }
    
    return phi;
}

/* the standard 3-pt kernel */
double stnd3pt(double r){

    r = fabs(r);
    double r2 = r*r;
    double phi;
    
    if (r<0.5)
    {
        phi = (1 + sqrt(1-3*r2)) / 3.;
    }
    else if (r>=0.5 && r<1.5)
    {
        phi = (5 - 3*r - sqrt(-3*(1-r)*(1-r)+1)) / 6.;
    }
    else 
    {
        phi = 0.;
    }
    
    return phi;

}
