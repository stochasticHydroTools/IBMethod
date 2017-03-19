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
#include <string.h>
#include <mex.h>

#define sgn(x) ( (x > 0) ? 1 : ((x < 0) ? -1 : 0) )


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

/* derivative of the standard 4-pt kernel */
double stnd4pt_d(double r){

    int s = (r>0) ? 1 : -1;    
    double r2 = r*r;
    double phi;   
 
    if (fabs(r) < 1)
    {
        phi = ((2*s-4*r)/sqrt(1+4*fabs(r)-4*r2)-2*s) / 8.; 
    }
    else if (fabs(r)>=1 && fabs(r)<2)
    {
        phi = (-(6*s-4*r)/sqrt(-7+12*fabs(r)-4*r2)-2*s) / 8.;
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


/* the 6-pt B-spline kernel */
double bspline6pt(double r){
  
    r = fabs(r);
    double r2 = r*r;
    double r3 = r2*r;
    double r4 = r2*r2;
    double r5 = r2*r3;
    double phi;
    
    if (r <= 1)
    {
       phi = 11./20-1./2*r2+1./4*r4-1./12*r5;
    }
    else if (r >1 && r<=2)
    {
       phi = 17./40+5./8*r-7./4*r2+5./4*r3-3./8*r4+1./24*r5;
    }
    else if (r > 2 && r<=3)
    {
       phi = 81./40-27./8*r+9./4*r2-3./4*r3+1./8*r4-1./120*r5;
    }
    else 
    {
        phi = 0;
    }
    
    return phi;
}


/* derivative of the 6-pt B-spline kernel */
double bspline6pt_d(double r){
  
    double phi;
    double r2 = r*r;
    double r3 = r2*r;
    double r4 = r2*r2;

    if (r<=-2 && r>-3)
    {
        phi = 27./8+9./2*r+9./4*r2+1./2*r3+5./120*r4;
    }
    else if (r<=-1 && r>-2)
    {
        phi = -5./8-7./2*r-15./4*r2-3./2*r3-5./24*r4;
    }
    else if (r<=0 && r>-1)
    {
        phi = -r+r3 + 5./12*r4;
    }
    else if (r<=1 && r>0)
    {
        phi = -r+r3 - 5./12*r4;
    }
    else if (r<=2 && r>1)
    {
        phi = 5./8-7./2*r+15./4*r2-3./2*r3+5./24*r4;
    }
    else if (r<=3 && r>2)
    {
        phi = -27./8+9./2*r-9./4*r2+1./2*r3-5./120*r4;
    }
    else
    {
        phi = 0;
    }
  
    return phi;
}

/*5pt C^3 kernel */
double flex5pt(double x, double K)
{
    
    double r,r2,r3,r4,r6,K2,phi,val;
    K2=K*K;
    
    if (fabs(x) < 0.5){ r = fabs(x); }
    else if (fabs (x) < 1.5){ r = fabs(x)-1; }
    else if (fabs (x) < 2.5){ r = fabs(x)-2; }
    else { return 0; };
    
    r2=r*r; 
    r3=r*r2;
    r4=r2*r2;
    r6=r3*r3;
    
    phi=(136 - 40*K - 40*r2 + sqrt(2)*sqrt(3123 - 6840*K + 3600*K2 - 12440*r2 + 25680*K*r2 - 12600*K2*r2 + 8080*r4 - 8400*K*r4 - 1400*r6))/280.;
    
    if (fabs(x) < 0.5) 
    { 
        val = phi; 
    }
    else if (fabs (x) < 1.5)
    { 
        val = (4 - 4*phi - K - 4*r + 3*K*r - r2 + r3)/6; 
    }
    else if (fabs (x) < 2.5)
    { 
        val = (-2 + 2*phi + 2*K + r - 3*K*r + 2*r2 - r3)/12;    
    }
    
    return val;   
}

/* derivative of C^3 5-pt kernel */
double flex5pt_d(double x, double K)
{
    
    double r,r2,r3,r4,r5,r6,K2,val,dphi,disc;
    K2=K*K;
    
    if (fabs(x) < 0.5){ r = fabs(x); }
    else if (fabs (x) < 1.5){ r = fabs(x)-1; }
    else if (fabs (x) < 2.5){ r = fabs(x)-2; }
    else { return 0; };
    
    r2=r*r; 
    r3=r*r2;
    r4=r2*r2;
    r5=r2*r3;
    r6=r3*r3;
    
    disc = sqrt(3123 - 6840*K + 3600*K2 - 12440*r2 + 25680*K*r2 - 12600*K2*r2 + 8080*r4 - 8400*K*r4 - 1400*r6);
    dphi = ( -40*2*r + sqrt(2)*(1./2./disc)*(-24880*r + 25680*K*2*r - 12600*K2*2*r + 8080*4*r3 - 8400*K*4*r3 - 1400*6*r5) )/280.;  
    
    if (fabs(x) < 0.5) 
    { 
        val = dphi*sgn(x); 
    }
    else if (fabs (x) < 1.5)
    { 
        val = (-4*dphi - 4 + 3*K - 2*r + 3*r2)/6. * sgn(x); 
    }
    else if (fabs (x) < 2.5)
    { 
        val = ( 2*dphi + 1 - 3*K + 4*r - 3*r2)/12. * sgn(x);    
    }
    
    return val;   
}


/*************************************************************************/
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
    else if (strcmp(kernel, "bspline6pt") == 0)
    {	
        for (i=0; i<6; i++)
        {
            r[i] = bspline6pt(r[i]);
        }
    }
    else if (strcmp(kernel, "bspline6pt_d") == 0)
    {	
        for (i=0; i<6; i++)
        {
            r[i] = bspline6pt_d(r[i]);
        }
    }
    
    else if (strcmp(kernel, "stnd4pt") == 0)
    {
        for (i=1; i<5; i++)
        {
            r[i] = stnd4pt(r[i]);
        }
        r[0] = 0.; r[5] = 0.;
    }
    else if (strcmp(kernel, "stnd4pt_d") == 0)
    {
        for (i=1; i<5; i++)
        {
            r[i] = stnd4pt_d(r[i]);
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
    else if (strcmp(kernel, "flex5pt") == 0)
    {
        double R = -r[2];  /* for odd number kernel, we need to check r<1/2 or r>1/2 */
        if (R <= 0.5)
        {
            for(i=0; i<5; i++)
            {
                r[i] = flex5pt(r[i],K);
            }
            r[5]=0.;
            
        } else {
            for (i=1; i<6; i++)
            {
                r[i] = flex5pt(r[i],K); 
            }
            r[0]=0.;            
        }
        
    }
    else if (strcmp(kernel, "flex5pt_d") == 0)
    {
        double R = -r[2];  /* for odd number kernel, we need to check r<1/2 or r>1/2 */
        if (R <= 0.5)
        {
            for(i=0; i<5; i++)
            {
                r[i] = flex5pt_d(r[i],K);
            }
            r[5]=0.;
            
        } else {
            for (i=1; i<6; i++)
            {
                r[i] = flex5pt_d(r[i],K); 
            }
            r[0]=0.;            
        }
        
    }
}
