#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <complex.h> 
#include <string.h> 

void main()
{

double vh,vc,gmix,g0,gex,Whc;
void PSEoS( double ,double )

/* EOS After Pitzer and Sterner, 1994 - API, The Journal of Chemical Physics */
if (strcmp( name, "H2O") == 0 || strcmp( name, "CO2") == 0 ){

	double p_bar = 1000.*P; //in bar
	double c1, c2, c3, c4, c5, c6, c7, c8, c9, c10; 

	if (strcmp( name, "H2O") == 0){
		c1  =  0.24657688e6 / T + 0.51359951e2;
		c2  =  0.58638965e0 / T - 0.28646939e-2 + 0.31375577e-4 * T;
		c3  = -0.62783840e1 / T + 0.14791599e-1 + 0.35779579e-3 * T +  0.15432925e-7 * pow(T,2.0);
		c4  = -0.42719875e0 - 0.16325155e-4 * T;
		c5  =  0.56654978e4 / T - 0.16580167e2 + 0.76560762e-1 * T;
		c6  =  0.10917883e0;
		c7  =  0.38878656e13 / pow(T,4.0) - 0.13494878e9 / pow(T,2.0) + 0.30916564e6 / T + 0.75591105e1;
		c8  = -0.65537898e5 / T + 0.18810675e3;
		c9  = -0.14182435e14 / pow(T,4.0) + 0.18165390e9 / pow(T,2.0) - 0.19769068e6 / T - 0.23530318e2;
		c10 =  0.92093375e5 / T + 0.12246777e3;
	}
	else {	// can only be CO2
		c1  =  0.18261340e7 / T + 	0.79224365e2;
		c2  =  						0.66560660e-4 	+ 0.57152798e-5 * T + 0.30222363e-9 * pow(T,2.0);
		c3 	= 						0.59957845e-2 	+ 0.71669631e-4 * T + 0.62416103e-8 * pow(T,2.0);
		c4  = -0.13270279e1 / T +  -0.15210731e0  	+ 0.53654244e-3 * T - 0.71115142e-7 * pow(T,2.0);
		c5  =  0.12456776e0 / T +   0.49045367e1    + 0.98220560e-2 * T + 0.55962121e-5 * pow(T,2.0);
		c6  = 				     	0.75522299e0;
		c7  = -0.39344644e+12 / pow(T,4.0) + 0.90918237e8 / pow(T,2.0) + 0.42776716e6 / T - 0.22347856e2;
		c8 	=  0.40282608e3 / T +   0.11971627e3;
		c9  =  0.22995650e8 / pow(T,2.0) - 0.78971817e5 / T - 0.63376456e2;
		c10 =  0.95029765e5 / T + 0.18038071e2;
	}

	
	/* solve for volume at P, T */
	int    err,  k;
	double vsub, yr;
	double R1     = 83.144;
	double data[] = {R1,T,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,p_bar};
	
	double x1     = 3.0;
	double x2     = R1*T/P;
	
	double e      = 1e-14;
	int maxiter   = 500;
	int mode      = 0;               														/** Mode is used to send the right *data (see root_finding.c) */
	
	vsub          =  BrentRoots(x1,x2,data,e,mode,maxiter, &yr, &k, &err);
	
	double r      =   1.0/vsub;
	double Ares   =   R1*T*( c1*r + (1.0/(c2 + c3*r + c4*pow(r, 2.0) + c5*pow(r, 3.0) + c6*pow(r, 4.0)) - 1.0/c2) - c7/c8*(exp(-c8*r) - 1.0) - c9/c10*(exp(-c10*r) - 1.0) );
	vterm         =   (Ares + p_bar*vsub + R1*T*(log( R1*T / vsub ) - 1.0)) * 1e-4;	
}



}
