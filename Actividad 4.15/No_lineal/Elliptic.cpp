//	###############################################################################################################
// Evaluation of elliptic integrals and functions, using approximations of Fenton and Gardener-Garden for m > 1/2
//	###############################################################################################################

#define	Subprograms
#include "Cnoidal.h"

void Elliptic_integrals(void)
{
int i;

m1 = 1.-m;
m4 = pow(m,0.25);

if(m1 > m1_limit)	K  =  2./pow(1+m4,2)*log(2*(1+m4)/(1-m4));

Kd = 2*pi/pow(1.+m4,2);

q1 = exp(-pi*K/Kd);
e = (2.-m)/3.+pi/2/K/Kd+2.*pow(pi/Kd,2)*(-1./24.+q1*q1/pow(1-q1*q1,2));

ee[1] = e;
mm[1] = m;

for(i=2; i<=5 ; ++i)
	{
	ee[i] = ee[i-1]*e;
	mm[i] = mm[i-1]*m;
	}

return;
}

//**********************************************
// Elliptic functions
//**********************************************

// Bring back to range [-2K,+2K]

double Shift(double u)
{
int N;
double uu;
N = trunc(u/(4*K));
uu = u-N*4*K;
if (uu < -2*K) uu+=4*K;
if (uu > 2*K) uu -= 4*K;
return uu;
}

double sn(double u)
{
double factor, term, w;
if(m1 > m1_limit)	factor = pow(m,-0.25);
else					factor = 1.;
w = 0.5*pi*Shift(u)/Kd;
term = factor*(sinh(w)-q1*q1*sinh(3*w))/(cosh(w)+q1*q1*cosh(3*w));
return(term);
}

double cn(double u)
{
double factor, term, w;
if(m1 > m1_limit)	factor = 0.5*pow(m1/m/q1,0.25);
else					factor = 1.;
w = 0.5*pi*Shift(u)/Kd;
term = factor*(1-2*q1*cosh(2*w))/(cosh(w)+q1*q1*cosh(3*w));
return(term);
}

double dn(double u)
{
double factor, term, w;
if(m1 > m1_limit)	factor = 0.5*pow(m1/q1,0.25);
else					factor = 1.;
w = 0.5*pi*Shift(u)/Kd;
term = factor*(1+2*q1*cosh(2*w))/(cosh(w)+q1*q1*cosh(3*w));
return(term);
}
