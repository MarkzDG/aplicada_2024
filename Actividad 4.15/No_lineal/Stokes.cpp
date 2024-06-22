
// Stokes theory calculations

#include <math.h>
#include <stdio.h>
#include <process.h>
#include <string.h>
#include <conio.h>
#include <stdlib.h>
#define	ANSI
#include "../Allocation.h"

#define	Main
#define	Char	char
#define	Int		int
#define	Double	double
#include "../Allocation.h"
#include "../Headers.h"


Double
	kH, skd, ckd, tkd, SU;
Double
	ss[9], t[9], C[6], D[6], E[6], e[6];

char Diagname[30], Theory[10];

// Main program

int main(void)
{
int	i, Read_data(void), iter, Iter_limit=40;

double	F(double), kd1, kd2, kFM, omega, delta, accuracy=1.e-6, F1, F2, Fd;
void	CDE(double), AB(void), Output(void);

Input1 = fopen("../Data.dat","r");
strcpy(Convergence_file,"../Convergence.dat");
strcpy(Points_file,"../Points.dat");
monitor = stdout;
strcpy(Theory,"Stokes");
strcpy(Diagname,"../Catalogue.res");

Read_data();

z = dvector(0,2*n+10);
Y = dvector(0,n);
B = dvector(0,n);
coeff = dvector(0,n);
Tanh = dvector(0,n);

monitor = stdout;

H = MaxH;

iff(Case,Wavelength)
	{
	kd = 2. * pi / L;
	kH = kd * H ;
	CDE(kd);
	}

// If period is specified, solve dispersion relation using secant method

// Until February 2015 the bisection method was used for this.
// I found that in an extreme case (large current) the bracketting
// of the solution was not correct, and the program failed,
// without printing out a proper error message.

iff(Case,Period)
	{
	fprintf(monitor,"\n# Period has been specified.\n# Now solving for L/d iteratively, printing to check convergence\n\n");
	omega = 2*pi/T;
	// Fenton & McKee for initial estimate
	kFM = omega*omega*pow(1/tanh(pow(omega,1.5)),(2./3.));
	kd1=kFM;
	kd2=kFM*1.01;
	CDE(kd2);
	F2=F(kd2);
	for(iter=1; iter<=Iter_limit; ++iter)
		{
		CDE(kd1);
		F1=F(kd1);
		Fd=(F2-F1)/(kd2-kd1);
		delta=F1/Fd;
		kd2=kd1;
		kd1=kd1-delta;
		fprintf(monitor,"%8.4f\n",2*pi/kd1);
		if (fabs(delta/kd1) < accuracy) break;
		F2=F1;
		if (iter >= Iter_limit)
			{
			printf("\n\nSecant for solution of wavenumber has not converged");
			printf("\nContact John Fenton johndfenton@gmail.com");
			getch();
			exit(1);
			}
		}
	kd=kd1;
	kH = kd * H ;
	}

z[1] = kd;
z[2] = kH;

SU = 0.5*kH/pow(kd,3);
printf("\n# Stokes-Ursell no.: %7.3f",SU);
if( SU > 0.5)
	printf(" > 1/2. Results are unreliable");
else
	printf(" < 1/2, Stokes theory should be valid");

e[1] = 0.5 * kH;
for ( i=2 ; i<=n ; i++ ) e[i] = e[i-1] * e[1];

// Calculate coefficients

AB();

z[7] = C[0] + e[2]*C[2] + e[4] * C[4]; // ubar
z[8] = - e[2]*D[2] - e[4]*D[4];
z[9] = 0.5 * C[0]*C[0] + e[2]*E[2] + e[4] * E[4];

if(Current_criterion==1)
	{
	z[5] = Current*sqrt(kd);
	z[4] = z[7] + z[5];
	z[6] = z[4] + z[8]/kd - z[7];
	}

if(Current_criterion==2)
	{
	z[6] = Current*sqrt(kd);
	z[4] = z[6] - z[8]/kd + z[7];
	z[5] = z[4] - z[7];
	}

iff(Case,Wavelength) z[3] = 2*pi/z[4];
iff(Case,Period) z[3] = T * sqrt(kd);

for (i=1; i<=n; i++ )
	Tanh[i] = tanh(i*z[1]);

//	Output results and picture of wave

Solution=fopen("Solution.res","w");
Elevation = fopen("Surface.res","w");
Flowfield = fopen("Flowfield.res","w");

Output();

fflush(NULL);
printf("\nTouch key to continue "); getch();
printf("\n\nFinished\n");
}
