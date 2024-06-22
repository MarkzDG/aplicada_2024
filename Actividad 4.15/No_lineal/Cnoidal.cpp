
// Cnoidal theory - C++ version

#define	Main
#include "Cnoidal.h"

// Main program

int main(void)
{
FILE		*monitor, *in, *in2, *Solution, *Flowfield, *Surface;
int		i, ii;
int		Points, Surface_points, Nprofiles;
int 		flushall(void);
char 		Current1[10]="Euler", Current2[10]="Stokes";
double	eta, x, y, u, v;
void		Solve(void), Input_Title_block(FILE*), Title_block(FILE*), Output(FILE*);

// This is the limiting value for m1:
m1_limit = 1.e-8;
// If m is c1oser to 1 than this, a different procedure is used,
// whereby m is set to 1 and K is calculated iteratively from the data

// Read data

in = fopen("../Data.dat","r");
monitor = stdout;
strcpy(Theory,"Cnoidal");

Readtext(in,Heading);
Read(in,H,lf);
fscanf(in,"%s", Known); Skip(in);

iff(Known,Wavelength)
	{
	Read(in,L,lf);
	if(L < 10.) {printf("\nThe dimensionless wavelength is less than 10.\nCnoidal theory should probably not be applied");}
	}

iff(Known,Period)
	{
	Read(in,T,lf);
	if(T < 10.) {printf("\nThe dimensionless period is less than 10.\nCnoidal theory should probably not be applied"); }
	}

Read(in,Case,d);
Read(in,Current,lf);
Read(in,Order,d);
if(Order>6) Order=6;
if(Order<6) sprintf(Method, "# Solution by %d-order cnoidal theory", Order);
else sprintf(Method, "# Solution by 5th order cnoidal theory with Aitken convergence enhancement");
fclose(in);

if (Case==1)
	{ce = Current;strcpy(Currentname,Current1);}
if (Case==2)
	{cs = Current;strcpy(Currentname,Current2);}

Input_Title_block(monitor);

// Solving for parameter m or K iteratively

printf("\n# Solving for parameter m or K iteratively\n\n      m          K\n");

// Initial estimate

iff(Known, Period)
   {
	m1 = 16.*exp(-sqrt(0.75*H*T*T));
	m = 1. - m1;
   }
iff(Known, Wavelength)
   {
	m1 = 16.*exp(-sqrt(0.75*H*L*L));
	m = 1.-m1;
   }

// Use direct iteration to solve for m, or K in the case of very long waves

Solve();

// Highest wave - eqn (32) of Fenton (1990)

Highest = (0.0077829*L*L*L+0.0095721*L*L+0.141063*L)/(0.0093407*L*L*L+0.0317567*L*L+0.078834*L+1);

// Evaluate all the global quantities

q = exp(-pi*Kd/K);

h = hoverd(H);
epsilon = H/h;
alpha = Alpha(epsilon);
delta = 4./3*alpha*alpha;
Ubar_d = Ubar_h(epsilon)*sqrt(h);
Q_d = Q_h(epsilon)*pow(h,1.5);
if (Case==1)
   {
	c = ce+Ubar_d;
   cs = c-Q_d;
   }
if (Case==2)
	{
   c = cs+Q_d;
   ce = c-Ubar_d;
	}

iff(Known,Wavelength) T = L/c;
iff(Known,Period) L = lambda_d(H);
R_d = R_h(epsilon)*h;
U1 = H*L*L;
U2 = U1/8/pi/pi;
if(U2 < 0.5) printf("\n\n# The program is being applied to a wave with Stokes-Ursell number < 1/2.\n# Results are unreliable");

// Output solution summary

Solution = fopen("Solution.res", "w");
fprintf(monitor,"\n\n# Solution summary:\n\n");
Title_block(monitor);
Title_block(Solution);
Output(Solution);
fclose(Solution);

// Output free surface

Surface = fopen("Surface.res","w");
fprintf(Surface,"#  %s\n",Heading);
in2 = fopen("../Points.dat","r");
Readtext(in2,dummy);
Read(in2,Surface_points,d);
Read(in2,Nprofiles,d);
Read(in2,Points,d);
fclose(in2);

fprintf(Surface, "\n%s",Method);
fprintf(Surface, "\n");
fprintf(Surface, "\n# Surface of wave - trough-crest-trough,");
fprintf(Surface, " note quadratic point spacing clustered around crest");
fprintf(Surface, "\n# Non-dimensionalised with respect to depth");
fprintf(Surface, "\n# X/d, eta/d, & check of surface pressure\n");

for (i=-Surface_points/2; i<=Surface_points/2; ++i)
	{
	x = (double)i/Surface_points;
	x = L*2*x*fabs(x);
	eta = eta_h(x/h);
	eta_d = eta*h;
	u = u_h(x/h, eta);
	v = v_h(x/h, eta);
	fprintf(Surface,"\n%8.4f\t%8.4f\t%8.0e", x, eta_d, 0.5*(u*u+v*v)+eta-R_d/h);
	}

fclose(Surface);

// Output flowfield summary

Flowfield = fopen("Flowfield.res", "w");
fprintf(Flowfield,"#  %s", Heading);
fprintf(Flowfield,"\n\n%s", Method);
fprintf(Flowfield,"\n\n# Velocity profiles\n");
fprintf(Flowfield,"\n# All quantities are dimensionless with respect to g and/or d\n");
fprintf(Flowfield,"\n#*********************");
fprintf(Flowfield,"\n# y        u       v");
fprintf(Flowfield,"\n# -     -------------");
fprintf(Flowfield,"\n# d        sqrt(gd)");
fprintf(Flowfield,"\n#*********************");

For(ii,0,Nprofiles)
	{
	x = L*0.5*(double)ii/Nprofiles;
	eta = eta_h(x/h);
	eta_d = eta*h;
	fprintf(Flowfield,"\n\n# x/d = %8.4f, Phase = %6.1f°\n", x, x/L*360);
	For(i,0,Points)
		{
		y = i*eta/Points;
		u = u_h(x/h, y)*sqrt(h)+c;
		v = v_h(x/h, y)*sqrt(h);
		fprintf(Flowfield,"\n%7.4f\t%7.4f\t%7.4f",y*h,u,v);
		}
	}

fclose(Flowfield);
fflush(NULL);
printf("\nTouch key to continue "); getch();
printf("\n\nFinished\n");
} // End main program

//	Two title blocks - at input and output

void Input_Title_block(FILE* file)
{
fprintf(file,"# %s", Heading);
fprintf(file,"\n\n# Printing input data here to check");
fprintf(file,"\n\n# Height/Depth:%6.3f", H);
iff(Known,Wavelength)
	{
	fprintf(file,"\n# Length/Depth:%7.2f", L);
	}
iff(Known,Period)
	{
	fprintf(file,"\n# Dimensionless Period T*sqrt(g/d):%7.2f", T);
	}
fprintf(file,"\n# Current criterion: %s,  Dimensionless value:%6.3lf\n", Currentname, Current);
fprintf(file,"\n%s\n", Method);
}

void Title_block(FILE* file)
{
fprintf(file,"# %s", Heading);
fprintf(file,"\n\n%s", Method);
fprintf(file,"\n\n# Height/Depth:%6.3f, %3.0lf\%% of the maximum of H/d =%6.3f for this length:",
	H,H/Highest*100., Highest);
fprintf(file,"\n# Length/Depth:%7.2f", L);
fprintf(file,"\n# Dimensionless Period T*sqrt(g/d):%7.2f", T);
fprintf(file,"\n# Current criterion: %s,  Dimensionless value:%6.3lf\n", Currentname, Current);
}
