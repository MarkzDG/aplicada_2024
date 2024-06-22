#define	ANSI
#include	<math.h>
#include	<stdio.h>
#include	<conio.h>
#include	<stdlib.h>
#include	<string.h>

#define Readtext(stream,x)	fgets(x,400,stream); x[strlen(x)-1] = '\0'
#define iff(x,y)		if(strcmp(x,#y)==0)
#define For(x,y,z)		for((x)=(y);(x)<=(z);(x)++)
#define Skip(stream)	fgets(dummy,400,stream)
#define Read(stream,x,y) {fscanf(stream,"%"#y, &x);Skip(stream);}
#define pi				3.14159265358979324

#ifdef Main
#define Char char
#define Int int
#define Double double
#define Void void
#endif
#ifdef Subprograms
#define Char extern char
#define Int extern int
#define Double extern double
#define Void extern void
#endif

Char
	Heading[100], dummy[100], Known[20], Currentname[10];
Int
	Order, Current_criterion, Case;
Double
	alpha, c, ce, cs, Current, delta, e, epsilon, H, Highest, h, K, Kd, KK, L, m, m1, m1_limit, m4, q, q1, Q_d, R_d, T, U1, U2, Ubar_d;
Double
	ee[7], mm[7], Wavelength[7];
Double
	Aitken(double *, int), Alpha(double), cn(double), dn(double), E(double), eta_d, eta_h(double),
	hoverd(double H), lambda_d(double), lambda_series(double), Q_h(double), R_h(double), sn(double),
	u_h(double,double), Ubar_h(double), v_h(double,double);
Char
	Theory[10], Method[100];
Void
	Elliptic_integrals(void);
