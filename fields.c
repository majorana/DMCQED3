#include <stdlib.h>
#include <math.h>
#ifndef M_PI
#define M_PI    3.14159265358979323846f
#endif
#include <complex.h>
#include "lattice.h"
#include "linalg.h"
#include "rand/ranlxd.h"
#include "fields.h"


double At[GRIDPOINTS], Ax[GRIDPOINTS], Ay[GRIDPOINTS];         //Non-compact real-valued gauge fields
complex double Ut[GRIDPOINTS], Ux[GRIDPOINTS], Uy[GRIDPOINTS];   //Compact lattice gauge fields: link = exp(I*gauge)

double S_Gtx(int i)
{
	return -beta0*cos(At[i] + Ax[tp[i]] - At[xp[i]] - Ax[i]);
}

double S_Gty(int i)
{
	return -beta0*cos(At[i] + Ay[tp[i]] - At[yp[i]] - Ay[i]);
}

double S_Gxy(int i)
{
	return -beta*cos(Ax[i] + Ay[xp[i]] - Ax[yp[i]] - Ay[i]);
}


double S_G(int i)
{
	//return (-beta0*(cos(At[i] + Ax[tp[i]] - At[xp[i]] - Ax[i]) 
	//		+ cos(At[i] + Ay[tp[i]] - At[yp[i]] - Ay[i])) 
	//		- beta*cos(Ax[i] + Ay[xp[i]] - Ax[yp[i]] - Ay[i]) );
	return S_Gtx(i) + S_Gty(i) + S_Gxy(i);
}

double localSt(int i)
{
	return S_Gtx(i) + S_Gtx(xm[i]) + S_Gty(i) + S_Gty(ym[i]);
}

double localSxy(int i)
{
	return S_Gtx(i) + S_Gtx(tm[i]) + S_Gty(tm[i]) + S_Gxy(i) + S_Gxy(xm[i]) + S_Gxy(ym[i]);
}

void coldstart()
{
 	int i;
 	for(i=0; i<GRIDPOINTS; i++)
 	{
  		At[i]     = 0.0;
  		Ax[i]     = 0.0;
		Ay[i]     = 0.0;

 	};
 	calculatelinkvars();
}

void hotstart()
{
 	int i;
 	double r[GRIDPOINTS*3];
 	ranlxd(r, GRIDPOINTS*3);
 	for(i=0; i<GRIDPOINTS; i++)
 	{
  		At[i]=2*M_PI*r[i]-M_PI;
  		Ax[i]=2*M_PI*r[i+GRIDPOINTS]-M_PI;
		Ay[i]=2*M_PI*r[i+2*GRIDPOINTS]-M_PI;
 	}
 	calculatelinkvars();
}

void calculatelinkvars()
{
 	int i;
 	for(i=0; i<GRIDPOINTS; i++)
 	{
  		Ut[i] = cos(At[i]) + I*sin(At[i]);
  		Ux[i] = cos(Ax[i]) + I*sin(Ax[i]);
		Uy[i] = cos(Ay[i]) + I*sin(Ay[i]);
 	};
}
