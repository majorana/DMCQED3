#ifndef _FIELDS_H
#define _FIELDS_H

#include "linalg.h"  
#include "lattice.h"
#include "rand/ranlxd.h"


/* declared in qed.c */
extern double g_mu;
extern double g_t;
extern double dt;
extern double beta0, beta;

extern double At[GRIDPOINTS], Ax[GRIDPOINTS], Ay[GRIDPOINTS];         //Non-compact real-valued gauge fields
extern complex double Ut[GRIDPOINTS], Ux[GRIDPOINTS], Uy[GRIDPOINTS];   //Compact lattice gauge fields: link = exp(I*gauge)


void coldstart(); //Cold start for the gauge fields
void hotstart();  //Hot (random) start for the gauge fields

void calculatelinkvars(); //Sets the values of the compact gauge fields from noncompact ones

double S_G(int i);   //Action of the gauge fields
double localS(int i); 

double S_Gtx(int i);
double S_Gxy(int i);

#endif
