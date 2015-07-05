#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "complex/complex.h"
#include "rand/ranlxd.h"
#include "rand/gauss.h"
#include "linalg.h"
#include "fields.h"
#include "lattice.h"
#include "fermion.h"
#include "mc.h"

int R;
int mc_iter;
double range;

complex double Minv[GRIDPOINTS][GRIDPOINTS];
complex double temp[GRIDPOINTS][GRIDPOINTS];

void mc_init()
{
	up_counter = 0;
	hard_inverse(Minv);
	range = 1.7;
}

void propose(const int i)
{
	double rA[3];
	double phi;
	ranlxd(rA,3);
	//At[i] = 2*M_PI*rA[0];
	//Ut[i] = cos(At[i]) + I*sin(At[i]);
	//Ax[i] = 2*M_PI*rA[1];
	//Ux[i] = cos(Ax[i]) + I*sin(Ax[i]);
	//Ay[i] = 2*M_PI*rA[2];
	//Uy[i] = cos(Ay[i]) + I*sin(Ay[i]);
	phi = range*(2.0*rA[0] - 1.0);
	Ut[i] *= (1.0+I*phi)/sqrt(1.0+phi*phi);
	At[i] += atan(phi);
	phi = range*(2.0*rA[1] - 1.0);
	Ux[i] *= (1.0+I*phi)/sqrt(1.0+phi*phi);
	Ax[i] += atan(phi);
	phi = range*(2.0*rA[2] - 1.0);
	Uy[i] *= (1.0+I*phi)/sqrt(1.0+phi*phi);
	Ay[i] += atan(phi);
}

void mc_update() //Basic MC update step
{
 	int i, acc;
	double saveAt, saveAx, saveAy, prob, S0, S1;
	complex double saveUt, saveUx, saveUy, rdet;
  	double r[1];

	acc = 0;
	for(i = 0; i<GRIDPOINTS; i++) 
	{
		// propose a change
		saveAt = At[i]; saveUt = Ut[i];
		saveAx = Ax[i]; saveUx = Ux[i];
		saveAy = Ay[i]; saveUy = Uy[i];

		S0 = localS(i);

		propose(i);
	
		S1 = localS(i);

		rdet = det_ratio(i, Minv);
		//rdet = 1.0;
		prob = cconj(rdet)*rdet*exp(S0-S1);
		//printf("%.8f\n", prob);
 	 	if(prob >= 1) {
    		R += 1;
			update_inverse(i, Minv, temp);
			//matrix_print(Minv);
  		}
  		else {
    		ranlxd(r, 1);
    		if(r[0] < prob) {
      			R += 1;
				update_inverse(i, Minv, temp);
				//matrix_print(Minv);
    		}
    		else {
      		// reject the change, get the old values for A
				At[i] = saveAt; Ut[i] = saveUt;
			    Ax[i] = saveAx; Ux[i] = saveUx;
				Ay[i] = saveAy; Uy[i] = saveUy;
			}
		}
	}
 	 
 	//Increase the counter of the total number of MC updates...
 	mc_iter++;
}



