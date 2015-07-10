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

void mc_init()
{
	emu = exp(-g_mu*dt);
	up_counter = 0;
	hard_inverse(Minv);
	range = 3.0;
}

void mc_update() //Basic MC update step
{
 	int i, acc;
	double saveAt, saveAx, saveAy, prob, S0, S1, phi;
	complex double saveUt, saveUx, saveUy, rdet;
  	double r[1], rA[1];

	acc = 0;
	for(i = 0; i<GRIDPOINTS; i++) 
	{
		// propose a change
		saveAt = At[i]; saveUt = Ut[i];

		S0 = localS(i);
		
		ranlxd(rA, 1);
		phi = range*(2.0*rA[0] - 1.0);
		Ut[i] *= (1.0+I*phi)/sqrt(1.0+phi*phi);
		At[i] += atan(phi);

	
		S1 = localS(i);

		rdet = det_ratio_t(i, Minv);
		//rdet = 1.0;
		prob = cconj(rdet)*rdet*exp(S0-S1);
		//printf("%.8f\n", cabs(rdet));
 	 	if(prob >= 1) {
    		R += 1;
			acc++;
			update_inverse(i, Minv);
			//matrix_print(Minv);
  		}
  		else {
    		ranlxd(r, 1);
    		if(r[0] < prob) {
      			R += 1;
				acc++;
				update_inverse(i, Minv);
				//matrix_print(Minv);
    		}
    		else {
      		// reject the change, get the old values for A
				At[i] = saveAt; Ut[i] = saveUt;
			    //Ax[i] = saveAx; Ux[i] = saveUx;
				//Ay[i] = saveAy; Uy[i] = saveUy;
			}
		}
	}
	
	for(i = 0; i<GRIDPOINTS; i++) 
	{
		// propose a change
		//saveAt = At[i]; saveUt = Ut[i];
		saveAx = Ax[i]; saveUx = Ux[i];

		S0 = localS(i);

		ranlxd(rA, 1);
		phi = range*(2.0*rA[0] - 1.0);
		Ux[i] *= (1.0+I*phi)/sqrt(1.0+phi*phi);
		Ax[i] += atan(phi);

	
		S1 = localS(i);

		rdet = det_ratio_xy(i, Minv);
		//rdet = 1.0;
		prob = cconj(rdet)*rdet*exp(S0-S1);
		//printf("%.8f\n", cabs(rdet));
 	 	if(prob >= 1) {
    		R += 1;
			acc++;
			update_inverse(i, Minv);
			//matrix_print(Minv);
  		}
  		else {
    		ranlxd(r, 1);
    		if(r[0] < prob) {
      			R += 1;
				acc++;
				update_inverse(i, Minv);
				//matrix_print(Minv);
    		}
    		else {
      		// reject the change, get the old values for A
			    Ax[i] = saveAx; Ux[i] = saveUx;
			}
		}
	}
	
	for(i = 0; i<GRIDPOINTS; i++) 
	{
		// propose a change
		//saveAt = At[i]; saveUt = Ut[i];
		saveAy = Ay[i]; saveUy = Uy[i];

		S0 = localS(i);
		
		ranlxd(rA, 1);
		phi = range*(2.0*rA[0] - 1.0);
		Uy[i] *= (1.0+I*phi)/sqrt(1.0+phi*phi);
		Ay[i] += atan(phi); 

	
		S1 = localS(i);

		rdet = det_ratio_xy(i, Minv);
		//rdet = 1.0;
		prob = cconj(rdet)*rdet*exp(S0-S1);
		//printf("%.8f\n", cabs(rdet));
 	 	if(prob >= 1) {
    		R += 1;
			acc++;
			update_inverse(i, Minv);
			//matrix_print(Minv);
  		}
  		else {
    		ranlxd(r, 1);
    		if(r[0] < prob) {
      			R += 1;
				acc++;
				update_inverse(i, Minv);
				//matrix_print(Minv);
    		}
    		else {
      		// reject the change, get the old values for A
				Ay[i] = saveAy; Uy[i] = saveUy;
			}
		}
	}
 	printf("Acceptance rate per update: %d, %d, %.5f\n", acc, 3*GRIDPOINTS, (double)acc/(3.0*GRIDPOINTS)); 
	fflush(stdout);
 	//Increase the counter of the total number of MC updates...
 	mc_iter++;
}



