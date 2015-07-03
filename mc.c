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
#include "test.h"

int R;
int mc_iter;

void update() //Basic HMC update step
{
 	int i, acc;
	double saveAt, saveAx, saveAy, prob, S0, S1;
	complex double saveUt, saveUx, saveUy, d;
  	double r[1];

	acc = 0;
	for(i = 0; i<GRIDPOINTS; i++) 
	{
		saveAt = At[i];
		saveUt = Ut[i];
		// propose a change
		S0 = localSt(i);

		ranlxd(r,1);
		At[i] = 2*M_PI*r[0];
		Ut[i] = cos(At[i]) + I*sin(At[i]);
		S1 = localSt(i);

		d = 1.0; //det_ratio_At(i);
		prob = cconj(d)*d*exp(S0-S1);
 	 	if(prob >= 1) {
    		R += 1;
			//update_inverse_At(i);
  		}
  		else {
    		ranlxd(r,1);
    		if(r[0] < prob) {
      			R += 1;
				//update_inverse_At(i);
    		}
    		else {
      		// reject the change, get the old values for A
			    At[i] = saveAt;
				Ut[i] = saveUt;
			}
		}

		// propose a change
		saveAx = Ax[i];
		saveUx = Ux[i];
		saveAy = Ay[i];
		saveUy = Uy[i];

		S0 = localSxy(i);

		ranlxd(r,1);
		Ax[i] = 2*M_PI*r[0];
		Ux[i] = cos(Ax[i]) + I*sin(Ax[i]);

		ranlxd(r,1);
		Ay[i] = 2*M_PI*r[0];
		Uy[i] = cos(Ay[i]) + I*sin(Ay[i]);
		
		S1 = localSxy(i);

		d = 1.0; //det_ratio_Axy(i);
		prob = cconj(d)*d*exp(S0-S1);
 	 	if(prob >= 1) {
    		R += 1;
			//update_inverse_Axy(i);
  		}
  		else {
    		ranlxd(r,1);
    		if(r[0] < prob) {
      			R += 1;
				//update_inverse_Axy(i);
    		}
    		else {
      		// reject the change, get the old values for A
			    Ax[i] = saveAx;
				Ux[i] = saveUx;
				Ay[i] = saveAy;
				Uy[i] = saveUy;
			}
		}
	}
 	 
 	//Increase the counter of the total number of MC updates...
 	mc_iter++;
 	//Return 1 if the configuration was accepted, 0 otherwise
}



