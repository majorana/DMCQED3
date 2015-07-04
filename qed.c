#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#include "lattice.h"
#include "mc.h"
#include "fields.h"

#include "fermion.h"
#include "linalg.h"

/* global variables */
double g_mu = 1.0;
double g_t = 1.0;

int g_thermalize   = 10;   //Number of MC updates for thermalization
int g_measurements = 400;    //Number of measurements (statistically independent configurations)
int g_intermediate =  2;    //Number of MC updates between the measurements

/* extern in fields.h   */

double beta0  = 1.0;
double beta   = 1.0;        //Coupling constant for the gauge field, allow anisotropy between space and time. This is a non-relativistic system.

void echo_sim_params();

int main(int argc, char **argv) 
{
	int i, l;
  	int accepted = 0;        //Total number of accepted configurations
  	int total_updates = 0;   //Total number of updates
	complex double detr;

  	/* Initialize the random number generator */
  	rlxd_init(2, time(NULL)); 
  	/* Initialize the lattice geometry */
  	init_lattice(Lx, Ly, Lt);
  	/* Initialize the fields */
  	hotstart();
  	/* Print out the run parameters */
  	echo_sim_params();
	
	hard_inverse(Minv1);
	i = 20;
	// inverse is stored in Minv
	At[i] = 3.8;
	Ax[i] = 2.5;
	Ay[i] = 0.0;
	calculatelinkvars();
	detr = det_ratio(i, Minv1);
	quick_update_inverse(i, Minv1, Minv2); // Minv2 is temporary storage
	printf("%.12f+ I*%.12f\n", creal(detr), cimag(detr));
	hard_inverse(Minv2);
	//printf("%.12f+ I*%.12f\n", creal(det2/det1), cimag(det2/det1));
	printf("%.12f\n", matrix_diff(Minv1, Minv2));
	return 0;

	mc_init();
  	/* thermalization */
  	mc_iter = 0; //Counts the total number of calls to the update() routine
  	printf("\n Thermalization: \n\n");
  	for(i=0; i<g_thermalize; i++)
  	{
   		mc_update();
		printf("\t Step %04i\n", i);
  	};
	return 0;

	/* measure the iterations only during real simulation, not thermalization */
  	R              = 0; //Counts the total number of accepted configurations
  	mc_iter       = 0; //Counts the total number of calls to the update() routine

  	printf("\n Generation: \n\n");
  	for(i=0; i<g_measurements; i++) 
  	{
   	/* do g_intermediate updates before measurement */
   		for (l=0; l<g_intermediate; l++)
   		{
    		mc_update();
   		};
   		mc_update();
   	/* Measurements should go here... */
  
  	};
 
  	/* Some output for diagnostics */
  	total_updates = g_measurements*(g_intermediate + 1)*GRIDPOINTS;
  	printf("\n\n Algorithm performance:\n");
  	printf("\t Acceptance rate:             %.3f\n", (double)R/(double)total_updates);

  	return 0;
}

void echo_sim_params()
{
 	printf("Quantum Monte-Carlo for U(1) gauge theory with spinor Fermi surface\n\n");
 	printf("Run parameters (adjust in qed.c !!!):\n");
 	printf("\t Beta0:                            %2.4lf\n",  beta0);
	printf("\t Beta:                            %2.4lf\n",  beta);
 	printf("\t Lattice size:                    %i x %i x %i\n", Lx, Ly, Lt);
 	printf("\t Thermalization steps:            %i\n",      g_thermalize);
 	printf("\t Number of measurements:          %i\n",      g_measurements);
 	printf("\t MC updates between measurements: %i\n",      g_intermediate);
	printf("\n\n");
}
