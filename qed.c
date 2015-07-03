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

int g_thermalize   = 200;   //Number of MC updates for thermalization
int g_measurements = 400;    //Number of measurements (statistically independent configurations)
int g_intermediate =  2;    //Number of MC updates between the measurements

/* extern in fields.h   */

double beta0  = 1.0;
double beta   = 1.0;        //Coupling constant for the gauge field, allow anisotropy between space and time. This is a non-relativistic system.

void echo_sim_params();

int main(int argc, char **argv) 
{
	int i, j, l;
  	int accepted = 0;        //Total number of accepted configurations
  	int total_updates = 0;   //Total number of updates
	complex double det1, det2, detr;

  	/* Initialize the random number generator */
  	rlxd_init(2, time(NULL)); 
  	/* Initialize the lattice geometry */
  	init_lattice(Lx, Ly, Lt);
  	/* Initialize the fields */
  	hotstart();
  	/* Print out the run parameters */
  	echo_sim_params();
	
	det1 = get_fermion_mat();
	At[10] = 0.8;
	Ax[10] = 2.5;
	Ay[10] = 1.5;
	calculatelinkvars();
	detr = det_ratio(10);
	printf("%.12f+ I*%.12f\n", creal(detr), cimag(detr));
	det2 = get_fermion_mat();
	printf("%.12f+ I*%.12f\n", creal(det2/det1), cimag(det2/det1));

	return 0;

  	/* thermalization */
  	mc_iter = 0; //Counts the total number of calls to the update() routine
  	printf("\n Thermalization: \n\n");
  	for(i=0; i<g_thermalize; i++)
  	{
   		update();
		//printf("\t Step %04i,\t mp = %2.4lf,\t pl = %2.4lf,\t cc = %2.4lf\n", i, mean_plaquette(), polyakov_loop());
  	};
  	
	/* measure the iterations only during real simulation, not thermalization */
  	R              = 0; //Counts the total number of accepted configurations
  	mc_iter       = 0; //Counts the total number of calls to the update() routine

  	printf("\n Generation: \n\n");
  	for(i=0; i<g_measurements; i++) 
  	{
   	/* do g_intermediate updates before measurement */
   		for (l=0; l<g_intermediate; l++)
   		{
    		update();
   		};
   		update();
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
