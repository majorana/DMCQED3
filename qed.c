#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#include "lattice.h"
#include "mc.h"
#include "fields.h"
#include "measurement.h"

#include "fermion.h"
#include "linalg.h"

/* ***************************************************************************************************************** */
// global variables 

double g_mu = 1.0;
double g_t = 1.0;
double dt = 8.0/(double)Lt;
double beta0  = 1.0;
double beta   = 2.0;        //Coupling constant for the gauge field, allow anisotropy between space and time. This is a non-relativistic system.

int g_thermalize   = 0;   //Number of MC updates for thermalization; a few hundreds
int g_measurements = 100;    //Number of measurements (statistically independent configurations)
int g_intermediate =  1;    //Number of MC updates between the measurements

/* ***************************************************************************************************************** */

int measure_iter = 0;

void echo_sim_params();

void output_measurement();

void test();

int main(int argc, char **argv) 
{
	int i, l;
  	int accepted = 0;        //Total number of accepted configurations
  	int total_updates = 0;   //Total number of updates

  	/* Initialize the random number generator */
  	rlxd_init(2, time(NULL)); 
  	/* Initialize the lattice geometry */
  	init_lattice(Lx, Ly, Lt);
  	/* Initialize the fields */
  	hotstart();
  	/* Print out the run parameters */
  	echo_sim_params();
	
	//test();
	//return 0;

	mc_init();
  	/* thermalization */
  	mc_iter = 0; //Counts the total number of calls to the update() routine
  	printf("\n Thermalization: \n\n");
  	for(i=0; i<g_thermalize; i++)
  	{
   		mc_update();
		//printf("\t Step %04i\n", i);
  	};
	/* measure the iterations only during real simulation, not thermalization */
  	R              = 0; //Counts the total number of accepted configurations
  	mc_iter       = 0; //Counts the total number of calls to the update() routine

  	printf("\n Generation: \n\n");
	measurement_init();
 	density(Minv);
	printf("Average density: \t %.5f %.5f\n", creal(m_density[measure_iter]), cimag(m_density[measure_iter]));
	printf("Wilson plaquette: \t %.5f\n", mean_plaq());

  	for(i=0; i<g_measurements; i++) 
  	{
   	/* do g_intermediate updates before measurement */
   		for (l=0; l<g_intermediate; l++)
   		{
    		mc_update();
   		};
   		mc_update();

		/* doing measurement */
  		density(Minv);
		//density_correlation(Minv);
		//wilson_loop_xy();
		printf("Average density: \t %.5f %.5f\n", creal(m_density[measure_iter]), cimag(m_density[measure_iter]));
		printf("Wilson plaquette: \t %.5f\n", mean_plaq());
		fflush(stdout);
		measure_iter++;
  	};
	output_measurement();
 	measurement_finish();
	
  	/* Some output for diagnostics */
  	total_updates = g_measurements*(g_intermediate + 1)*GRIDPOINTS;
  	printf("\n\n Algorithm performance:\n");
  	printf("\t Acceptance rate:             %.4f\n", (double)R/(double)total_updates);

  	return 0;
}

void echo_sim_params()
{
 	printf("\n Quantum Monte-Carlo for U(1) gauge theory with spinor Fermi surface\n\n");
 	printf("Run parameters:\n");
	printf("\t Chemical potential:              %2.4lf\n",   g_mu);
 	printf("\t Beta0:                            %2.4lf\n",  beta0);
	printf("\t Beta:                            %2.4lf\n",  beta);
	printf("\t Temperature:                     %2.4lf\n", dt*Lt);
 	printf("\t Lattice size:                    %i x %i x %i\n", Lx, Ly, Lt);
 	printf("\t Thermalization steps:            %i\n",      g_thermalize);
 	printf("\t Number of measurements:          %i\n",      g_measurements);
 	printf("\t MC updates between measurements: %i\n",      g_intermediate);
	printf("\n\n");
	fflush(stdout);
}

void output_measurement()
{
	int i, j, k;
	FILE *fp;

	fp = fopen("density.dat", "w");
	for(i = 0; i < g_measurements; i++)
		fprintf(fp, "%.5f %.5f\n", creal(m_density[i]), cimag(m_density[i]));
	fclose(fp);

	fp = fopen("wilsonxy.dat", "w");
	for(i = 0; i < g_measurements; i++)
	{
		for(j = 0; j < Lx/2; j++)
			for(k = 0; k < Ly/2; k++)
				fprintf(fp, "\t %.5f %.5f", creal(m_wilson_xy[i][j][k]), cimag(m_wilson_xy[i][j][k]));
		fprintf(fp, "\n");
	}
	fclose(fp);
	
	fp = fopen("density_corr.dat", "w");
	for(i = 0; i < g_measurements; i++)
	{
		for(j = 0; j < Lx; j++)
			for(k = 0; k < Ly; k++)
				fprintf(fp, "\t %.5f %.5f", creal(m_density_corr[i][j][k]), cimag(m_density_corr[i][j][k]) );
		fprintf(fp, "\n");
	}
	fclose(fp);
}

void test()
{
	int i;
	complex double detr, det1, det2;
	double rA[2];


	//det2 = hard_inverse(Minv1);
	for(i = 0; i < 60; i++) {
		//At[i] += 1.8;
		ranlxd(rA, 2);
		At[i] += rA[0];
		//Ax[i] += rA[0];
		//Ay[i] += rA[1];
		calculatelinkvars();
		det1 = det2;
		//det2 = hard_inverse(Minv2);
		detr = det_ratio_t(i, Minv1);
		update_inverse(i, Minv1); 
		printf("%.12f+ I*%.12f, %.12f\n", creal(detr), cimag(detr), cabs(detr));
		printf("%.12f+ I*%.12f, %.12f\n\n", creal(det2/det1), cimag(det2/det1), cabs(det2/det1));
		//printf("%g\n", matrix_diff(Minv1, Minv2));
	}

}
