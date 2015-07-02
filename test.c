#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "complex/complex.h"
#include "linalg.h"
#include "fields.h"
#include "fermion.h"
#include "hmc.h"
#include "test.h"

void test_gauge_force(int n) {
	int i;
	double dA;
	double sg0, sg1;

	dA = 0.001;

	sg0 = 0.0;
 	for (i=0; i<GRIDPOINTS; i++)
 	{
  		sg0 += S_G(i);
 	};
	
	Ax[n] = Ax[n] + dA;
	sg1 = 0.0;
	for (i=0; i<GRIDPOINTS; i++)
 	{
  		sg1 += S_G(i);
 	};
	printf("%f\n", (sg1 - sg0)/dA);
	printf("%f\n", DS_Gx(n));

}

void fprint_fermion_mat() {
	int i, j;
	complex double x;
	complex double basis[GRIDPOINTS];
	complex double out[GRIDPOINTS];
	complex double temp[GRIDPOINTS];
	FILE *fp;
	
	fp = fopen("fmat_real.dat", "w");
	
	printf("\n Output fermion determinant...\n");
	set_zero(basis);
	for(i = 0; i<GRIDPOINTS; i++) 
	{
		basis[i] = 1.0;
		fermion_fp(out, temp, basis);
		//printf("{");
		for(j = 0; j < GRIDPOINTS-1; j++)
		{
			x = out[j];
			fprintf(fp, "%f  ", creal(x));
		}
		fprintf(fp, "%f\n", creal(out[GRIDPOINTS-1]));
		//printf("},");
		basis[i] = 0.0;
	}
	fclose(fp);

	fp = fopen("fmat_imag.dat", "w");

	for(i = 0; i<GRIDPOINTS; i++) 
	{
		basis[i] = 1.0;
		fermion_fp(out, temp, basis);
		//printf("{");
		for(j = 0; j < GRIDPOINTS-1; j++)
		{
			x = out[j];
			fprintf(fp, "%f  ", cimag(x));
		}
		fprintf(fp, "%f\n", cimag(out[GRIDPOINTS-1]));
		//printf("},");
		basis[i] = 0.0;
	}
	fclose(fp);
}

void test_fermion_force(int n) {
	int i, j;
	double squnrm;
	complex double x;
	double dA, f;

	for(i=0; i<GRIDPOINTS; i++)
 	{
  		g_R[i] = (gauss() + I*gauss())/sqrt(2); //Pseudofermion fields times M^{-1} 
 	};
	squnrm = square_norm(g_R);
  	fermion(g_fermion, g_R); //g_fermion the pseudofermion field, i.e. phi = M R
  	ham_old = squnrm;

	g_cgiterations1 += cg(g_eta, g_fermion, ITER_MAX, DELTACG, &fermion_sqr);
	f = fermion_forcex(n);

	dA = 0.0001;
	Ax[n] += dA;
	calculatelinkvars();
	g_cgiterations1 += cg(g_eta, g_fermion, ITER_MAX, DELTACG, &fermion_sqr);

	ham += scalar_prod_r(g_fermion, g_eta);

	printf("Fermion force: %f\n", f);
	printf("%f %f\n", f, (ham-ham_old)/dA);
	return;
}

double stupid_fermion_force_x(const int i) {
	double dA = 0.0001;
	double s0 = 0.0, s1 = 0.0;

	g_cgiterations1 += cg(g_temp2, g_fermion, ITER_MAX, DELTACG, &fermion_sqr);
	s0 = scalar_prod_r(g_fermion, g_temp2);

	Ax[i] = Ax[i] + dA;
	calculatelinkvars();

	g_cgiterations1 += cg(g_temp2, g_fermion, ITER_MAX, DELTACG, &fermion_sqr);
	s1 = scalar_prod_r(g_fermion, g_temp2);

	Ax[i] = Ax[i] - dA;
	calculatelinkvars();
	return (s1 - s0)/dA;
}

void calculate_fermion_force() 
{
  	int i;
	double fft[GRIDPOINTS];
	double ffx[GRIDPOINTS];
	double ffy[GRIDPOINTS];

  	g_cgiterations1 += cg(g_eta, g_fermion, ITER_MAX, DELTACG, &fermion_sqr);
  	for(i = 0; i < GRIDPOINTS; i++) {
		//printf("%d \n", i);
		//printf("Brute-force %f\n", stupid_fermion_force_x(i));
		//printf("Smart %f\n", fermion_forcex(i));
		ffx[i] = stupid_fermion_force_x(i);
		fft[i] = fermion_forcet(i);
		ffy[i] = fermion_forcey(i);
  	}
	print_vector_r(ffx);
	print_vector_r(fft);
	printf("Max x fermion force: %f\n", max_r(ffx));
	printf("Max t fermion force: %f\n", max_r(fft));
  	return;
}


