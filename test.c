#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "complex/complex.h"
#include "linalg.h"
#include "fields.h"
#include "fermion.h"
#include "hmc.h"
#include "test.h"



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




