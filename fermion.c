#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#ifndef M_PI
#define M_PI    3.14159265358979323846f
#endif
#include "fields.h"
#include "complex/complex.h"
#include "lattice.h"
#include "linalg.h"
#include "rand/ranlxd.h"
#include "fermion.h"

complex double fdet[GRIDPOINTS][GRIDPOINTS];
complex double mtemp[GRIDPOINTS][GRIDPOINTS];
complex double Minv1[GRIDPOINTS][GRIDPOINTS];
complex double Minv2[GRIDPOINTS][GRIDPOINTS];

// pointer to the 2D array
// should always point to the inverse matrix with the *current* configuration of gauge fields
// except when the matrices are updated
complex double (*Minv)[GRIDPOINTS];

// pointer to the spare 2D array
complex double (*Minv_spare)[GRIDPOINTS];

complex double det_ratio_At(const int i) {
// the i-th row
	complex double rdet;
	int j;
	//fdet[i][tp[i]] = exp(-g_mu)*Ut[i];

	rdet = 0.0;
	for(j=0;j<GRIDPOINTS;j++) {
		//printf("%.5f, %.5f\n", Minv[j][i]);
		//rdet = rdet + fdet[i][j]*Minv[j][i];
	}
	rdet = exp(-g_mu)*Ut[i]*Minv[tp[i]][i] - Minv[i][i] - g_t*(Ux[i]*Minv[xp[i]][i] + Uy[i]*Minv[yp[i]][i] + cconj(Ux[xm[i]])*Minv[xm[i]][i] + cconj(Uy[ym[i]])*Minv[ym[i]][i]);
	return rdet;
}

complex double det_ratio_Axy(const int i) {
	int j, k;
	complex double rdet1, rdet2, x;
	complex double col[GRIDPOINTS];
	j = xp[i];
	rdet1 = exp(-g_mu)*Ut[i]*Minv[tp[i]][i] - Minv[i][i] - g_t*(Ux[i]*Minv[xp[i]][i] + Uy[i]*Minv[yp[i]][i] + cconj(Ux[xm[i]])*Minv[xm[i]][i] + cconj(Uy[ym[i]])*Minv[ym[i]][i]);
	x = exp(-g_mu)*Ut[i]*Minv[tp[i]][j] - Minv[i][j] - g_t*(Ux[i]*Minv[xp[i]][j] + Uy[i]*Minv[yp[i]][j] + cconj(Ux[xm[i]])*Minv[xm[i]][j] + cconj(Uy[ym[i]])*Minv[ym[i]][j]);
	// need to find the j-th column of the updated inverse
	for(k = 0; k < GRIDPOINTS; k++) 
	{
		col[k] = Minv[k][j] - x/rdet1*Minv[k][i];
	}
	rdet2 = exp(-g_mu)*Ut[j]*col[tp[i]] - g_t*(Ux[i]*col[xp[i]] + Uy[i]*col[yp[i]] + cconj(Ux[xm[i]])*col[xm[i]] + cconj(Uy[ym[i]])*col[ym[i]]);
	return rdet1*rdet2;
}

void update_row(const int i) {
	int j, k;
	complex double x, T;
	complex double (*temp)[GRIDPOINTS];

	// B_i A^i
	T = exp(-g_mu)*Ut[i]*Minv[tp[i]][i] - Minv[i][i] - g_t*(Ux[i]*Minv[xp[i]][i] + Uy[i]*Minv[yp[i]][i] + cconj(Ux[xm[i]])*Minv[xm[i]][i] + cconj(Uy[ym[i]])*Minv[ym[i]][i]);


	for(j = 0; j < GRIDPOINTS; j++) 
	{
		if (j == i) 
		{
			for(k = 0; k < GRIDPOINTS; k++)
			{
				// (new A)^i = A^i/T
				Minv_spare[k][j] = 1.0/T*Minv[k][i];
			}
		}
		else
		{
			// x = B_i A^j
			x = exp(-g_mu)*Ut[i]*Minv[tp[i]][j] - Minv[i][j] - g_t*(Ux[i]*Minv[xp[i]][j] + Uy[i]*Minv[yp[i]][j] + cconj(Ux[xm[i]])*Minv[xm[i]][j] + cconj(Uy[ym[i]])*Minv[ym[i]][j]);
			for(k = 0; k < GRIDPOINTS; k++)
			{
				// (new A)^j = A^j - (B_i A^j)/T A^i
				Minv_spare[k][j] = Minv[k][j] - x/T*Minv[k][i];
			}
		}
	}

	// need to exchange Minv and Minv_spare
	temp = Minv;
	Minv = Minv_spare;
	Minv_spare = temp;
}

void update_col(const int i) {
	int j, k;
	complex double x, T;
	complex double (*temp)[GRIDPOINTS];

	// B^i A_i
	T = exp(-g_mu)*Ut[tm[i]]*Minv[i][tm[i]] - Minv[i][i] - g_t*(Ux[xm[i]]*Minv[i][xm[i]] + Uy[ym[i]]*Minv[i][ym[i]] + cconj(Ux[i])*Minv[i][xp[i]] + cconj(Uy[i])*Minv[i][yp[i]]);


	for(j = 0; j < GRIDPOINTS; j++) 
	{
		if (j == i) 
		{
			for(k = 0; k < GRIDPOINTS; k++)
			{
				// (new A)_i = A_i/T
				Minv_spare[j][k] = 1.0/T*Minv[i][k];
			}
		}
		else
		{
			// x = B^i A_j
			x = exp(-g_mu)*Ut[tm[i]]*Minv[j][tm[i]] - Minv[j][i] - g_t*(Ux[xm[i]]*Minv[j][xm[i]] + Uy[ym[i]]*Minv[j][ym[i]] + cconj(Ux[i])*Minv[j][xp[i]] + cconj(Uy[i])*Minv[j][yp[i]]);
			for(k = 0; k < GRIDPOINTS; k++)
			{
				// (new A)_j = A_j - (B^i A_j)/T A_i
				Minv_spare[j][k] = Minv[j][k] - x/T*Minv[i][k];
			}
		}
	}

	// need to exchange Minv and Minv_spare
	temp = Minv;
	Minv = Minv_spare;
	Minv_spare = temp;
}

void update_inverse_At(const int i) {
	update_row(i);
	return;
}

void update_inverse_Axy(const int i) {
	update_row(i);
	update_col(i);
	return;
}

void fermion(complex double *out, complex double *in) 
{
	int i;
  	for(i=0; i<GRIDPOINTS; i++) {
		out[i] = exp(-g_mu)*Ut[i]*in[tp[i]] - in[i] 
			- g_t*(Ux[i]*in[xp[i]] + cconj(Ux[xm[i]])*in[xm[i]] + Uy[i]*in[yp[i]] + cconj(Uy[ym[i]])*in[ym[i]]);
	}
	return;
}


void init_M()
{
	// M1 = matrix_inv
	Minv = Minv1;
	Minv_spare = Minv2;
}

/* ****************************************************************************************
 * Test routines 
 * ****************************************************************************************
 */
void print_fermion_mat() {
	int i, j;
	complex double x;
	complex double basis[GRIDPOINTS];
	complex double out[GRIDPOINTS];
	complex double temp[GRIDPOINTS];
	//FILE *fp;
	
	//fp = fopen("fmat_real.dat", "w");
	
	printf("\n Output fermion determinant...\n Real part:\n");
	set_zero(basis);
	for(i = 0; i<GRIDPOINTS; i++) 
	{
		basis[i] = 1.0;
		fermion(out, basis);
		//printf("{");
		for(j = 0; j < GRIDPOINTS-1; j++)
		{
			x = out[j];
			//fprintf(fp, "%f  ", creal(x));
			printf("%.3f  ", creal(x));

		}
		//fprintf(fp, "%f\n", creal(out[GRIDPOINTS-1]));
		printf("%f\n", creal(out[GRIDPOINTS-1]));
		//printf("},");
		basis[i] = 0.0;
	}
	//fclose(fp);

	//fp = fopen("fmat_imag.dat", "w");
	
	printf("Imaginary part:\n");
	for(i = 0; i<GRIDPOINTS; i++) 
	{
		basis[i] = 1.0;
		fermion(out, basis);
		//printf("{");
		for(j = 0; j < GRIDPOINTS-1; j++)
		{
			x = out[j];
			//fprintf(fp, "%f  ", cimag(x));
			printf("%.3f  ", cimag(x));
		}
		//fprintf(fp, "%f\n", cimag(out[GRIDPOINTS-1]));
		printf("%f\n", cimag(out[GRIDPOINTS-1]));

		//printf("},");
		basis[i] = 0.0;
	}
	//fclose(fp);
}


complex double get_fermion_mat()
{
	int i, j;
	complex double r;
	complex double basis[GRIDPOINTS];
	set_zero(basis);
	for(i = 0; i<GRIDPOINTS;i++)
	{
		set_zero(fdet[i]);
		fdet[i][tp[i]] = exp(-g_mu)*Ut[i];
		fdet[i][i] = -1.0;
		fdet[i][xp[i]] = -g_t;
		fdet[i][xm[i]] = -g_t;
		fdet[i][yp[i]] = -g_t;
		fdet[i][ym[i]] = -g_t;
	}

	for(i = 0; i<GRIDPOINTS;i++) {
		//printf("{");
		for(j = 0; j<GRIDPOINTS; j++) {
			Minv1[i][j] =  fdet[i][j];
			//printf("%.3f,  ", creal(fdet[i][j]));
		}
		//printf("%.3f", creal(fdet[i][j]));
		//printf("},\n");
	}
	printf("------------------\n");
	for(i = 0; i<GRIDPOINTS;i++) {
		//printf("{");
		for(j = 0; j<GRIDPOINTS-1; j++) {
			//printf("%.3f,  ", cimag(fdet[i][j]));
		}
		//printf("%.3f", cimag(fdet[i][j]));
		//printf("},\n");
	}
	r = matrix_det(*fdet);
	printf("Det: %.12f, %.12f\n", creal(r), cimag(r));
	matrix_inverse(*Minv1);
	for(i = 0; i<GRIDPOINTS;i++) {
		for(j = 0; j<GRIDPOINTS; j++) {
			//printf("%.3f  ", creal(Minv1[i][j]));
		}
		//printf("\n");
	}
	Minv = Minv1;
	return r;
}
