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

double matrix_diff(complex double (*A1)[GRIDPOINTS],  complex double (*A2)[GRIDPOINTS]) 
{
	int i, j;
	double nm;
	nm = 0.0;
	for(i=0; i<GRIDPOINTS; i++) 
		for(j=0; j<GRIDPOINTS; j++)
			nm += cabs(A1[i][j] - A2[i][j])*cabs(A1[i][j] - A2[i][j]);
	return sqrt(nm);
}

void matrix_print(complex double (*A)[GRIDPOINTS]) 
{
	int i, j;
	for(i = 0; i<GRIDPOINTS;i++) {
		for(j = 0; j<GRIDPOINTS; j++) {
			printf("%.2f ", creal(A[i][j]));
		}
		printf("\n");
	}
	printf("---------------------------------------------\n");
	for(i = 0; i<GRIDPOINTS;i++) {
		for(j = 0; j<GRIDPOINTS; j++) {
			printf("%.2f ", cimag(A[i][j]));
		}
		printf("\n");
	}
	return;
}

complex double prod_row_col(int i, int j) 
{
	return exp(-g_mu)*Ut[i]*Minv[tp[i]][j] - Minv[i][j] - g_t*(Ux[i]*Minv[xp[i]][j] + Uy[i]*Minv[yp[i]][j] + cconj(Ux[xm[i]])*Minv[xm[i]][j] + cconj(Uy[ym[i]])*Minv[ym[i]][j]);
}

complex double prod_col_row(int i, int j) 
{
	return exp(-g_mu)*Ut[tm[i]]*Minv[j][tm[i]] - Minv[j][i] - g_t*(Ux[xm[i]]*Minv[j][xm[i]] + Uy[ym[i]]*Minv[j][ym[i]] + cconj(Ux[i])*Minv[j][xp[i]] + cconj(Uy[i])*Minv[j][yp[i]]);
}

complex double prod_row_vec(int i, complex double *v)
{
	return exp(-g_mu)*Ut[i]*v[tp[i]] - v[i] - g_t*(Ux[i]*v[xp[i]] + Uy[i]*v[yp[i]] + cconj(Ux[xm[i]])*v[xm[i]] + cconj(Uy[ym[i]])*v[ym[i]]);
}

complex double det_ratio(const int i) {
	int j, l, k;
	complex double rdet1, rdet2, rdet3, r, x, y;
	complex double colx[GRIDPOINTS], coly[GRIDPOINTS];
	//rdet1 = exp(-g_mu)*Ut[i]*Minv[tp[i]][i] - Minv[i][i] - g_t*(Ux[i]*Minv[xp[i]][i] + Uy[i]*Minv[yp[i]][i] + cconj(Ux[xm[i]])*Minv[xm[i]][i] + cconj(Uy[ym[i]])*Minv[ym[i]][i]);
	
	// det M'/det M 
	rdet1 = prod_row_col(i, i);
	//x = exp(-g_mu)*Ut[i]*Minv[tp[i]][j] - Minv[i][j] - g_t*(Ux[i]*Minv[xp[i]][j] + Uy[i]*Minv[yp[i]][j] + cconj(Ux[xm[i]])*Minv[xm[i]][j] + cconj(Uy[ym[i]])*Minv[ym[i]][j]);

	// calculate the xp[i] and yp[i] rows of the inv of M'
	j = xp[i];
	x = prod_row_col(i, j);
	// need to find the j-th column of the inverse of M' 
	for(k = 0; k < GRIDPOINTS; k++) 
		colx[k] = Minv[k][j] - x/rdet1*Minv[k][i];
	

	l = yp[i];
	y = prod_row_col(i, l);
	for(k = 0; k < GRIDPOINTS; k++) 
		coly[k] = Minv[k][l] - y/rdet1*Minv[k][i];

	
	// calculate det M''/ det M'
	//rdet2 = exp(-g_mu)*Ut[j]*col[tp[j]] - col[j] - g_t*(Ux[j]*col[xp[j]] + Uy[j]*col[yp[j]] + cconj(Ux[xm[j]])*col[xm[j]] + cconj(Uy[ym[j]])*col[ym[j]]);
	rdet2 = prod_row_vec(j, colx);
	
	y = prod_row_vec(xp[i], coly);
	for(k = 0; k < GRIDPOINTS; k++) 
		coly[k] = coly[k] - y/rdet2*colx[k];

	// calculate det M'''/det M''
	rdet3 = prod_row_vec(l, coly);
	r = rdet1*rdet2*rdet3;
	return rdet1*rdet2*rdet3;
}

void update_row(const int i, complex double rdet) {
	int j, k;
	complex double x;
	complex double (*temp)[GRIDPOINTS];
	
	// just for test
	Minv_spare = Minv2;

	for(j = 0; j < GRIDPOINTS; j++) 
	{
		if (j == i) 
		{
			for(k = 0; k < GRIDPOINTS; k++)
			{
				// (new A)^i = A^i/T
				Minv_spare[k][i] = 1.0/rdet*Minv[k][i];
			}
		}
		else
		{
			// x = B_i A^j
			//x = exp(-g_mu)*Ut[i]*Minv[tp[i]][j] - Minv[i][j] - g_t*(Ux[i]*Minv[xp[i]][j] + Uy[i]*Minv[yp[i]][j] + cconj(Ux[xm[i]])*Minv[xm[i]][j] + cconj(Uy[ym[i]])*Minv[ym[i]][j]);
			x = prod_row_col(i, j);
			for(k = 0; k < GRIDPOINTS; k++)
			{
				// (new A)^j = A^j - (B_i A^j)/T A^i
				Minv_spare[k][j] = Minv[k][j] - x/rdet*Minv[k][i];
			}
		}
	}

	// need to exchange Minv and Minv_spare
	//temp = Minv;
	//Minv = Minv_spare;
	//Minv_spare = temp;
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
	//temp = Minv;
	//Minv = Minv_spare;
	//Minv_spare = temp;
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
		fdet[i][xp[i]] = -g_t*Ux[i];
		fdet[i][xm[i]] = -g_t*cconj(Ux[xm[i]]);
		fdet[i][yp[i]] = -g_t*Uy[i];
		fdet[i][ym[i]] = -g_t*cconj(Uy[ym[i]]);
	}

	for(i = 0; i<GRIDPOINTS;i++) {
		//printf("{");
		for(j = 0; j<GRIDPOINTS; j++) {
			Minv1[i][j] =  fdet[i][j];
		}
	}
	
	r = matrix_det(*fdet);

	printf("Det: %.12f, %.12f\n", creal(r), cimag(r));

	matrix_inverse(*Minv1);
	//matrix_print(Minv1);
	
	Minv = Minv1;
	return r;
}
