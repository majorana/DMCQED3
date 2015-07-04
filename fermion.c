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
complex double Minv1[GRIDPOINTS][GRIDPOINTS];
complex double Minv2[GRIDPOINTS][GRIDPOINTS];
complex double Minv3[GRIDPOINTS][GRIDPOINTS];

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

complex double prod_row_col(int i, int j, fmat A) 
{
	return exp(-g_mu)*Ut[i]*A[tp[i]][j] - A[i][j] - g_t*(Ux[i]*A[xp[i]][j] + Uy[i]*A[yp[i]][j] + cconj(Ux[xm[i]])*A[xm[i]][j] + cconj(Uy[ym[i]])*A[ym[i]][j]);
}

complex double prod_col_row(int i, int j, fmat A) 
{
	return exp(-g_mu)*Ut[tm[i]]*A[j][tm[i]] - A[j][i] - g_t*(Ux[xm[i]]*A[j][xm[i]] + Uy[ym[i]]*A[j][ym[i]] + cconj(Ux[i])*A[j][xp[i]] + cconj(Uy[i])*A[j][yp[i]]);
}

complex double prod_row_vec(int i, complex double *v)
{
	return exp(-g_mu)*Ut[i]*v[tp[i]] - v[i] - g_t*(Ux[i]*v[xp[i]] + Uy[i]*v[yp[i]] + cconj(Ux[xm[i]])*v[xm[i]] + cconj(Uy[ym[i]])*v[ym[i]]);
}

complex double det_ratio(const int i, fmat A) {
	int j, l, k;
	complex double rdet1, rdet2, rdet3, r, x, y;
	complex double colx[GRIDPOINTS], coly[GRIDPOINTS];
	//rdet1 = exp(-g_mu)*Ut[i]*Minv[tp[i]][i] - Minv[i][i] - g_t*(Ux[i]*Minv[xp[i]][i] + Uy[i]*Minv[yp[i]][i] + cconj(Ux[xm[i]])*Minv[xm[i]][i] + cconj(Uy[ym[i]])*Minv[ym[i]][i]);
	
	// det M'/det M 
	rdet1 = prod_row_col(i, i, A);
	//x = exp(-g_mu)*Ut[i]*Minv[tp[i]][j] - Minv[i][j] - g_t*(Ux[i]*Minv[xp[i]][j] + Uy[i]*Minv[yp[i]][j] + cconj(Ux[xm[i]])*Minv[xm[i]][j] + cconj(Uy[ym[i]])*Minv[ym[i]][j]);
	// calculate the xp[i] and yp[i] rows of the inv of M'
	j = xp[i];
	x = prod_row_col(i, j, A);
	// need to find the j-th column of the inverse of M' 
	for(k = 0; k < GRIDPOINTS; k++) 
		colx[k] = A[k][j] - x/rdet1*A[k][i];
	

	l = yp[i];
	y = prod_row_col(i, l, A);
	for(k = 0; k < GRIDPOINTS; k++) 
		coly[k] = A[k][l] - y/rdet1*A[k][i];

	
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

void update_row(const int i, fmat out, fmat in) {
	int j, k;
	complex double x, rdet;
	
	rdet = prod_row_col(i, i, in);
	for(j = 0; j < GRIDPOINTS; j++) 
	{
		if (j == i) 
		{
			// (new A)^i = A^i/T
			for(k = 0; k < GRIDPOINTS; k++)
				out[k][i] = 1.0/rdet*in[k][i];
		}
		else
		{
			// x = B_i A^j
			x = prod_row_col(i, j, in);
			// (new A)^j = A^j - (B_i A^j)/T A^i
			for(k = 0; k < GRIDPOINTS; k++)
				out[k][j] = in[k][j] - x/rdet*in[k][i];
		}
	}
}

void update_col(const int i, fmat out, fmat in) {
	int j, k;
	complex double x, T;

	// B^i A_i
	T = prod_col_row(i, i, in);

	for(j = 0; j < GRIDPOINTS; j++) 
	{
		if (j == i) 
		{
			// (new A)_i = A_i/T
			for(k = 0; k < GRIDPOINTS; k++)
				out[j][k] = 1.0/T*in[i][k];
		}
		else
		{
			// x = B^i A_j
			x = prod_col_row(i, j, in);
				//exp(-g_mu)*Ut[tm[i]]*Minv[j][tm[i]] - Minv[j][i] - g_t*(Ux[xm[i]]*Minv[j][xm[i]] + Uy[ym[i]]*Minv[j][ym[i]] + cconj(Ux[i])*Minv[j][xp[i]] + cconj(Uy[i])*Minv[j][yp[i]]);
			for(k = 0; k < GRIDPOINTS; k++)
			{
				// (new A)_j = A_j - (B^i A_j)/T A_i
				out[j][k] = in[j][k] - x/T*in[i][k];
			}
		}
	}
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



/* ****************************************************************************************
 * Test routines 
 * ****************************************************************************************
 */
void print_fermion_mat() {
	int i, j;
	complex double x;
	complex double basis[GRIDPOINTS];
	complex double out[GRIDPOINTS];
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


complex double get_fermion_mat(fmat M)
{
	int i, j;
	complex double r;
	
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

	for(i = 0; i<GRIDPOINTS;i++) 
		for(j = 0; j<GRIDPOINTS; j++) 
			Minv1[i][j] =  fdet[i][j];
	
	r = matrix_det(*fdet);

	printf("Det: %.12f + I* %.12f\n", creal(r), cimag(r));

	matrix_inverse(*M);
	
	return r;
}
