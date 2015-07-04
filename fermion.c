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
//complex double Minv3[GRIDPOINTS][GRIDPOINTS];

int up_counter;
int g_inverse = 6; //Number of quick inverse updates between hard-core inverse

double matrix_diff(fmat A1, fmat A2) 
{
	int i, j;
	double nm;
	nm = 0.0;
	for(i=0; i<GRIDPOINTS; i++) 
		for(j=0; j<GRIDPOINTS; j++)
			nm += cabs(A1[i][j] - A2[i][j])*cabs(A1[i][j] - A2[i][j]);
	return sqrt(nm);
}

void matrix_print(fmat A) 
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

complex double prod_row_col(const int i, const int j, fmat A) 
{
	return exp(-g_mu)*Ut[i]*A[tp[i]][j] - A[i][j] - g_t*(Ux[i]*A[xp[i]][j] + Uy[i]*A[yp[i]][j] + cconj(Ux[xm[i]])*A[xm[i]][j] + cconj(Uy[ym[i]])*A[ym[i]][j]);
}

complex double prod_col_row(const int i, const int j, fmat A) 
{
	return exp(-g_mu)*Ut[tm[i]]*A[j][tm[i]] - A[j][i] - g_t*(Ux[xm[i]]*A[j][xm[i]] + Uy[ym[i]]*A[j][ym[i]] + cconj(Ux[i])*A[j][xp[i]] + cconj(Uy[i])*A[j][yp[i]]);
}

complex double prod_row_vec(const int i, complex double *v)
{
	return exp(-g_mu)*Ut[i]*v[tp[i]] - v[i] - g_t*(Ux[i]*v[xp[i]] + Uy[i]*v[yp[i]] + cconj(Ux[xm[i]])*v[xm[i]] + cconj(Uy[ym[i]])*v[ym[i]]);
}

complex double det_ratio(const int i, fmat A) {
	int j, l, k;
	complex double rdet1, rdet2, rdet3, x, y;
	complex double colx[GRIDPOINTS], coly[GRIDPOINTS];
	
	// det M'/det M 
	rdet1 = prod_row_col(i, i, A);
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
	rdet2 = prod_row_vec(j, colx);
	
	y = prod_row_vec(xp[i], coly);
	for(k = 0; k < GRIDPOINTS; k++) 
		coly[k] = coly[k] - y/rdet2*colx[k];

	// calculate det M'''/det M''
	rdet3 = prod_row_vec(l, coly);

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
			for(k = 0; k < GRIDPOINTS; k++)
				// (new A)_j = A_j - (B^i A_j)/T A_i
				out[j][k] = in[j][k] - x/T*in[i][k];
		}
	}
}

// quick update the inverse after At, Ax, Ay at site i are modified
void quick_update_inverse(const int i, fmat in, fmat temp)
{
	update_row(i, temp, in);
	update_col(i, in, temp);
}


void hard_inverse(fmat M)
{
	int i, j;
	
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
			M[i][j] =  fdet[i][j];
	
	//r = matrix_det(*fdet);

	//printf("Det: %.12f + I* %.12f\n", creal(r), cimag(r));

	matrix_inverse(*M);
	matrix_print(M);
}

void update_inverse(int i, fmat M, fmat temp)
{
	if (up_counter < g_inverse) 
	{
		up_counter++;
		quick_update_inverse(i, M, temp);
	}
	else 
	{
		up_counter = 0; //reset the counter
		hard_inverse(M);
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


