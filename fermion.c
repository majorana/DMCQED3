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

complex double temp1[GRIDPOINTS][GRIDPOINTS];
complex double temp2[GRIDPOINTS][GRIDPOINTS];


int up_counter;
int g_inverse = REINVERSE; //Number of quick inverse updates between hard-core inverse

double emu;

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

complex double apb(const int i)
// matrix element <i|M|i+dtau>
{
	int it, ix, iy;
	coordiate(i, &it, &ix, &iy);
	if (it == Lt - 1)
		return -1.0;
	else
		return 1.0;
}

complex double prod_row_col(const int i, const int j, fmat A) 
{
	//int k;
	//complex double r = 0.0;
	return -emu*apb(i)*Ut[i]*A[tp[i]][j] + (1.0)*A[i][j] + dt*Ut[i]*apb(i)*(Ux[tp[i]]*A[tp[xp[i]]][j] + Uy[tp[i]]*A[tp[yp[i]]][j] + cconj(Ux[tp[xm[i]]])*A[tp[xm[i]]][j] + cconj(Uy[tp[ym[i]]])*A[tp[ym[i]]][j]);
	//for(k = 0;k<GRIDPOINTS;k++)
	//	r += fdet[i][k]*A[k][j];
	//return r;
}

complex double prod_col_row(const int i, const int j, fmat A) 
{
	//int k;
	//complex double r = 0.0;
	//for(k = 0;k<GRIDPOINTS;k++)
	//	r += fdet[k][i]*A[j][k];
	//return r;
	return -emu*Ut[tm[i]]*apb(tm[i])*A[j][tm[i]] + (1.0)*A[j][i] + dt*apb(tm[i])*(Ut[tm[xm[i]]]*Ux[xm[i]]*A[j][tm[xm[i]]] + Ut[tm[ym[i]]]*Uy[ym[i]]*A[j][tm[ym[i]]] + Ut[tm[xp[i]]]*cconj(Ux[i])*A[j][tm[xp[i]]] + Ut[tm[yp[i]]]*cconj(Uy[i])*A[j][tm[yp[i]]]);
}

complex double prod_row_vec(const int i, complex double *v)
{
	return -emu*apb(i)*Ut[i]*v[tp[i]] + (1.0)*v[i] + dt*Ut[i]*apb(i)*(Ux[tp[i]]*v[tp[xp[i]]] + Uy[tp[i]]*v[tp[yp[i]]] + cconj(Ux[tp[xm[i]]])*v[tp[xm[i]]] + cconj(Uy[tp[ym[i]]])*v[tp[ym[i]]]);
}

complex double det_ratio_t(const int i, fmat A) {
	// det M'/det M, update row i
	return prod_row_col(i, i, A);
}

complex double det_ratio_xy(const int i, fmat A) {
	int j, l, k;
	complex double rdet1, rdet2, rdet3, x, y;
	complex double colx[GRIDPOINTS], coly[GRIDPOINTS];
	
	// det M'/det M, update row tm[i]
	rdet1 = prod_row_col(tm[i], tm[i], A);
	// calculate the j=tm[xp[i]] row of the inv of M'
	j = tm[xp[i]];
	x = prod_row_col(tm[i], j, A);
	// need to find the j-th column of the inverse of M' 
	for(k = 0; k < GRIDPOINTS; k++) 
		colx[k] = A[k][j] - x/rdet1*A[k][tm[i]];
	
	l = tm[yp[i]];
	y = prod_row_col(tm[i], l, A);
	for(k = 0; k < GRIDPOINTS; k++) 
		coly[k] = A[k][l] - y/rdet1*A[k][tm[i]];


	// calculate det M''/ det M'
	rdet2 = prod_row_vec(j, colx);
	
	y = prod_row_vec(j, coly);
	for(k = 0; k < GRIDPOINTS; k++) 
		coly[k] = coly[k] - y/rdet2*colx[k];

	// calculate det M'''/det M''
	rdet3 = prod_row_vec(l, coly);

	return rdet1*rdet2*rdet3;
}



complex double det_ratio_old(const int i, fmat A) {
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

// quick update the inverse after Ax, Ay at site i are modified
void quick_update_inverse(const int i, fmat in, fmat temp1, fmat temp2)
{
	update_row(i, temp1, in);
	update_row(tm[i], temp2, temp1);
	update_col(i, in, temp2);
}


void hard_inverse(fmat M)
{
	int i, j, info;
	complex double r;
	complex double s;
	
	for(i = 0; i<GRIDPOINTS;i++)
	{
		set_zero(M[i]);
		s = apb(i);
		M[i][tp[i]] = -emu*s*Ut[i];
		M[i][i] = 1.0;
		M[i][tp[xp[i]]] = dt*Ut[i]*s*Ux[tp[i]];
		M[i][tp[xm[i]]] = dt*Ut[i]*s*cconj(Ux[tp[xm[i]]]);
		M[i][tp[yp[i]]] = dt*Ut[i]*s*Uy[tp[i]];
		M[i][tp[ym[i]]] = dt*Ut[i]*s*cconj(Uy[tp[ym[i]]]);
	}

	//for(i = 0;i<GRIDPOINTS;i++)
	//	for(j = 0;j < GRIDPOINTS;j++)
	//		fdet[i][j] = M[i][j];

	//r = matrix_det(*fdet);

	//printf("Det: %.12f + I* %.12f\n", creal(r), cimag(r));

	info = matrix_inverse(*M);
	//matrix_print(M);
	//return r;
}

void update_inverse(int i, fmat M)
{
	if (up_counter < g_inverse) 
	{
		up_counter++;
		quick_update_inverse(i, M, temp1, temp2);
	}
	else 
	{
		up_counter = 0; //reset the counter
		hard_inverse(M);
	}
}


