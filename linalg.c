#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#ifdef __APPLE__ // on laptop
	#include "lapack/lapacke.h"
#elif __linux // on linux machine with MKL, Q1 cluster 
	#include <mkl.h>
#endif

#include "complex/complex.h"
#include "lattice.h"
#include "linalg.h"

void add(complex double *Q, complex double *R, complex double *S)
{
  	int ix;
  	for (ix=0;ix<GRIDPOINTS;ix++) {
    	Q[ix] = R[ix] + S[ix];
  	}
}

double square_norm(complex double *S)
{
  	int ix;
  	static double ds;

  	ds=0.0;
  	for (ix=0;ix<GRIDPOINTS;ix++) {
    	ds += creal(cconj(S[ix])*S[ix]);
  	}
  	return ds;
}

void assign(complex double *R, complex double *S)
{
  	int ix;
  	for (ix=0;ix<GRIDPOINTS;ix++){
    	R[ix] = S[ix];
  	}
}

void assign_add_mul(complex double *P, complex double *Q, complex double c)
{
  	int ix;
  	for (ix=0;ix<GRIDPOINTS;ix++){
		P[ix] = P[ix] + c*Q[ix];
  	  	//(*r).s1 = (*r).s1 + c*(*s).s1;
  	}
}

// P = P + c Q, c is real
void assign_add_mul_r(complex double *P, complex double *Q, double c)
{
  	int ix;
  	static double fact;

  	fact=c;
   
  	for (ix=0;ix<GRIDPOINTS;ix++){
    	P[ix] = P[ix] + c*Q[ix];
  	}
}

void assign_diff_mul(complex double *R, complex double *S, complex double c){
  	int ix;
  	for (ix=0;ix<GRIDPOINTS;ix++) {
    	R[ix]=R[ix]-c*S[ix];
  	}
}

void assign_mul_add_r(complex double *R, complex double *S, double c)
{
  	int ix;
  	static double fact;
  
  	fact=c;
  
  	for (ix=0;ix<GRIDPOINTS;ix++){
    	R[ix] = c*R[ix] + S[ix];
    //(*r).s2=fact*(*r).s2+(*s).s2;
  	}
}

void diff(complex double *Q, complex double *R, complex double *S)
{
  	int ix;
  	for (ix=0;ix<GRIDPOINTS;ix++){
		Q[ix] = R[ix]-S[ix];
    //(*q).s2=(*r).s2-(*s).s2;
  	}
}

void mul_r(complex double *R, double c, complex double *S)
{
  	int ix;

  	for (ix=0;ix<GRIDPOINTS;ix++){
    	R[ix]=c*S[ix]; 
    //(*r).s2=c*(*s).s2;
  	}
}

void mul_c(complex double *R, complex double c, complex double *S)
{
  	int ix;
  	for (ix=0;ix<GRIDPOINTS;ix++){
    	R[ix]=c*S[ix]; 
    //(*r).s2=c*(*s).s2;
  	}
}

complex double scalar_prod(complex double *S, complex double *R)
{
  	int ix;
  	static complex double ds;
  
  	/* Real Part */

  	ds=0.0 + I*0.0;
  
  	for (ix=0;ix<GRIDPOINTS;ix++){
    	ds+=cconj(S[ix])*R[ix]; //+cconj((*s).s2)*(*r).s2;
  	}

  	return(ds);
}

double scalar_prod_r(complex double *S, complex double *R)
{
  	int ix;
  	static double ds;
  
  	/* Real Part */

  	ds=0.0;
  
  	for (ix=0;ix<GRIDPOINTS;ix++){
    //ds+=conj((*s).s1)*(*r).s1+conj((*s).s2)*(*r).s2;
    	ds += creal(cconj(S[ix])*R[ix]);
  	}

  	return(ds);
}

void print_vector(complex double *v)
{
	int i;
	for(i = 0; i<GRIDPOINTS; i++) {
		printf("(%f,%f),", creal(v[i]), cimag(v[i]));
	}
	printf("\n\n");
}

void print_vector_r(double *v)
{
	int i;
	for(i = 0; i<GRIDPOINTS; i++) {
		printf("%.3f,", v[i]);
	}
	printf("\n");
}

double max_r(double *v)
{
	int i;
	double m;
	m = abs(v[0]);
	for(i=0; i<GRIDPOINTS;i++) {
		if (abs(v[i]) > m) {
			m = abs(v[i]);
		}
	}
	return(m);
}

/************ Conjugate gradient ****************/
/***    Solves the equation f*P = Q           ***/
/************************************************/

int cg(complex double *P, complex double *Q, int max_iter, double eps_sq, matrix_mult f)
{
 	double normsq, pro, err, alpha_cg, beta_cg;
 	int iteration;
 	complex double r[GRIDPOINTS], p[GRIDPOINTS];
 	complex double x[GRIDPOINTS], q2p[GRIDPOINTS];
 	complex double tmp1[GRIDPOINTS];
 

 	/* Initial guess for the solution - zero works well here */ 
 	set_zero(x);
  
 	/* initialize residue r and search vector p */
 	assign(r, Q); /* r = Q - f*x, x=0 */
 	assign(p, r);
 	normsq = square_norm(r);
  
 	/* main loop */
#ifdef MONITOR_CG_PROGRESS
  	printf("\n\n Starting CG iterations...\n");
#endif
 	for(iteration=0; iteration<max_iter; iteration++)
 	{
  		f(q2p, tmp1, p);
	
  		pro = scalar_prod_r(p, q2p);
  		/*  Compute alpha_cg(i+1)   */
  		alpha_cg = normsq/pro;
  		/*  Compute x_(i+1) = x_i + alpha_cg(i+1) p_i    */
  		assign_add_mul_r(x, p,  alpha_cg);
  		/*  Compute r_(i+1) = r_i - alpha_cg(i+1) Q2p_i   */
  		assign_add_mul_r(r, q2p, -alpha_cg);
  		/* Check whether the precision is reached ... */
  		err=square_norm(r);
#ifdef MONITOR_CG_PROGRESS
  		printf("\t CG iteration %i, alpha = %.4f, |r|^2 = %.2f\n", iteration, alpha_cg, err);
#endif
  		if(err <= eps_sq)
  		{
#ifdef MONITOR_CG_PROGRESS
   			printf("Required precision reached, stopping CG iterations...\n\n");
#endif
   			assign(P, x);
   			return(iteration);
  		};
     
  		/* Compute beta_cg(i+1) */
  		beta_cg = err/normsq;
  		/* Compute p_(i+1) = r_i+1 + beta_(i+1) p_i     */
  		assign_mul_add_r(p, r, beta_cg);
  		normsq = err;
 	}
	fprintf(stderr, "WARNING: CG didn't converge after %d iterations!\n", max_iter);
	return (-1);
}


void set_zero(complex double *P)
{
 	int i;
 	for(i=0; i<GRIDPOINTS; i++)
 	{
  		P[i] = 0.0 + I*0.0;
 	};
}

/* *************************************************************
 * LAPACKE routines wrapper: matrix inverse and determinant
 * *************************************************************
 */

// Notice that all these LAPACK routines use the array passed to them to store results or do some other things. Lots of side-effects.

// We let LAPACKE to handle workspace allocation.

// Need dereference if a 2D array, e.g. a[][] is passed, i.e. matrix_inverse_r(*a)
int matrix_inverse_r(double *mat)
{
   	lapack_int n, lda, info1, info2;
   	lapack_int ipiv[GRIDPOINTS];

   	n = GRIDPOINTS;
   	lda = GRIDPOINTS;

   	info1 = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, n, n, mat, lda, ipiv);
	info2 = LAPACKE_dgetri(LAPACK_ROW_MAJOR, n, mat, lda, ipiv);
	
	if (info1 || info2) 
		return 1; // singular?
	else
		return 0; // correct
}

int matrix_inverse(complex double *mat)
{
   	lapack_int n, lda, info1, info2;
   	lapack_int ipiv[GRIDPOINTS];

   	n = GRIDPOINTS;
   	lda = GRIDPOINTS;

   	info1 = LAPACKE_zgetrf(LAPACK_ROW_MAJOR, n, n, mat, lda, ipiv);
	info2 = LAPACKE_zgetri(LAPACK_ROW_MAJOR, n, mat, lda, ipiv);

	if (info1 || info2) 
		return 1;
	else
		return 0;
}

double matrix_det_r(double *mat)
// overall sign is not calculated
{
	int i;
	double det = 1.0;
	lapack_int n, lda, info;
   	lapack_int ipiv[GRIDPOINTS];

   	n = GRIDPOINTS;
   	lda = GRIDPOINTS;

   	info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, n, n, mat, lda, ipiv);

	for(i = 0; i<n; i++) {
		det *= mat[i*GRIDPOINTS + i];
	}
	return det;
}

complex double matrix_det(complex double *mat)
// overall sign is not calculated
{
	int i;
	complex double det = 1.0;
	lapack_int n, lda, info;
   	lapack_int ipiv[GRIDPOINTS];

   	n = GRIDPOINTS;
   	lda = GRIDPOINTS;

   	info = LAPACKE_zgetrf(LAPACK_ROW_MAJOR, n, n, mat, lda, ipiv);

	for(i = 0; i<n; i++) {
		//printf("(%.4f)+I*(%.4f),\n", creal(mat[i*GRIDPOINTS + i]), cimag(mat[i*GRIDPOINTS + i]));
		det = det*mat[i*GRIDPOINTS + i];
	}
	return det;
}
