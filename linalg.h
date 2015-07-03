#ifndef _LINALG_H
#define _LINALG_H

#include "complex/complex.h"

//If this is #defined, CG will print the norm of each residue
#undef MONITOR_CG_PROGRESS

/***********************************************************/
/**** Implementation of linear algebra on complex double fields ****/
/***********************************************************/


//Function type for linear operators acting on complex double vectors 
typedef void (*matrix_mult) (complex double * const, complex double * const, complex double * const);

void set_zero(complex double *P);                                     //Set the field P to zero

void assign_add_mul(complex double *P, complex double *Q, complex double c);  // P = P + c Q
void assign(complex double *R, complex double *S);                            // R = S
double scalar_prod_r(complex double *S, complex double *R);                   // Re(S*, R)
complex double scalar_prod(complex double *S, complex double *R);             // (S*, R)
void assign_mul_add_r(complex double *R, complex double *S, double c);        // R = c R + S, c is real
void assign_add_mul_r(complex double *P, complex double *Q, double c);        // P = P + c Q, c is real
void assign_add_mul(complex double *P, complex double *Q, complex double c);  // P = P + c Q
void assign_diff_mul(complex double *R, complex double *S, complex double c); // R = R - c S
void diff(complex double *Q, complex double *R, complex double *S);                   // Q = R - S
void mul_r(complex double *R, double c, complex double *S);                   // R = c S, c is real
void mul_c(complex double *R, complex double c, complex double *S);           // R = c S, c is complex
double square_norm(complex double *P);                                // (P, P*)
void add(complex double *Q, complex double *R, complex double *S);                    // Q = R + S

void print_vector(complex double *v);
void print_vector_r(double *v);
double max_r(double *v);

//Conjugate gradient method...
int cg(complex double *P, complex double *Q, int max_iter, double eps_sq, matrix_mult f);

//LAPACK wrapper
int matrix_inverse(complex double *mat);
int matrix_inverse_r(double *mat);
double matrix_det_r(double *mat);

#endif
