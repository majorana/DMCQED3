#ifndef _FERMION_H
#define _FERMION_H

/**********************************************************************************************/
/*** Implementation of the fermion single-particle operator
 * and quick updates of the propagator matrix
**********************************************************************************************/

#include "linalg.h"

#define REINVERSE 40

typedef complex double (*fmat)[GRIDPOINTS];

extern double emu;

void fermion(complex double *out, complex double *in);

complex double det_ratio(const int i, fmat A); 
void update_row(const int i, fmat out, fmat in);
void update_col(const int i, fmat out, fmat in);
void quick_update_inverse(int i, fmat in, fmat temp);
void hard_inverse(fmat M);
void update_inverse(int i, fmat M, fmat temp);

extern int up_counter;

// for test
extern complex double Minv1[GRIDPOINTS][GRIDPOINTS];
extern complex double Minv2[GRIDPOINTS][GRIDPOINTS];
extern complex double Minv3[GRIDPOINTS][GRIDPOINTS];

void matrix_print(fmat A);

#endif
