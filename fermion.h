#ifndef _FERMION_H
#define _FERMION_H

/**********************************************************************************************/
/*** Implementation of the fermion single-particle operator
 * and quick updates of the propagator matrix
**********************************************************************************************/

#include "linalg.h"

#define REINVERSE 80

typedef complex double (*fmat)[GRIDPOINTS];

extern double emu;

void fermion(complex double *out, complex double *in);

complex double det_ratio_t(const int i, fmat A); 
complex double det_ratio_xy(const int i, fmat A); 

//void quick_update_inverse(int i, fmat in, fmat temp);
void hard_inverse(fmat M);
void update_inverse(int i, fmat M);

extern int up_counter;

// for test
extern complex double Minv1[GRIDPOINTS][GRIDPOINTS];
extern complex double Minv2[GRIDPOINTS][GRIDPOINTS];
extern complex double Minv3[GRIDPOINTS][GRIDPOINTS];

void matrix_print(fmat A);

#endif
