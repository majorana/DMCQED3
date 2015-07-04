#ifndef _FERMION_H
#define _FERMION_H

/**********************************************************************************************/
/*** Implementation of the fermion single-particle operator
 * and quick updates of the propagator matrix
**********************************************************************************************/

#include "linalg.h"

typedef complex double (*fmat)[GRIDPOINTS];

void fermion(complex double *out, complex double *in);

complex double det_ratio(const int i, fmat A); 
void update_row(const int i, fmat out, fmat in);
void update_col(const int i, fmat out, fmat in);
void quick_update_inverse(int i, fmat in, fmat temp);

void hard_inverse(fmat M);

extern complex double Minv1[GRIDPOINTS][GRIDPOINTS];
extern complex double Minv2[GRIDPOINTS][GRIDPOINTS];
extern complex double Minv3[GRIDPOINTS][GRIDPOINTS];

extern int up_counter;

#endif
