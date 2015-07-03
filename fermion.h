#ifndef _FERMION_H
#define _FERMION_H

/**********************************************************************************************/
/*** Implementation of the fermion single-particle operator
 * and quick updates of the propagator matrix
**********************************************************************************************/

#include "linalg.h"

complex double (*Minv)[GRIDPOINTS];
complex double (*Minv_spare)[GRIDPOINTS];


extern void fermion(complex double *out, complex double *in);

complex double det_ratio(const int i); 

complex double get_fermion_mat();

#endif
