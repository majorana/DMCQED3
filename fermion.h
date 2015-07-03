#ifndef _FERMION_H
#define _FERMION_H

/**********************************************************************************************/
/*** Implementation of the fermion single-particle operator
 * and quick updates of the propagator matrix
**********************************************************************************************/

#include "linalg.h"

extern void fermion(complex double *out, complex double *in);

complex double det_ratio_At(const int i);

complex double det_ratio_Axy(const int i); 

void update_inverse_At(const int i); 

void update_inverse_Axy(const int i); 

complex double get_fermion_mat();

#endif
