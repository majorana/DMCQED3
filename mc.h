#ifndef _HMC_H
#define _HMC_H
#include "complex/complex.h"
#include "lattice.h"

/***********************************************************************************/
/*** This unit implement the basic MC update step and the necessary procedures ****/
/***********************************************************************************/

/* Maximal number of iterations for CG method */
#define ITER_MAX 1000
/* Tolerance for CG method */
#define DELTACG 1.e-22

extern int R;  // Counter of all accepted configurations
extern int mc_iter;
extern complex double Minv[GRIDPOINTS][GRIDPOINTS]; // inverse of the fermion determinant, i.e. propagator

void mc_init();
void mc_update(); //Basic MC update step

#endif
