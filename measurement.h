#ifndef _H_MEASUREMENT
#define _H_MEASUREMENT

#include "complex/complex.h"
#include "lattice.h"
#include "fields.h"
#include "fermion.h"


/* Measurement of gauge-invariant quantites
 * Wilson loop, Polyakov loop
 * Fermion density, density-density correlation function
 * More general four point correlation functions: T_{ij}=c_i^\dag U_{ij} c_j, <T_{ij}T^\dag_{kl}>
 * 2k_F Friedel oscillations?
 * Pair correlation function c_{i,up}c_{i,down} is a gauge-invariant quantity because up and down have opposite charges. 
 */

extern int g_measurements;
extern int measure_iter;

void measurement_init();
void measurement_finish();

complex double Wilsonloop(int i, int lx, int ly);

void density(fmat G);

#endif
