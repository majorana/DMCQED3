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
 * Pair correlation function c_{i,up}c_{i,down} is a gauge-invariant quantity because up and down have opposite charges. <c_{i up}^\dag c_{i down}^\dag c_{j down} c_{j up}>=<c_{i up}^\dag c_{j up}><c_{i down}^\dag c_{j down}> = |<c_{i up}^\dag c_{j up}>|^2
 */

extern int g_measurements;
extern int measure_iter;

extern complex double (*m_density);
extern complex double (*m_density_corr)[Lx][Ly];
extern complex double (*m_wilson)[Lx/2];


void measurement_init();
void measurement_finish();

void wilson_loop(int nt);
void density(fmat G);
void density_correlation(fmat G);
double mean_plaq();

#endif
