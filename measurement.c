#include <stdlib.h>
#include <stdio.h>
#include "measurement.h"
#include "lattice.h"
#include "fields.h"
#include "fermion.h"

complex double (*m_density)[Lx*Ly];
complex double (*m_polyakov)[Lx*Ly];

complex double polyakov_loop(int x, int y)
{
	// x, y are spatial coordinates
	int i;
	complex double Px;
	Px = 1.0 + 0.0*I;
	for(i = 0; i<Lt; i++)
		Px *= Ut[idx(i, x, y)];
	return Px;
}

complex double wilson_loop(int nx, int nt)
{
	// nx: spatial extension, say x direction
	// nt: time extension
	int i, j, kx, kt;
	complex double w, avgw;
	avgw = 0.0 + 0.0*I;
	for(i = 0; i < GRIDPOINTS; i++)
	{
		j = i;
		w = 1.0 + 0.0*I;
		// going from j to j+ nx*e_x
		for(kx = 0; kx < nx; kx++) {
			w *= Ux[j];
			j = xp[j];
		}
		// going from j+nx*e_x to j+nx*ex+nt*et
		for(kt = 0; kt < nt; kt++) {
			w *= Ut[j];
			j = tp[j];
		}
		// going from j+nx*ex+nt*et back to j+nt*et
		for(kx = 0; kx < nx; kx++) {
			w *= cconj(Ux[j]);
			j = xm[j];
		}
		// going from j+nt*et back to j
		for(kt = 0; kt < nt; kt++) {
			w *= cconj(Ut[j]);
			j = tm[j];
		}
		avgw += w;
	}
	return avgw/GRIDPOINTS;
}

void measurement_init()
{
	measure_iter = 0;
	m_density = malloc(g_measurements*Lx*Ly*sizeof(complex double));
}

void measurement_finish()
{
	free(m_density);
}

void density(fmat G)
// average over time
{
	int i, j, k, s;
	complex double avg;
	avg = 0.0 + I*0.0;
	for(i = 0; i<Lx; i++)
	{
		for(j = 0; j<Ly; j++) 
		{
			s = idx(0, i, j);
			m_density[measure_iter][s] = 0.0 + I*0.0;
			for (k=0; k<Lt;k++)
				m_density[measure_iter][s] += G[idx(k, i, j)][idx(k, i, j)];
			m_density[measure_iter][s] /= Lt;
			//printf("%.4f+I*%.4f, ", creal(m_density[measure_iter][s]), cimag(m_density[measure_iter][s]));
			avg += m_density[measure_iter][s];
		}
	}
	printf("\n Average over space-time: %.4f+I*%.4f ", creal(avg/(Lx*Ly)), cimag(avg/(Lx*Ly)));
	printf("\n");
}

// <n_i n_j> = <c_i^\dag c_i c_j^\dag c_j> = <n_i><n_j> - <c_i^\dag c_j><c_j^\dag c_i> 

void density_correlation(fmat G)
{
	return;
}

