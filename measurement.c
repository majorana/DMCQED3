#include <stdlib.h>
#include <stdio.h>
#include "measurement.h"
#include "lattice.h"
#include "fields.h"
#include "fermion.h"

complex double (*m_density)[Lx*Ly];
complex double (*m_polyakov)[Lx*Ly];

void polyakov_loop(int x, int y)
{
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

