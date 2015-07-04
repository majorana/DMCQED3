#include "measurement.h"
#include "lattice.h"
#include "fields.h"
#include "fermion.h"

complex double (*m_density)[GRIDPOINTS];

void measurement_init()
{
	measure_inter = 0;
	m_density = (fmat)malloc(g_measurements*GRIDPOINTS*sizeof(int));
}

complex double Wilsonloop(int i, int lx, int ly)
{
}

double density(fmat G)
{
	int i = 0;
	for(i = 0; i<GRIDPOINTS; i++)
	{
		m_density[measure_iter][i] = G[i, i];
	}
}

double density-density(int i, int j)
{

}

