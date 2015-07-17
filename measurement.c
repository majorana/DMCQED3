#include <stdlib.h>
#include <stdio.h>
#include "measurement.h"
#include "lattice.h"
#include "fields.h"
#include "fermion.h"

complex double (*m_density_profile)[Lx*Ly];
complex double (*m_density);
complex double (*m_density_corr)[Lx][Ly];
complex double (*m_polyakov)[Lx*Ly];
complex double (*m_wilson)[Lx];
complex double (*m_wilson_xy)[Lx][Ly];

void measurement_init()
{
	measure_iter = 0;
	m_density_profile = malloc(g_measurements*Lx*Ly*sizeof(complex double));
	m_density = malloc(g_measurements*sizeof(complex double));
	m_density_corr = malloc(g_measurements*Lx*Ly*sizeof(complex double));
	m_wilson = malloc(g_measurements*(Lx/2)*sizeof(complex double));
	m_wilson_xy = malloc(g_measurements*(Lx/2)*(Ly/2)*sizeof(complex double));
}

void measurement_finish()
{
	free(m_density_profile);
	free(m_density);
	free(m_density_corr);
	free(m_wilson);
	free(m_wilson_xy);
}

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


complex double wilson_loop_xy1(int nx, int ny)
{
	// nx: spatial extension, say x direction
	// nt: time extension
	int i, j, kx, ky;
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
		// going from j+nx*e_x to j+nx*ex+ny*ey
		for(ky = 0; ky < ny; ky++) {
			w *= Uy[j];
			j = yp[j];
		}
		// going from j+nx*ex+ny*ey back to j+ny*ey
		for(kx = 0; kx < nx; kx++) {
			j = xm[j];
			w *= cconj(Ux[j]);
		}
		// going from j+nt*et back to j
		for(ky = 0; ky < ny; ky++) {
			j = ym[j];
			w *= cconj(Uy[j]);
		}
		avgw += w;
	}
	avgw /= GRIDPOINTS;
	//printf("Wilson loop: \t %.4f+I*%.4f\n", creal(avgw), cimag(avgw));
	return avgw;
}

complex double wilson_loop1(int nx, int nt)
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
			j = xm[j];
			w *= cconj(Ux[j]);
		}
		// going from j+nt*et back to j
		for(kt = 0; kt < nt; kt++) {
			j = tm[j];
			w *= cconj(Ut[j]);
		}
		avgw += w;
	}
	avgw /= GRIDPOINTS;
	//printf("Wilson loop: \t %.4f+I*%.4f\n", creal(avgw), cimag(avgw));
	return avgw;
}

// Wilson loop in the t-x plane
void wilson_loop(int nt) 
{
	int i;
	for (i = 0; i<Lx; i++)
	{
		m_wilson[measure_iter][i] = wilson_loop1(i, nt);
	}
}

// Wilson loop in the x-y plane
void wilson_loop_xy() 
{
	int i, j;
	for (i = 0; i<Lx; i++)
		for (j = 0; j<Ly; j++)
			m_wilson_xy[measure_iter][i][j] = wilson_loop_xy1(i, j);
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
			m_density_profile[measure_iter][s] = 0.0 + I*0.0;
			for (k=0; k<Lt;k++)
				m_density_profile[measure_iter][s] += G[idx(k, i, j)][idx(k, i, j)];
			m_density_profile[measure_iter][s] /= Lt;
			//printf("%.4f+I*%.4f, ", creal(m_density[measure_iter][s]), cimag(m_density[measure_iter][s]));
			avg += m_density_profile[measure_iter][s];
		}
	}
	m_density[measure_iter] = avg/(Lx*Ly);
	//printf("Average density: \t %.4f+I*%.4f\n", creal(avg/(Lx*Ly)), cimag(avg/(Lx*Ly)));
}

double mean_plaq()
{
	int i;
	double mp = 0.0;
	for (i = 0; i < GRIDPOINTS; i++) 
	{
		mp += S_Gxy(i);
	}
	return -mp/(beta*dt*GRIDPOINTS);
}

// <n_i n_j> = <c_i^\dag c_i c_j^\dag c_j> = <n_i><n_j> - <c_i^\dag c_j><c_j^\dag c_i>. The second term is the connected component 
// calculate 1/N\sum_x n(x) n(x+y), also average over time
// update(07/09): should get the full dynamical density-density response function?
void density_correlation(fmat G)
{
	int it, ix, iy, i, j;
	int nx, ny;
	complex double c;
	for(nx = 0; nx < Lx; nx++)
	{
		for(ny = 0; ny < Ly; ny++)
		{
			c = 0.0 + 0.0*I;
			for(i = 0; i<GRIDPOINTS; i++)
			{
				coordinate(i, &it, &ix, &iy);
				j = idx(it, ix + nx, iy + ny);
				c += G[i][j]*G[j][i];
			}
			m_density_corr[measure_iter][nx][ny] = c/(Lx*Ly*Lt);
		}
	}
	return;
}



