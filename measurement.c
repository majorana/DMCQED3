#include <stdlib.h>
#include <stdio.h>
#include "measurement.h"
#include "lattice.h"
#include "fields.h"
#include "fermion.h"
#include "mc.h"

complex double (*m_density_profile)[Lx*Ly];
complex double (*m_density);
complex double (*m_density_corr)[Lx][Ly];
complex double (*m_wilson_xy)[Lx][Ly];
complex double (*m_A)[3][GRIDPOINTS];

void wilson_loop_xy();
void density(fmat G);
void density_correlation(fmat G);

void measurement_init()
{
	measure_iter = 0;
	m_density_profile = malloc(g_measurements*Lx*Ly*sizeof(complex double));
	m_density = malloc(g_measurements*sizeof(complex double));
	m_density_corr = malloc(g_measurements*Lx*Ly*sizeof(complex double));
	m_wilson_xy = malloc(g_measurements*(Lx)*(Ly)*sizeof(complex double));
	m_A = malloc(g_measurements*Lx*Ly*3*sizeof(complex double));
}

void measurement_finish()
{
	free(m_density_profile);
	free(m_density);
	free(m_density_corr);
	free(m_wilson_xy);
	free(m_A);
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
void density_corr(fmat G)
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

void measure()
{
	int i;
	density(Minv);
	density_corr(Minv);
	for(i=0; i<GRIDPOINTS; i++) 
	{
		m_A[measure_iter][0][i] = At[i];
		m_A[measure_iter][1][i] = Ax[i];
		m_A[measure_iter][2][i] = Ay[i];
	}
	measure_iter++;
}

void output_measurement()
{
	int i, j, k;
	FILE *fp;

	fp = fopen("density.dat", "w");
	for(i = 0; i < g_measurements; i++)
		fprintf(fp, "%.5f %.5f\n", creal(m_density[i]), cimag(m_density[i]));
	fclose(fp);

	fp = fopen("wilsonxy.dat", "w");
	for(i = 0; i < g_measurements; i++)
	{
		for(j = 0; j < Lx/2; j++)
			for(k = 0; k < Ly/2; k++)
				fprintf(fp, "\t %.5f %.5f", creal(m_wilson_xy[i][j][k]), cimag(m_wilson_xy[i][j][k]));
		fprintf(fp, "\n");
	}
	fclose(fp);
	
	fp = fopen("density_corr.dat", "w");
	for(i = 0; i < g_measurements; i++)
	{
		for(j = 0; j < Lx; j++)
			for(k = 0; k < Ly; k++)
				fprintf(fp, "\t %.5f %.5f", creal(m_density_corr[i][j][k]), cimag(m_density_corr[i][j][k]) );
		fprintf(fp, "\n");
	}
	fclose(fp);
}


