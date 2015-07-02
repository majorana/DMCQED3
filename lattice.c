#include <stdlib.h>
#include "lattice.h"

#define MOD(a, b) (((a) % (b)) + (b)) % (b)

/* index arrays for neighbours */
int * xp;
int * xm;
int * yp;
int * ym;
int * tp;
int * tm;

/* routine to initialise the geometry index routines */
/* takes care of the periodic boundary conditions    */
int init_lattice(const int Nx, const int Ny, const int Nt) {
  	int ix, iy, it, s;
  	div_t dg;
  	int N = Nx*Ny*Nt;

  	xp = (int *)malloc(N*sizeof(int));
  	xm = (int *)malloc(N*sizeof(int));
  	yp = (int *)malloc(N*sizeof(int));
  	ym = (int *)malloc(N*sizeof(int));
  	tp = (int *)malloc(N*sizeof(int));
  	tm = (int *)malloc(N*sizeof(int));

  	for(it = 0; it < Nt; it++) {
	  	for(ix = 0; ix < Nx; ix++) {
		  	for(iy = 0; iy < Ny; iy++) {
			  	s = it*Nx*Ny + ix*Ny + iy;
			  	xp[s] = it*Nx*Ny + MOD(ix+1, Nx)*Ny + iy;
				xm[s] = it*Nx*Ny + MOD(ix-1, Nx)*Ny + iy;
				yp[s] = it*Nx*Ny + ix*Ny + MOD(iy+1,Ny);
				ym[s] = it*Nx*Ny + ix*Ny + MOD(iy-1,Ny);
				tp[s] = MOD(it+1, Nt)*Nx*Ny + ix*Ny + iy;
				tm[s] = MOD(it-1, Nt)*Nx*Ny + ix*Ny + iy;
		  }
	  }
  }
  return(0);
}



