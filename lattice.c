#include <stdlib.h>
#include "lattice.h"

#define MOD(a, b) (((a) % (b)) + (b)) % (b)

/* index arrays for neighbours */
int xp[GRIDPOINTS];
int xm[GRIDPOINTS];
int yp[GRIDPOINTS];
int ym[GRIDPOINTS]; 
int tp[GRIDPOINTS];
int tm[GRIDPOINTS];

/* routine to initialise the geometry index routines */
/* takes care of the periodic boundary conditions    */
int init_lattice() 
{
  	int ix, iy, it, s;

  	for(it = 0; it < Lt; it++) {
	  	for(ix = 0; ix < Lx; ix++) {
		  	for(iy = 0; iy < Ly; iy++) {
			  	s = it*Lx*Ly + ix*Ly + iy;
			  	xp[s] = it*Lx*Ly + MOD(ix+1, Lx)*Ly + iy;
				xm[s] = it*Lx*Ly + MOD(ix-1, Lx)*Ly + iy;
				yp[s] = it*Lx*Ly + ix*Ly + MOD(iy+1, Ly);
				ym[s] = it*Lx*Ly + ix*Ly + MOD(iy-1, Ly);
				tp[s] = MOD(it+1, Lt)*Lx*Ly + ix*Ly + iy;
				tm[s] = MOD(it-1, Lt)*Lx*Ly + ix*Ly + iy;
		  }
	  }
  }
  return(0);
}

int idx(int it, int ix, int iy) {
	return MOD(it, Lt)*Lx*Ly + MOD(ix, Lx)*Ly + MOD(iy, Ly);
}


