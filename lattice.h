#ifndef _LATTICE_H
#define _LATTICE_H

/***********************************************/
/***** This unit defines lattice geometry  *****/
/***********************************************/    

#define Lx (6)                                       
#define Ly (6)                                      
#define Lt (6)
#define GRIDPOINTS (Lx*Ly*Lt)                           //Total number of lattice sites

extern int xp[GRIDPOINTS];
extern int xm[GRIDPOINTS];
extern int yp[GRIDPOINTS];
extern int ym[GRIDPOINTS]; 
extern int tp[GRIDPOINTS];
extern int tm[GRIDPOINTS];


int init_lattice();      //This procedure initializes the above arrays

int idx(int it, int ix, int iy);

void coordinate(int i, int *it, int *ix, int *iy);

#endif
