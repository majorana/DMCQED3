#ifndef _LATTICE_H
#define _LATTICE_H

/***********************************************/
/***** This unit defines lattice geometry  *****/
/***********************************************/    

#define Lx (6)                                       
#define Ly (6)                                      
#define Lt (6)
#define GRIDPOINTS (Lx*Ly*Lt)                           //Total number of lattice sites


extern int *xp;
extern int *xm;
extern int *yp;
extern int *ym;
extern int *tp;
extern int *tm;


extern int  init_lattice(const int Nx, const int Ny, const int Nt);      //This procedure initializes the above arrays

extern int idx(int it, int ix, int iy);

#endif
