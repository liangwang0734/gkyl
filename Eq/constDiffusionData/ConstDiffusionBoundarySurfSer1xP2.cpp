#include <ConstDiffusionModDecl.h> 
void ConstDiffusionBoundarySurf1xSerP2_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[1]:      Cell-center coordinates.
  // dx[1]:     Cell spacing.
  // idx[1]:    current grid index.
  // nu[1]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 4.0*nu[0]/(dxl[0]*dxl[0]); 
  double rdxFnur = 4.0*nu[0]/(dxr[0]*dxr[0]); 

  double incr1[3]; 
  double incr2[3]; 

  if (idxr[0] == 1) {

  incr2[1] = 1.936491673103708*fr[2]-1.5*fr[1]+0.8660254037844386*fr[0]; 
  incr2[2] = (-7.5*fr[2])+5.809475019311125*fr[1]-3.354101966249685*fr[0]; 

  outr[1] += incr2[1]*rdxFnur; 
  outr[2] += incr2[2]*rdxFnur; 

  } else {

  incr2[1] = 1.936491673103708*fl[2]+1.5*fl[1]+0.8660254037844386*fl[0]; 
  incr2[2] = (-7.5*fl[2])-5.809475019311125*fl[1]-3.354101966249685*fl[0]; 

  outl[1] += -1.0*incr2[1]*rdxFnul; 
  outl[2] += incr2[2]*rdxFnul; 

  }

} 
void ConstHyperDiffusion4BoundarySurf1xSerP2_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[1]:      Cell-center coordinates.
  // dx[1]:     Cell spacing.
  // idx[1]:    current grid index.
  // nu[1]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 16.0*nu[0]/(dxl[0]*dxl[0]*dxl[0]*dxl[0]); 
  double rdxFnur = 16.0*nu[0]/(dxr[0]*dxr[0]*dxr[0]*dxr[0]); 

  double incr1[3]; 
  double incr2[3]; 
  double incr3[3]; 
  double incr4[3]; 

  if (idxr[0] == 1) {


  incr2[1] = 2.8125*fr[1]-5.083290641897234*fr[2]; 
  incr2[2] = 19.6875*fr[2]-10.89276566120836*fr[1]; 



  outr[1] += -1.0*incr2[1]*rdxFnur; 
  outr[2] += -1.0*incr2[2]*rdxFnur; 

  } else {


  incr2[1] = 5.083290641897234*fl[2]+2.8125*fl[1]; 
  incr2[2] = 19.6875*fl[2]+10.89276566120836*fl[1]; 



  outl[1] += -1.0*incr2[1]*rdxFnul; 
  outl[2] += -1.0*incr2[2]*rdxFnul; 

  }

} 
void ConstHyperDiffusion6BoundarySurf1xSerP2_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[1]:      Cell-center coordinates.
  // dx[1]:     Cell spacing.
  // idx[1]:    current grid index.
  // nu[1]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 64.0*nu[0]/(dxl[0]*dxl[0]*dxl[0]*dxl[0]*dxl[0]*dxl[0]); 
  double rdxFnur = 64.0*nu[0]/(dxr[0]*dxr[0]*dxr[0]*dxr[0]*dxr[0]*dxr[0]); 

  double incr1[3]; 
  double incr2[3]; 
  double incr3[3]; 
  double incr4[3]; 
  double incr5[3]; 
  double incr6[3]; 

  if (idxr[0] == 1) {


  incr2[1] = 19.06233990711463*fr[2]-4.921875*fr[1]; 
  incr2[2] = 19.06233990711463*fr[1]-73.828125*fr[2]; 





  outr[1] += incr2[1]*rdxFnur; 
  outr[2] += incr2[2]*rdxFnur; 

  } else {


  incr2[1] = (-19.06233990711463*fl[2])-4.921875*fl[1]; 
  incr2[2] = (-73.828125*fl[2])-19.06233990711463*fl[1]; 





  outl[1] += incr2[1]*rdxFnul; 
  outl[2] += incr2[2]*rdxFnul; 

  }

} 
