#include <VlasovModDecl.h> 
__host__ __device__ void VlasovSurfStream1x2vSer_X_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells.
  double rdxl2 = 2.0/dxvl[0]; 
  double rdxr2 = 2.0/dxvr[0]; 

  double incr[8]; 

  if (wr[1]>0) { 
  incr[0] = dxvl[1]*(0.25*fl[4]+0.1443375672974065*fl[2])+(0.8660254037844386*fl[1]+0.5*fl[0])*wl[1]; 
  incr[1] = dxvl[1]*((-0.4330127018922193*fl[4])-0.25*fl[2])+((-1.5*fl[1])-0.8660254037844386*fl[0])*wl[1]; 
  incr[2] = wl[1]*(0.8660254037844386*fl[4]+0.5*fl[2])+dxvl[1]*(0.25*fl[1]+0.1443375672974065*fl[0]); 
  incr[3] = dxvl[1]*(0.25*fl[7]+0.1443375672974065*fl[6])+wl[1]*(0.8660254037844386*fl[5]+0.5*fl[3]); 
  incr[4] = wl[1]*((-1.5*fl[4])-0.8660254037844386*fl[2])+dxvl[1]*((-0.4330127018922193*fl[1])-0.25*fl[0]); 
  incr[5] = dxvl[1]*((-0.4330127018922193*fl[7])-0.25*fl[6])+wl[1]*((-1.5*fl[5])-0.8660254037844386*fl[3]); 
  incr[6] = wl[1]*(0.8660254037844386*fl[7]+0.5*fl[6])+dxvl[1]*(0.25*fl[5]+0.1443375672974065*fl[3]); 
  incr[7] = wl[1]*((-1.5*fl[7])-0.8660254037844386*fl[6])+dxvl[1]*((-0.4330127018922193*fl[5])-0.25*fl[3]); 

  outr[0] += incr[0]*rdxr2; 
  outr[1] += incr[1]*rdxr2; 
  outr[2] += incr[2]*rdxr2; 
  outr[3] += incr[3]*rdxr2; 
  outr[4] += incr[4]*rdxr2; 
  outr[5] += incr[5]*rdxr2; 
  outr[6] += incr[6]*rdxr2; 
  outr[7] += incr[7]*rdxr2; 

  outl[0] += -1.0*incr[0]*rdxl2; 
  outl[1] += incr[1]*rdxl2; 
  outl[2] += -1.0*incr[2]*rdxl2; 
  outl[3] += -1.0*incr[3]*rdxl2; 
  outl[4] += incr[4]*rdxl2; 
  outl[5] += incr[5]*rdxl2; 
  outl[6] += -1.0*incr[6]*rdxl2; 
  outl[7] += incr[7]*rdxl2; 
  } else { 
  incr[0] = dxvr[1]*(0.1443375672974065*fr[2]-0.25*fr[4])+(0.5*fr[0]-0.8660254037844386*fr[1])*wr[1]; 
  incr[1] = dxvr[1]*(0.4330127018922193*fr[4]-0.25*fr[2])+(1.5*fr[1]-0.8660254037844386*fr[0])*wr[1]; 
  incr[2] = wr[1]*(0.5*fr[2]-0.8660254037844386*fr[4])+dxvr[1]*(0.1443375672974065*fr[0]-0.25*fr[1]); 
  incr[3] = dxvr[1]*(0.1443375672974065*fr[6]-0.25*fr[7])+wr[1]*(0.5*fr[3]-0.8660254037844386*fr[5]); 
  incr[4] = wr[1]*(1.5*fr[4]-0.8660254037844386*fr[2])+dxvr[1]*(0.4330127018922193*fr[1]-0.25*fr[0]); 
  incr[5] = dxvr[1]*(0.4330127018922193*fr[7]-0.25*fr[6])+wr[1]*(1.5*fr[5]-0.8660254037844386*fr[3]); 
  incr[6] = wr[1]*(0.5*fr[6]-0.8660254037844386*fr[7])+dxvr[1]*(0.1443375672974065*fr[3]-0.25*fr[5]); 
  incr[7] = wr[1]*(1.5*fr[7]-0.8660254037844386*fr[6])+dxvr[1]*(0.4330127018922193*fr[5]-0.25*fr[3]); 

  outr[0] += incr[0]*rdxr2; 
  outr[1] += incr[1]*rdxr2; 
  outr[2] += incr[2]*rdxr2; 
  outr[3] += incr[3]*rdxr2; 
  outr[4] += incr[4]*rdxr2; 
  outr[5] += incr[5]*rdxr2; 
  outr[6] += incr[6]*rdxr2; 
  outr[7] += incr[7]*rdxr2; 

  outl[0] += -1.0*incr[0]*rdxl2; 
  outl[1] += incr[1]*rdxl2; 
  outl[2] += -1.0*incr[2]*rdxl2; 
  outl[3] += -1.0*incr[3]*rdxl2; 
  outl[4] += incr[4]*rdxl2; 
  outl[5] += incr[5]*rdxl2; 
  outl[6] += -1.0*incr[6]*rdxl2; 
  outl[7] += incr[7]*rdxl2; 
  } 
} 
