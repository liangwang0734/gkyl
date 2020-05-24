#include <VlasovModDecl.h> 
double VlasovSurfElcMag3x3vMax_VX_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *EM, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w: Cell-center coordinates. dxv[NDIM]: Cell spacing. amax: amax in global lax flux. E: EM field. fl/fr: Distribution function in left/right cells 
// outl/outr: output distribution function in left/right cells 
// returns abs(amid) for use in determining amax in cfl and global lax flux 
  double dv10l = 2/dxvl[3]; 
  double dv10r = 2/dxvr[3]; 
  const double *E0 = &EM[0]; 

  const double dv1 = dxvr[3], wv1 = wr[3]; 
  const double dv2 = dxvr[4], wv2 = wr[4]; 
  const double dv3 = dxvr[5], wv3 = wr[5]; 

  const double *B0 = &EM[12]; 
  const double *B1 = &EM[16]; 
  const double *B2 = &EM[20]; 

  double Ghat[7]; 

  double alpha[7]; 

  double favg[7]; 

  favg[0] = 1*fr[0]+fl[0]; 
  favg[1] = 1*fr[1]+fl[1]; 
  favg[2] = 1*fr[2]+fl[2]; 
  favg[3] = 1*fr[3]+fl[3]; 
  favg[4] = -1*fr[4]+fl[4]; 
  favg[5] = 1*fr[5]+fl[5]; 
  favg[6] = 1*fr[6]+fl[6]; 
  double fjump[7]; 

  fjump[0] = amax*(1*fr[0]-fl[0]); 
  fjump[1] = amax*(1*fr[1]-fl[1]); 
  fjump[2] = amax*(1*fr[2]-fl[2]); 
  fjump[3] = amax*(1*fr[3]-fl[3]); 
  fjump[4] = amax*(-1*fr[4]-fl[4]); 
  fjump[5] = amax*(1*fr[5]-fl[5]); 
  fjump[6] = amax*(1*fr[6]-fl[6]); 
  alpha[0] = 2.828427124746191*(B2[0]*wv2+E0[0])-2.828427124746191*B1[0]*wv3; 
  alpha[1] = 2.828427124746191*(B2[1]*wv2+E0[1])-2.828427124746191*B1[1]*wv3; 
  alpha[2] = 2.828427124746191*(B2[2]*wv2+E0[2])-2.828427124746191*B1[2]*wv3; 
  alpha[3] = 2.828427124746191*(B2[3]*wv2+E0[3])-2.828427124746191*B1[3]*wv3; 
  alpha[5] = 0.8164965809277261*B2[0]*dv2; 
  alpha[6] = -0.8164965809277261*B1[0]*dv3; 
  const double amid = 0.125*alpha[0]; 
  Ghat[0] = 0.0625*(alpha[6]*favg[6]+alpha[5]*favg[5])-0.8660254037844386*fjump[4]+alpha[0]*(0.1082531754730548*favg[4]+0.0625*favg[0])+0.0625*(alpha[3]*favg[3]+alpha[2]*favg[2]+alpha[1]*favg[1])-0.5*fjump[0]; 
  Ghat[1] = alpha[1]*(0.1082531754730548*favg[4]+0.0625*favg[0])-0.5*fjump[1]+0.0625*alpha[0]*favg[1]; 
  Ghat[2] = alpha[2]*(0.1082531754730548*favg[4]+0.0625*favg[0])-0.5*fjump[2]+0.0625*alpha[0]*favg[2]; 
  Ghat[3] = alpha[3]*(0.1082531754730548*favg[4]+0.0625*favg[0])-0.5*fjump[3]+0.0625*alpha[0]*favg[3]; 
  Ghat[5] = (-0.5*fjump[5])+0.0625*alpha[0]*favg[5]+(0.1082531754730548*favg[4]+0.0625*favg[0])*alpha[5]; 
  Ghat[6] = (-0.5*fjump[6])+0.0625*alpha[0]*favg[6]+(0.1082531754730548*favg[4]+0.0625*favg[0])*alpha[6]; 

  outr[0] += 0.5*Ghat[0]*dv10r; 
  outr[1] += 0.5*Ghat[1]*dv10r; 
  outr[2] += 0.5*Ghat[2]*dv10r; 
  outr[3] += 0.5*Ghat[3]*dv10r; 
  outr[4] += -0.8660254037844386*Ghat[0]*dv10r; 
  outr[5] += 0.5*Ghat[5]*dv10r; 
  outr[6] += 0.5*Ghat[6]*dv10r; 

  outl[0] += -0.5*Ghat[0]*dv10l; 
  outl[1] += -0.5*Ghat[1]*dv10l; 
  outl[2] += -0.5*Ghat[2]*dv10l; 
  outl[3] += -0.5*Ghat[3]*dv10l; 
  outl[4] += -0.8660254037844386*Ghat[0]*dv10l; 
  outl[5] += -0.5*Ghat[5]*dv10l; 
  outl[6] += -0.5*Ghat[6]*dv10l; 
return std::abs(amid); 
} 
double VlasovSurfElcMag3x3vMax_VY_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *EM, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w: Cell-center coordinates. dxv[NDIM]: Cell spacing. amax: amax in global lax flux. E: EM field. fl/fr: Distribution function in left/right cells 
// outl/outr: output distribution function in left/right cells 
// returns abs(amid) for use in determining amax in cfl and global lax flux 
  double dv11l = 2/dxvl[4]; 
  double dv11r = 2/dxvr[4]; 
  const double *E1 = &EM[4]; 

  const double dv1 = dxvr[3], wv1 = wr[3]; 
  const double dv2 = dxvr[4], wv2 = wr[4]; 
  const double dv3 = dxvr[5], wv3 = wr[5]; 

  const double *B0 = &EM[12]; 
  const double *B1 = &EM[16]; 
  const double *B2 = &EM[20]; 

  double Ghat[7]; 

  double alpha[7]; 

  double favg[7]; 

  favg[0] = 1*fr[0]+fl[0]; 
  favg[1] = 1*fr[1]+fl[1]; 
  favg[2] = 1*fr[2]+fl[2]; 
  favg[3] = 1*fr[3]+fl[3]; 
  favg[4] = 1*fr[4]+fl[4]; 
  favg[5] = -1*fr[5]+fl[5]; 
  favg[6] = 1*fr[6]+fl[6]; 
  double fjump[7]; 

  fjump[0] = amax*(1*fr[0]-fl[0]); 
  fjump[1] = amax*(1*fr[1]-fl[1]); 
  fjump[2] = amax*(1*fr[2]-fl[2]); 
  fjump[3] = amax*(1*fr[3]-fl[3]); 
  fjump[4] = amax*(1*fr[4]-fl[4]); 
  fjump[5] = amax*(-1*fr[5]-fl[5]); 
  fjump[6] = amax*(1*fr[6]-fl[6]); 
  alpha[0] = 2.828427124746191*B0[0]*wv3-2.828427124746191*B2[0]*wv1+2.828427124746191*E1[0]; 
  alpha[1] = 2.828427124746191*B0[1]*wv3-2.828427124746191*B2[1]*wv1+2.828427124746191*E1[1]; 
  alpha[2] = 2.828427124746191*B0[2]*wv3-2.828427124746191*B2[2]*wv1+2.828427124746191*E1[2]; 
  alpha[3] = 2.828427124746191*B0[3]*wv3-2.828427124746191*B2[3]*wv1+2.828427124746191*E1[3]; 
  alpha[4] = -0.8164965809277261*B2[0]*dv1; 
  alpha[6] = 0.8164965809277261*B0[0]*dv3; 
  const double amid = 0.125*alpha[0]; 
  Ghat[0] = 0.0625*alpha[6]*favg[6]-0.8660254037844386*fjump[5]+alpha[0]*(0.1082531754730548*favg[5]+0.0625*favg[0])+0.0625*(alpha[4]*favg[4]+alpha[3]*favg[3]+alpha[2]*favg[2]+alpha[1]*favg[1])-0.5*fjump[0]; 
  Ghat[1] = alpha[1]*(0.1082531754730548*favg[5]+0.0625*favg[0])-0.5*fjump[1]+0.0625*alpha[0]*favg[1]; 
  Ghat[2] = alpha[2]*(0.1082531754730548*favg[5]+0.0625*favg[0])-0.5*fjump[2]+0.0625*alpha[0]*favg[2]; 
  Ghat[3] = alpha[3]*(0.1082531754730548*favg[5]+0.0625*favg[0])-0.5*fjump[3]+0.0625*alpha[0]*favg[3]; 
  Ghat[4] = alpha[4]*(0.1082531754730548*favg[5]+0.0625*favg[0])-0.5*fjump[4]+0.0625*alpha[0]*favg[4]; 
  Ghat[6] = (-0.5*fjump[6])+0.0625*alpha[0]*favg[6]+(0.1082531754730548*favg[5]+0.0625*favg[0])*alpha[6]; 

  outr[0] += 0.5*Ghat[0]*dv11r; 
  outr[1] += 0.5*Ghat[1]*dv11r; 
  outr[2] += 0.5*Ghat[2]*dv11r; 
  outr[3] += 0.5*Ghat[3]*dv11r; 
  outr[4] += 0.5*Ghat[4]*dv11r; 
  outr[5] += -0.8660254037844386*Ghat[0]*dv11r; 
  outr[6] += 0.5*Ghat[6]*dv11r; 

  outl[0] += -0.5*Ghat[0]*dv11l; 
  outl[1] += -0.5*Ghat[1]*dv11l; 
  outl[2] += -0.5*Ghat[2]*dv11l; 
  outl[3] += -0.5*Ghat[3]*dv11l; 
  outl[4] += -0.5*Ghat[4]*dv11l; 
  outl[5] += -0.8660254037844386*Ghat[0]*dv11l; 
  outl[6] += -0.5*Ghat[6]*dv11l; 
return std::abs(amid); 
} 
double VlasovSurfElcMag3x3vMax_VZ_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *EM, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w: Cell-center coordinates. dxv[NDIM]: Cell spacing. amax: amax in global lax flux. E: EM field. fl/fr: Distribution function in left/right cells 
// outl/outr: output distribution function in left/right cells 
// returns abs(amid) for use in determining amax in cfl and global lax flux 
  double dv12l = 2/dxvl[5]; 
  double dv12r = 2/dxvr[5]; 
  const double *E2 = &EM[8]; 

  const double dv1 = dxvr[3], wv1 = wr[3]; 
  const double dv2 = dxvr[4], wv2 = wr[4]; 
  const double dv3 = dxvr[5], wv3 = wr[5]; 

  const double *B0 = &EM[12]; 
  const double *B1 = &EM[16]; 
  const double *B2 = &EM[20]; 

  double Ghat[7]; 

  double alpha[7]; 

  double favg[7]; 

  favg[0] = 1*fr[0]+fl[0]; 
  favg[1] = 1*fr[1]+fl[1]; 
  favg[2] = 1*fr[2]+fl[2]; 
  favg[3] = 1*fr[3]+fl[3]; 
  favg[4] = 1*fr[4]+fl[4]; 
  favg[5] = 1*fr[5]+fl[5]; 
  favg[6] = -1*fr[6]+fl[6]; 
  double fjump[7]; 

  fjump[0] = amax*(1*fr[0]-fl[0]); 
  fjump[1] = amax*(1*fr[1]-fl[1]); 
  fjump[2] = amax*(1*fr[2]-fl[2]); 
  fjump[3] = amax*(1*fr[3]-fl[3]); 
  fjump[4] = amax*(1*fr[4]-fl[4]); 
  fjump[5] = amax*(1*fr[5]-fl[5]); 
  fjump[6] = amax*(-1*fr[6]-fl[6]); 
  alpha[0] = 2.828427124746191*(B1[0]*wv1+E2[0])-2.828427124746191*B0[0]*wv2; 
  alpha[1] = 2.828427124746191*(B1[1]*wv1+E2[1])-2.828427124746191*B0[1]*wv2; 
  alpha[2] = 2.828427124746191*(B1[2]*wv1+E2[2])-2.828427124746191*B0[2]*wv2; 
  alpha[3] = 2.828427124746191*(B1[3]*wv1+E2[3])-2.828427124746191*B0[3]*wv2; 
  alpha[4] = 0.8164965809277261*B1[0]*dv1; 
  alpha[5] = -0.8164965809277261*B0[0]*dv2; 
  const double amid = 0.125*alpha[0]; 
  Ghat[0] = (-0.8660254037844386*fjump[6])+alpha[0]*(0.1082531754730548*favg[6]+0.0625*favg[0])+0.0625*(alpha[5]*favg[5]+alpha[4]*favg[4]+alpha[3]*favg[3]+alpha[2]*favg[2]+alpha[1]*favg[1])-0.5*fjump[0]; 
  Ghat[1] = alpha[1]*(0.1082531754730548*favg[6]+0.0625*favg[0])-0.5*fjump[1]+0.0625*alpha[0]*favg[1]; 
  Ghat[2] = alpha[2]*(0.1082531754730548*favg[6]+0.0625*favg[0])-0.5*fjump[2]+0.0625*alpha[0]*favg[2]; 
  Ghat[3] = alpha[3]*(0.1082531754730548*favg[6]+0.0625*favg[0])-0.5*fjump[3]+0.0625*alpha[0]*favg[3]; 
  Ghat[4] = alpha[4]*(0.1082531754730548*favg[6]+0.0625*favg[0])-0.5*fjump[4]+0.0625*alpha[0]*favg[4]; 
  Ghat[5] = alpha[5]*(0.1082531754730548*favg[6]+0.0625*favg[0])-0.5*fjump[5]+0.0625*alpha[0]*favg[5]; 

  outr[0] += 0.5*Ghat[0]*dv12r; 
  outr[1] += 0.5*Ghat[1]*dv12r; 
  outr[2] += 0.5*Ghat[2]*dv12r; 
  outr[3] += 0.5*Ghat[3]*dv12r; 
  outr[4] += 0.5*Ghat[4]*dv12r; 
  outr[5] += 0.5*Ghat[5]*dv12r; 
  outr[6] += -0.8660254037844386*Ghat[0]*dv12r; 

  outl[0] += -0.5*Ghat[0]*dv12l; 
  outl[1] += -0.5*Ghat[1]*dv12l; 
  outl[2] += -0.5*Ghat[2]*dv12l; 
  outl[3] += -0.5*Ghat[3]*dv12l; 
  outl[4] += -0.5*Ghat[4]*dv12l; 
  outl[5] += -0.5*Ghat[5]*dv12l; 
  outl[6] += -0.8660254037844386*Ghat[0]*dv12l; 
return std::abs(amid); 
} 