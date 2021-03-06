#include <GyrokineticGenGeoModDecl.h> 
double EmGyrokineticGenGeoVol2x2vSerP1_Bvars_0(const double q_, const double m_, const double *w, const double *dxv, const double *Bmag, const double *jacobTotInv, const double *Gradpar, const double *geoX, const double *geoY, const double *geoZ, const double *Phi, const double *Apar, const double* dApardt, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_y = 2.0/dxv[1]; 
  double dfac_v = 2.0/dxv[2]; 
  double dfac_m = 2.0/dxv[3]; 
  double wx = w[0]; 
  double wy = w[1]; 
  double wv = w[2]; 
  double wm = w[3]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double m2 = m_*m_; 
  double hamil[16]; 
  hamil[0] = (0.6666666666666666*(3.0*dfac_v2*(m_*wv2+Bmag[0]*wm+Phi[0]*q_)+m_))/dfac_v2; 
  hamil[1] = 2.0*Phi[1]*q_; 
  hamil[2] = 2.0*Phi[2]*q_; 
  hamil[3] = (2.309401076758503*m_*wv)/dfac_v; 
  hamil[4] = (1.154700538379252*Bmag[0])/dfac_m; 
  hamil[5] = 2.0*Phi[3]*q_; 
  double BstarX_by_Bmag[16]; 
  BstarX_by_Bmag[0] = 0.8660254037844386*geoZ[0]*jacobTotInv[0]*Apar[2]*dfac_y; 
  BstarX_by_Bmag[1] = 0.8660254037844386*geoZ[0]*jacobTotInv[0]*Apar[3]*dfac_y; 

  double BstarY_by_Bmag[16]; 
  BstarY_by_Bmag[0] = -0.8660254037844386*geoZ[0]*jacobTotInv[0]*Apar[1]*dfac_x; 
  BstarY_by_Bmag[2] = -0.8660254037844386*geoZ[0]*jacobTotInv[0]*Apar[3]*dfac_x; 

  double cflFreq = 0.0; 
  double alphaL = 0.0; 
  double alphaR = 0.0; 
  double alphax[16]; 
  alphax[0] = 1.732050807568877*dfac_x*((0.25*BstarX_by_Bmag[0]*hamil[3]*dfac_v)/m_-(0.25*geoZ[0]*jacobTotInv[0]*hamil[2]*dfac_y)/q_); 
  alphax[1] = 1.732050807568877*dfac_x*((0.25*BstarX_by_Bmag[1]*hamil[3]*dfac_v)/m_-(0.25*geoZ[0]*jacobTotInv[0]*hamil[5]*dfac_y)/q_); 
#if cflType == SURFAVG 
  // evaluate surface-averaged alpha on left 
  alphaL = -0.125*(1.732050807568877*alphax[1]-1.0*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // evaluate surface-averaged alpha on right 
  alphaR = 0.125*(1.732050807568877*alphax[1]+alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#elif cflType == QUAD 
  // evaluate alpha at left surface quadrature points 
  alphaL = 0.0625*(0.25*alphax[0]-0.4330127018922193*alphax[1]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alphax[0]-0.4330127018922193*alphax[1]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alphax[0]-0.4330127018922193*alphax[1]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alphax[0]-0.4330127018922193*alphax[1]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alphax[0]-0.4330127018922193*alphax[1]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alphax[0]-0.4330127018922193*alphax[1]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alphax[0]-0.4330127018922193*alphax[1]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alphax[0]-0.4330127018922193*alphax[1]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // evaluate alpha at right surface quadrature points 
  alphaR = 0.0625*(0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#endif 

  double alphay[16]; 
  alphay[0] = 0.4330127018922193*dfac_y*((geoZ[0]*jacobTotInv[0]*hamil[1]*dfac_x)/q_+(BstarY_by_Bmag[0]*hamil[3]*dfac_v)/m_); 
  alphay[2] = 0.4330127018922193*dfac_y*((geoZ[0]*jacobTotInv[0]*hamil[5]*dfac_x)/q_+(BstarY_by_Bmag[2]*hamil[3]*dfac_v)/m_); 
#if cflType == SURFAVG 
  // evaluate surface-averaged alpha on left 
  alphaL = -0.125*(1.732050807568877*alphay[2]-1.0*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // evaluate surface-averaged alpha on right 
  alphaR = 0.125*(1.732050807568877*alphay[2]+alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#elif cflType == QUAD 
  // evaluate alpha at left surface quadrature points 
  alphaL = 0.0625*(0.25*alphay[0]-0.4330127018922193*alphay[2]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alphay[0]-0.4330127018922193*alphay[2]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alphay[0]-0.4330127018922193*alphay[2]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alphay[0]-0.4330127018922193*alphay[2]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alphay[0]-0.4330127018922193*alphay[2]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alphay[0]-0.4330127018922193*alphay[2]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alphay[0]-0.4330127018922193*alphay[2]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alphay[0]-0.4330127018922193*alphay[2]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // evaluate alpha at right surface quadrature points 
  alphaR = 0.0625*(0.4330127018922193*alphay[2]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alphay[2]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alphay[2]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alphay[2]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alphay[2]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alphay[2]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alphay[2]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alphay[2]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#endif 

  double alphav[16]; 
  alphav[0] = -(0.4330127018922193*dfac_v*(BstarY_by_Bmag[0]*hamil[2]*dfac_y+BstarX_by_Bmag[0]*hamil[1]*dfac_x))/m_; 
  alphav[1] = -(0.4330127018922193*dfac_v*(BstarY_by_Bmag[0]*hamil[5]*dfac_y+BstarX_by_Bmag[1]*hamil[1]*dfac_x))/m_; 
  alphav[2] = -(0.4330127018922193*dfac_v*(BstarY_by_Bmag[2]*hamil[2]*dfac_y+BstarX_by_Bmag[0]*hamil[5]*dfac_x))/m_; 
  alphav[5] = -(0.4330127018922193*hamil[5]*dfac_v*(BstarY_by_Bmag[2]*dfac_y+BstarX_by_Bmag[1]*dfac_x))/m_; 
#if cflType == SURFAVG 
  // evaluate surface-averaged alpha on left 
  alphaL = 0.125*alphav[0]; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // evaluate surface-averaged alpha on right 
  alphaR = 0.125*alphav[0]; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#elif cflType == QUAD 
  // evaluate alpha at left surface quadrature points 
  alphaL = 0.0625*(0.25*alphav[5]-0.25*(alphav[2]+alphav[1])+0.25*alphav[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*(alphav[1]+alphav[0])-0.25*(alphav[5]+alphav[2])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.25*alphav[5])+0.25*alphav[2]-0.25*alphav[1]+0.25*alphav[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.015625*(alphav[5]+alphav[2]+alphav[1]+alphav[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alphav[5]-0.25*(alphav[2]+alphav[1])+0.25*alphav[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*(alphav[1]+alphav[0])-0.25*(alphav[5]+alphav[2])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.25*alphav[5])+0.25*alphav[2]-0.25*alphav[1]+0.25*alphav[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.015625*(alphav[5]+alphav[2]+alphav[1]+alphav[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // evaluate alpha at right surface quadrature points 
  alphaR = 0.0625*(0.25*alphav[5]-0.25*(alphav[2]+alphav[1])+0.25*alphav[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.25*(alphav[1]+alphav[0])-0.25*(alphav[5]+alphav[2])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.25*alphav[5])+0.25*alphav[2]-0.25*alphav[1]+0.25*alphav[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.015625*(alphav[5]+alphav[2]+alphav[1]+alphav[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.25*alphav[5]-0.25*(alphav[2]+alphav[1])+0.25*alphav[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.25*(alphav[1]+alphav[0])-0.25*(alphav[5]+alphav[2])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.25*alphav[5])+0.25*alphav[2]-0.25*alphav[1]+0.25*alphav[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.015625*(alphav[5]+alphav[2]+alphav[1]+alphav[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#endif 

  out[1] += 0.4330127018922193*(alphax[1]*f[1]+alphax[0]*f[0]); 
  out[2] += 0.4330127018922193*(alphay[2]*f[2]+alphay[0]*f[0]); 
  out[3] += 0.4330127018922193*(alphav[5]*f[5]+alphav[2]*f[2]+alphav[1]*f[1]+alphav[0]*f[0]); 
  out[5] += 0.4330127018922193*((alphay[2]+alphax[1])*f[5]+alphax[0]*f[2]+alphay[0]*f[1]); 
  out[6] += 0.4330127018922193*(alphax[1]*f[6]+alphav[2]*f[5]+f[2]*alphav[5]+alphax[0]*f[3]+alphav[0]*f[1]+f[0]*alphav[1]); 
  out[7] += 0.4330127018922193*(alphay[2]*f[7]+alphav[1]*f[5]+f[1]*alphav[5]+alphay[0]*f[3]+alphav[0]*f[2]+f[0]*alphav[2]); 
  out[8] += 0.4330127018922193*(alphax[1]*f[8]+alphax[0]*f[4]); 
  out[9] += 0.4330127018922193*(alphay[2]*f[9]+alphay[0]*f[4]); 
  out[10] += 0.4330127018922193*(alphav[5]*f[12]+alphav[2]*f[9]+alphav[1]*f[8]+alphav[0]*f[4]); 
  out[11] += 0.4330127018922193*((alphay[2]+alphax[1])*f[11]+alphax[0]*f[7]+alphay[0]*f[6]+alphav[0]*f[5]+f[0]*alphav[5]+alphav[1]*f[2]+f[1]*alphav[2]); 
  out[12] += 0.4330127018922193*((alphay[2]+alphax[1])*f[12]+alphax[0]*f[9]+alphay[0]*f[8]); 
  out[13] += 0.4330127018922193*(alphax[1]*f[13]+alphav[2]*f[12]+alphax[0]*f[10]+alphav[5]*f[9]+alphav[0]*f[8]+alphav[1]*f[4]); 
  out[14] += 0.4330127018922193*(alphay[2]*f[14]+alphav[1]*f[12]+alphay[0]*f[10]+alphav[0]*f[9]+alphav[5]*f[8]+alphav[2]*f[4]); 
  out[15] += 0.4330127018922193*((alphay[2]+alphax[1])*f[15]+alphax[0]*f[14]+alphay[0]*f[13]+alphav[0]*f[12]+alphav[1]*f[9]+alphav[2]*f[8]+f[4]*alphav[5]); 
  return cflFreq; 
} 
double EmGyrokineticGenGeoStep2Vol2x2vSerP1(const double q_, const double m_, const double *w, const double *dxv, const double *dApardt, const double *f, double *out) 
{ 
  double dvInv = 1.0/dxv[2]; 
  double dfac_v = 2.0/dxv[2]; 
  out[3] += -(0.8660254037844386*(dApardt[3]*f[5]+dApardt[2]*f[2]+dApardt[1]*f[1]+dApardt[0]*f[0])*dfac_v*q_)/m_; 
  out[6] += -(0.8660254037844386*(dApardt[2]*f[5]+f[2]*dApardt[3]+dApardt[0]*f[1]+f[0]*dApardt[1])*dfac_v*q_)/m_; 
  out[7] += -(0.8660254037844386*(dApardt[1]*f[5]+f[1]*dApardt[3]+dApardt[0]*f[2]+f[0]*dApardt[2])*dfac_v*q_)/m_; 
  out[10] += -(0.8660254037844386*(dApardt[3]*f[12]+dApardt[2]*f[9]+dApardt[1]*f[8]+dApardt[0]*f[4])*dfac_v*q_)/m_; 
  out[11] += -(0.8660254037844386*(dApardt[0]*f[5]+f[0]*dApardt[3]+dApardt[1]*f[2]+f[1]*dApardt[2])*dfac_v*q_)/m_; 
  out[13] += -(0.8660254037844386*(dApardt[2]*f[12]+dApardt[3]*f[9]+dApardt[0]*f[8]+dApardt[1]*f[4])*dfac_v*q_)/m_; 
  out[14] += -(0.8660254037844386*(dApardt[1]*f[12]+dApardt[0]*f[9]+dApardt[3]*f[8]+dApardt[2]*f[4])*dfac_v*q_)/m_; 
  out[15] += -(0.8660254037844386*(dApardt[0]*f[12]+dApardt[1]*f[9]+dApardt[2]*f[8]+dApardt[3]*f[4])*dfac_v*q_)/m_; 
  double cflFreq = 0.0; 
  double alphaL = 0.0; 
  double alphaR = 0.0; 
#if cflType == SURFAVG 
  // evaluate surface-averaged alpha on left 
  alphaL = -(0.25*dApardt[0]*dfac_v*q_)/m_; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // evaluate surface-averaged alpha on right 
  alphaR = -(0.25*dApardt[0]*dfac_v*q_)/m_; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#elif cflType == QUAD 
  // evaluate alpha at left surface quadrature points 
  alphaL = -(0.0625*(0.5*dApardt[3]-0.5*(dApardt[2]+dApardt[1])+0.5*dApardt[0])*dfac_v*q_)/m_; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = -(0.0625*(0.5*(dApardt[1]+dApardt[0])-0.5*(dApardt[3]+dApardt[2]))*dfac_v*q_)/m_; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = -(0.0625*((-0.5*dApardt[3])+0.5*dApardt[2]-0.5*dApardt[1]+0.5*dApardt[0])*dfac_v*q_)/m_; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = -(0.03125*(dApardt[3]+dApardt[2]+dApardt[1]+dApardt[0])*dfac_v*q_)/m_; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = -(0.0625*(0.5*dApardt[3]-0.5*(dApardt[2]+dApardt[1])+0.5*dApardt[0])*dfac_v*q_)/m_; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = -(0.0625*(0.5*(dApardt[1]+dApardt[0])-0.5*(dApardt[3]+dApardt[2]))*dfac_v*q_)/m_; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = -(0.0625*((-0.5*dApardt[3])+0.5*dApardt[2]-0.5*dApardt[1]+0.5*dApardt[0])*dfac_v*q_)/m_; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = -(0.03125*(dApardt[3]+dApardt[2]+dApardt[1]+dApardt[0])*dfac_v*q_)/m_; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // evaluate alpha at right surface quadrature points 
  alphaR = -(0.0625*(0.5*dApardt[3]-0.5*(dApardt[2]+dApardt[1])+0.5*dApardt[0])*dfac_v*q_)/m_; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = -(0.0625*(0.5*(dApardt[1]+dApardt[0])-0.5*(dApardt[3]+dApardt[2]))*dfac_v*q_)/m_; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = -(0.0625*((-0.5*dApardt[3])+0.5*dApardt[2]-0.5*dApardt[1]+0.5*dApardt[0])*dfac_v*q_)/m_; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = -(0.03125*(dApardt[3]+dApardt[2]+dApardt[1]+dApardt[0])*dfac_v*q_)/m_; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = -(0.0625*(0.5*dApardt[3]-0.5*(dApardt[2]+dApardt[1])+0.5*dApardt[0])*dfac_v*q_)/m_; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = -(0.0625*(0.5*(dApardt[1]+dApardt[0])-0.5*(dApardt[3]+dApardt[2]))*dfac_v*q_)/m_; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = -(0.0625*((-0.5*dApardt[3])+0.5*dApardt[2]-0.5*dApardt[1]+0.5*dApardt[0])*dfac_v*q_)/m_; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = -(0.03125*(dApardt[3]+dApardt[2]+dApardt[1]+dApardt[0])*dfac_v*q_)/m_; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#endif 
return cflFreq; 
} 
double EmGyrokineticGenGeoVol2x2vSerP1_Bvars_1(const double q_, const double m_, const double *w, const double *dxv, const double *Bmag, const double *jacobTotInv, const double *Gradpar, const double *geoX, const double *geoY, const double *geoZ, const double *Phi, const double *Apar, const double* dApardt, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_y = 2.0/dxv[1]; 
  double dfac_v = 2.0/dxv[2]; 
  double dfac_m = 2.0/dxv[3]; 
  double wx = w[0]; 
  double wy = w[1]; 
  double wv = w[2]; 
  double wm = w[3]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double m2 = m_*m_; 
  double hamil[16]; 
  hamil[0] = (0.6666666666666666*(3.0*dfac_v2*(m_*wv2+Bmag[0]*wm+Phi[0]*q_)+m_))/dfac_v2; 
  hamil[1] = 2.0*(Bmag[1]*wm+Phi[1]*q_); 
  hamil[2] = 2.0*Phi[2]*q_; 
  hamil[3] = (2.309401076758503*m_*wv)/dfac_v; 
  hamil[4] = (1.154700538379252*Bmag[0])/dfac_m; 
  hamil[5] = 2.0*Phi[3]*q_; 
  hamil[8] = (1.154700538379252*Bmag[1])/dfac_m; 
  double BstarX_by_Bmag[16]; 
  BstarX_by_Bmag[0] = 0.8660254037844386*((geoZ[0]*jacobTotInv[1]+jacobTotInv[0]*geoZ[1])*Apar[3]+(geoZ[1]*jacobTotInv[1]+geoZ[0]*jacobTotInv[0])*Apar[2])*dfac_y; 
  BstarX_by_Bmag[1] = 0.1732050807568877*((9.0*geoZ[1]*jacobTotInv[1]+5.0*geoZ[0]*jacobTotInv[0])*Apar[3]+5.0*(geoZ[0]*jacobTotInv[1]+jacobTotInv[0]*geoZ[1])*Apar[2])*dfac_y; 

  double BstarY_by_Bmag[16]; 
  BstarY_by_Bmag[0] = -(0.8660254037844386*dfac_x*(2.0*jacobTotInv[0]*geoZ[1]*m_*wv+(2.0*Apar[1]*geoZ[1]*jacobTotInv[1]+jacobTotInv[0]*(Apar[0]*geoZ[1]+geoZ[0]*Apar[1]))*q_))/q_; 
  BstarY_by_Bmag[1] = -(0.8660254037844386*dfac_x*(2.0*geoZ[1]*jacobTotInv[1]*m_*wv+((Apar[0]*geoZ[1]+geoZ[0]*Apar[1])*jacobTotInv[1]+2.0*jacobTotInv[0]*Apar[1]*geoZ[1])*q_))/q_; 
  BstarY_by_Bmag[2] = -0.8660254037844386*((2.0*geoZ[1]*jacobTotInv[1]+geoZ[0]*jacobTotInv[0])*Apar[3]+jacobTotInv[0]*geoZ[1]*Apar[2])*dfac_x; 
  BstarY_by_Bmag[3] = -(1.0*jacobTotInv[0]*geoZ[1]*dfac_x*m_)/(dfac_v*q_); 
  BstarY_by_Bmag[5] = -0.8660254037844386*((geoZ[0]*jacobTotInv[1]+2.0*jacobTotInv[0]*geoZ[1])*Apar[3]+geoZ[1]*jacobTotInv[1]*Apar[2])*dfac_x; 
  BstarY_by_Bmag[6] = -(1.0*geoZ[1]*jacobTotInv[1]*dfac_x*m_)/(dfac_v*q_); 

  double cflFreq = 0.0; 
  double alphaL = 0.0; 
  double alphaR = 0.0; 
  double alphax[16]; 
  alphax[0] = 1.732050807568877*dfac_x*((0.25*BstarX_by_Bmag[0]*hamil[3]*dfac_v)/m_-(0.25*((geoZ[0]*jacobTotInv[1]+jacobTotInv[0]*geoZ[1])*hamil[5]+(geoZ[1]*jacobTotInv[1]+geoZ[0]*jacobTotInv[0])*hamil[2])*dfac_y)/q_); 
  alphax[1] = 1.732050807568877*dfac_x*((((-0.25*(geoZ[0]*jacobTotInv[0]*hamil[5]+(geoZ[0]*jacobTotInv[1]+jacobTotInv[0]*geoZ[1])*hamil[2]))-0.45*geoZ[1]*jacobTotInv[1]*hamil[5])*dfac_y)/q_+(0.25*BstarX_by_Bmag[1]*hamil[3]*dfac_v)/m_); 
#if cflType == SURFAVG 
  // evaluate surface-averaged alpha on left 
  alphaL = -0.125*(1.732050807568877*alphax[1]-1.0*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // evaluate surface-averaged alpha on right 
  alphaR = 0.125*(1.732050807568877*alphax[1]+alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#elif cflType == QUAD 
  // evaluate alpha at left surface quadrature points 
  alphaL = 0.0625*(0.25*alphax[0]-0.4330127018922193*alphax[1]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alphax[0]-0.4330127018922193*alphax[1]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alphax[0]-0.4330127018922193*alphax[1]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alphax[0]-0.4330127018922193*alphax[1]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alphax[0]-0.4330127018922193*alphax[1]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alphax[0]-0.4330127018922193*alphax[1]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alphax[0]-0.4330127018922193*alphax[1]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alphax[0]-0.4330127018922193*alphax[1]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // evaluate alpha at right surface quadrature points 
  alphaR = 0.0625*(0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#endif 

  double alphay[16]; 
  alphay[0] = 0.4330127018922193*dfac_y*((hamil[1]*(geoZ[1]*jacobTotInv[1]+geoZ[0]*jacobTotInv[0])*dfac_x)/q_+(BstarY_by_Bmag[0]*hamil[3]*dfac_v)/m_); 
  alphay[1] = 0.4330127018922193*dfac_y*((hamil[1]*(geoZ[0]*jacobTotInv[1]+jacobTotInv[0]*geoZ[1])*dfac_x)/q_+(BstarY_by_Bmag[1]*hamil[3]*dfac_v)/m_); 
  alphay[2] = 0.4330127018922193*dfac_y*(((geoZ[1]*jacobTotInv[1]+geoZ[0]*jacobTotInv[0])*hamil[5]*dfac_x)/q_+(BstarY_by_Bmag[2]*hamil[3]*dfac_v)/m_); 
  alphay[3] = (0.4330127018922193*BstarY_by_Bmag[3]*hamil[3]*dfac_v*dfac_y)/m_; 
  alphay[4] = (0.4330127018922193*(geoZ[1]*jacobTotInv[1]+geoZ[0]*jacobTotInv[0])*hamil[8]*dfac_x*dfac_y)/q_; 
  alphay[5] = 0.4330127018922193*dfac_y*(((geoZ[0]*jacobTotInv[1]+jacobTotInv[0]*geoZ[1])*hamil[5]*dfac_x)/q_+(hamil[3]*BstarY_by_Bmag[5]*dfac_v)/m_); 
  alphay[6] = (0.4330127018922193*hamil[3]*BstarY_by_Bmag[6]*dfac_v*dfac_y)/m_; 
  alphay[8] = (0.4330127018922193*(geoZ[0]*jacobTotInv[1]+jacobTotInv[0]*geoZ[1])*hamil[8]*dfac_x*dfac_y)/q_; 
#if cflType == SURFAVG 
  // evaluate surface-averaged alpha on left 
  alphaL = -0.125*(1.732050807568877*alphay[2]-1.0*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // evaluate surface-averaged alpha on right 
  alphaR = 0.125*(1.732050807568877*alphay[2]+alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#elif cflType == QUAD 
  // evaluate alpha at left surface quadrature points 
  alphaL = 0.0625*(0.25*(alphay[8]+alphay[6])+0.4330127018922193*alphay[5]-0.25*(alphay[4]+alphay[3])-0.4330127018922193*alphay[2]-0.25*alphay[1]+0.25*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.25*(alphay[8]+alphay[6]))-0.4330127018922193*alphay[5]-0.25*(alphay[4]+alphay[3])-0.4330127018922193*alphay[2]+0.25*(alphay[1]+alphay[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alphay[8]-0.25*alphay[6]+0.4330127018922193*alphay[5]-0.25*alphay[4]+0.25*alphay[3]-0.4330127018922193*alphay[2]-0.25*alphay[1]+0.25*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.25*alphay[8])+0.25*alphay[6]-0.4330127018922193*alphay[5]-0.25*alphay[4]+0.25*alphay[3]-0.4330127018922193*alphay[2]+0.25*(alphay[1]+alphay[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.25*alphay[8])+0.25*alphay[6]+0.4330127018922193*alphay[5]+0.25*alphay[4]-0.25*alphay[3]-0.4330127018922193*alphay[2]-0.25*alphay[1]+0.25*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alphay[8]-0.25*alphay[6]-0.4330127018922193*alphay[5]+0.25*alphay[4]-0.25*alphay[3]-0.4330127018922193*alphay[2]+0.25*(alphay[1]+alphay[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.25*(alphay[8]+alphay[6]))+0.4330127018922193*alphay[5]+0.25*(alphay[4]+alphay[3])-0.4330127018922193*alphay[2]-0.25*alphay[1]+0.25*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*(alphay[8]+alphay[6])-0.4330127018922193*alphay[5]+0.25*(alphay[4]+alphay[3])-0.4330127018922193*alphay[2]+0.25*(alphay[1]+alphay[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // evaluate alpha at right surface quadrature points 
  alphaR = 0.0625*(0.25*(alphay[8]+alphay[6])-0.4330127018922193*alphay[5]-0.25*(alphay[4]+alphay[3])+0.4330127018922193*alphay[2]-0.25*alphay[1]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.25*(alphay[8]+alphay[6]))+0.4330127018922193*alphay[5]-0.25*(alphay[4]+alphay[3])+0.4330127018922193*alphay[2]+0.25*(alphay[1]+alphay[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.25*alphay[8]-0.25*alphay[6]-0.4330127018922193*alphay[5]-0.25*alphay[4]+0.25*alphay[3]+0.4330127018922193*alphay[2]-0.25*alphay[1]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.25*alphay[8])+0.25*alphay[6]+0.4330127018922193*alphay[5]-0.25*alphay[4]+0.25*alphay[3]+0.4330127018922193*alphay[2]+0.25*(alphay[1]+alphay[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.25*alphay[8])+0.25*alphay[6]-0.4330127018922193*alphay[5]+0.25*alphay[4]-0.25*alphay[3]+0.4330127018922193*alphay[2]-0.25*alphay[1]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.25*alphay[8]-0.25*alphay[6]+0.4330127018922193*alphay[5]+0.25*alphay[4]-0.25*alphay[3]+0.4330127018922193*alphay[2]+0.25*(alphay[1]+alphay[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.25*(alphay[8]+alphay[6]))-0.4330127018922193*alphay[5]+0.25*(alphay[4]+alphay[3])+0.4330127018922193*alphay[2]-0.25*alphay[1]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.25*(alphay[8]+alphay[6])+0.4330127018922193*alphay[5]+0.25*(alphay[4]+alphay[3])+0.4330127018922193*alphay[2]+0.25*(alphay[1]+alphay[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#endif 

  double alphav[16]; 
  alphav[0] = -(0.4330127018922193*dfac_v*((BstarY_by_Bmag[1]*hamil[5]+BstarY_by_Bmag[0]*hamil[2])*dfac_y+BstarX_by_Bmag[0]*hamil[1]*dfac_x))/m_; 
  alphav[1] = -(0.4330127018922193*dfac_v*((BstarY_by_Bmag[0]*hamil[5]+BstarY_by_Bmag[1]*hamil[2])*dfac_y+BstarX_by_Bmag[1]*hamil[1]*dfac_x))/m_; 
  alphav[2] = -(0.4330127018922193*dfac_v*((BstarY_by_Bmag[5]*hamil[5]+BstarY_by_Bmag[2]*hamil[2])*dfac_y+BstarX_by_Bmag[0]*hamil[5]*dfac_x))/m_; 
  alphav[3] = -(0.4330127018922193*(hamil[5]*BstarY_by_Bmag[6]+hamil[2]*BstarY_by_Bmag[3])*dfac_v*dfac_y)/m_; 
  alphav[4] = -(0.4330127018922193*BstarX_by_Bmag[0]*hamil[8]*dfac_v*dfac_x)/m_; 
  alphav[5] = -(0.4330127018922193*dfac_v*((BstarY_by_Bmag[2]*hamil[5]+hamil[2]*BstarY_by_Bmag[5])*dfac_y+BstarX_by_Bmag[1]*hamil[5]*dfac_x))/m_; 
  alphav[6] = -(0.4330127018922193*(hamil[2]*BstarY_by_Bmag[6]+BstarY_by_Bmag[3]*hamil[5])*dfac_v*dfac_y)/m_; 
  alphav[8] = -(0.4330127018922193*BstarX_by_Bmag[1]*hamil[8]*dfac_v*dfac_x)/m_; 
#if cflType == SURFAVG 
  // evaluate surface-averaged alpha on left 
  alphaL = -0.125*(1.732050807568877*alphav[3]-1.0*alphav[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // evaluate surface-averaged alpha on right 
  alphaR = 0.125*(1.732050807568877*alphav[3]+alphav[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#elif cflType == QUAD 
  // evaluate alpha at left surface quadrature points 
  alphaL = 0.0625*(0.25*alphav[8]+0.4330127018922193*alphav[6]+0.25*alphav[5]-0.25*alphav[4]-0.4330127018922193*alphav[3]-0.25*(alphav[2]+alphav[1])+0.25*alphav[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.25*alphav[8])-0.4330127018922193*alphav[6]-0.25*(alphav[5]+alphav[4])-0.4330127018922193*alphav[3]-0.25*alphav[2]+0.25*(alphav[1]+alphav[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alphav[8]+0.4330127018922193*alphav[6]-0.25*(alphav[5]+alphav[4])-0.4330127018922193*alphav[3]+0.25*alphav[2]-0.25*alphav[1]+0.25*alphav[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.25*alphav[8])-0.4330127018922193*alphav[6]+0.25*alphav[5]-0.25*alphav[4]-0.4330127018922193*alphav[3]+0.25*(alphav[2]+alphav[1]+alphav[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.25*alphav[8])+0.4330127018922193*alphav[6]+0.25*(alphav[5]+alphav[4])-0.4330127018922193*alphav[3]-0.25*(alphav[2]+alphav[1])+0.25*alphav[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alphav[8]-0.4330127018922193*alphav[6]-0.25*alphav[5]+0.25*alphav[4]-0.4330127018922193*alphav[3]-0.25*alphav[2]+0.25*(alphav[1]+alphav[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.25*alphav[8])+0.4330127018922193*alphav[6]-0.25*alphav[5]+0.25*alphav[4]-0.4330127018922193*alphav[3]+0.25*alphav[2]-0.25*alphav[1]+0.25*alphav[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alphav[8]-0.4330127018922193*alphav[6]+0.25*(alphav[5]+alphav[4])-0.4330127018922193*alphav[3]+0.25*(alphav[2]+alphav[1]+alphav[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // evaluate alpha at right surface quadrature points 
  alphaR = 0.0625*(0.25*alphav[8]-0.4330127018922193*alphav[6]+0.25*alphav[5]-0.25*alphav[4]+0.4330127018922193*alphav[3]-0.25*(alphav[2]+alphav[1])+0.25*alphav[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.25*alphav[8])+0.4330127018922193*alphav[6]-0.25*(alphav[5]+alphav[4])+0.4330127018922193*alphav[3]-0.25*alphav[2]+0.25*(alphav[1]+alphav[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.25*alphav[8]-0.4330127018922193*alphav[6]-0.25*(alphav[5]+alphav[4])+0.4330127018922193*alphav[3]+0.25*alphav[2]-0.25*alphav[1]+0.25*alphav[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.25*alphav[8])+0.4330127018922193*alphav[6]+0.25*alphav[5]-0.25*alphav[4]+0.4330127018922193*alphav[3]+0.25*(alphav[2]+alphav[1]+alphav[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.25*alphav[8])-0.4330127018922193*alphav[6]+0.25*(alphav[5]+alphav[4])+0.4330127018922193*alphav[3]-0.25*(alphav[2]+alphav[1])+0.25*alphav[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.25*alphav[8]+0.4330127018922193*alphav[6]-0.25*alphav[5]+0.25*alphav[4]+0.4330127018922193*alphav[3]-0.25*alphav[2]+0.25*(alphav[1]+alphav[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.25*alphav[8])-0.4330127018922193*alphav[6]-0.25*alphav[5]+0.25*alphav[4]+0.4330127018922193*alphav[3]+0.25*alphav[2]-0.25*alphav[1]+0.25*alphav[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.25*alphav[8]+0.4330127018922193*alphav[6]+0.25*(alphav[5]+alphav[4])+0.4330127018922193*alphav[3]+0.25*(alphav[2]+alphav[1]+alphav[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#endif 

  out[1] += 0.4330127018922193*(alphax[1]*f[1]+alphax[0]*f[0]); 
  out[2] += 0.4330127018922193*(alphay[8]*f[8]+alphay[6]*f[6]+alphay[5]*f[5]+alphay[4]*f[4]+alphay[3]*f[3]+alphay[2]*f[2]+alphay[1]*f[1]+alphay[0]*f[0]); 
  out[3] += 0.4330127018922193*(alphav[8]*f[8]+alphav[6]*f[6]+alphav[5]*f[5]+alphav[4]*f[4]+alphav[3]*f[3]+alphav[2]*f[2]+alphav[1]*f[1]+alphav[0]*f[0]); 
  out[5] += 0.4330127018922193*(alphay[4]*f[8]+f[4]*alphay[8]+alphay[3]*f[6]+f[3]*alphay[6]+(alphay[2]+alphax[1])*f[5]+f[2]*(alphay[5]+alphax[0])+alphay[0]*f[1]+f[0]*alphay[1]); 
  out[6] += 0.4330127018922193*(alphav[4]*f[8]+f[4]*alphav[8]+(alphav[3]+alphax[1])*f[6]+f[3]*alphav[6]+alphav[2]*f[5]+f[2]*alphav[5]+alphax[0]*f[3]+alphav[0]*f[1]+f[0]*alphav[1]); 
  out[7] += 0.4330127018922193*(alphay[8]*f[13]+alphav[8]*f[12]+(alphav[6]+alphay[5])*f[11]+alphay[4]*f[10]+alphav[4]*f[9]+(alphav[3]+alphay[2])*f[7]+alphay[1]*f[6]+f[1]*alphay[6]+alphav[1]*f[5]+f[1]*alphav[5]+alphay[0]*f[3]+f[0]*alphay[3]+alphav[0]*f[2]+f[0]*alphav[2]); 
  out[8] += 0.4330127018922193*(alphax[1]*f[8]+alphax[0]*f[4]); 
  out[9] += 0.4330127018922193*(alphay[6]*f[13]+alphay[5]*f[12]+alphay[3]*f[10]+alphay[2]*f[9]+alphay[1]*f[8]+f[1]*alphay[8]+alphay[0]*f[4]+f[0]*alphay[4]); 
  out[10] += 0.4330127018922193*(alphav[6]*f[13]+alphav[5]*f[12]+alphav[3]*f[10]+alphav[2]*f[9]+alphav[1]*f[8]+f[1]*alphav[8]+alphav[0]*f[4]+f[0]*alphav[4]); 
  out[11] += 0.4330127018922193*(alphay[4]*f[13]+alphav[4]*f[12]+(alphav[3]+alphay[2]+alphax[1])*f[11]+alphay[8]*f[10]+alphav[8]*f[9]+(alphav[6]+alphay[5]+alphax[0])*f[7]+alphay[0]*f[6]+f[0]*alphay[6]+alphav[0]*f[5]+f[0]*alphav[5]+alphay[1]*f[3]+f[1]*alphay[3]+alphav[1]*f[2]+f[1]*alphav[2]); 
  out[12] += 0.4330127018922193*(alphay[3]*f[13]+(alphay[2]+alphax[1])*f[12]+alphay[6]*f[10]+(alphay[5]+alphax[0])*f[9]+alphay[0]*f[8]+f[0]*alphay[8]+alphay[1]*f[4]+f[1]*alphay[4]); 
  out[13] += 0.4330127018922193*((alphav[3]+alphax[1])*f[13]+alphav[2]*f[12]+(alphav[6]+alphax[0])*f[10]+alphav[5]*f[9]+alphav[0]*f[8]+f[0]*alphav[8]+alphav[1]*f[4]+f[1]*alphav[4]); 
  out[14] += 0.4330127018922193*((alphav[6]+alphay[5])*f[15]+(alphav[3]+alphay[2])*f[14]+alphay[1]*f[13]+alphav[1]*f[12]+alphay[0]*f[10]+alphav[0]*f[9]+(alphay[6]+alphav[5])*f[8]+f[6]*alphay[8]+f[5]*alphav[8]+(alphay[3]+alphav[2])*f[4]+f[3]*alphay[4]+f[2]*alphav[4]); 
  out[15] += 0.4330127018922193*((alphav[3]+alphay[2]+alphax[1])*f[15]+(alphav[6]+alphay[5]+alphax[0])*f[14]+alphay[0]*f[13]+alphav[0]*f[12]+alphay[1]*f[10]+alphav[1]*f[9]+(alphay[3]+alphav[2])*f[8]+f[3]*alphay[8]+f[2]*alphav[8]+alphay[4]*f[6]+f[4]*alphay[6]+alphav[4]*f[5]+f[4]*alphav[5]); 
  return cflFreq; 
} 
