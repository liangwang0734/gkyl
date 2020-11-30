#include <PassiveAdvectionModDecl.h> 
double PassiveAdvectionVol2xSerP2(const double *w, const double *dxv, double *positivityWeightByDir, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing.
  double dfac1 = 2.0/dxv[0]; 
  double w1 = w[0]; 
  double dfac2 = 2.0/dxv[1]; 
  double w2 = w[1]; 
  const double *v1 = &f[8]; 
  const double *v2 = &f[16]; 
  double cflRate = 0.0; 
  double alphaL = 0.0; 
  double alphaR = 0.0; 
  positivityWeightByDir[1] = 0.; 
#if cflType == SURFAVG 
  // evaluate surface-averaged alpha on left 
  alphaL = 0.25*(2.23606797749979*v1[4]-1.732050807568877*v1[1]+v1[0])*dfac1; 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  positivityWeightByDir[1] += -0.5*(alphaL-std::abs(alphaL)); 
  // evaluate surface-averaged alpha on right 
  alphaR = 0.25*(2.23606797749979*v1[4]+1.732050807568877*v1[1]+v1[0])*dfac1; 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  positivityWeightByDir[1] += 0.5*(alphaR+std::abs(alphaR)); 
#elif cflType == QUAD 
  // evaluate alpha at left surface quadrature points 
  alphaL = 0.25*((-1.936491673103709*(v1[7]+v1[6]))+1.118033988749895*(v1[5]+v1[4])+1.5*v1[3]-0.8660254037844386*(v1[2]+v1[1])+0.5*v1[0])*dfac1; 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  positivityWeightByDir[1] += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.25*(0.9682458365518543*v1[7]-0.5590169943749475*v1[5]+1.118033988749895*v1[4]-0.8660254037844386*v1[1]+0.5*v1[0])*dfac1; 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  positivityWeightByDir[1] += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.25*((-1.936491673103709*v1[7])+1.936491673103709*v1[6]+1.118033988749895*(v1[5]+v1[4])-1.5*v1[3]+0.8660254037844386*v1[2]-0.8660254037844386*v1[1]+0.5*v1[0])*dfac1; 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  positivityWeightByDir[1] += -0.5*(alphaL-std::abs(alphaL)); 
  // evaluate alpha at right surface quadrature points 
  alphaR = 0.25*(1.936491673103709*v1[7]-1.936491673103709*v1[6]+1.118033988749895*(v1[5]+v1[4])-1.5*v1[3]-0.8660254037844386*v1[2]+0.8660254037844386*v1[1]+0.5*v1[0])*dfac1; 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  positivityWeightByDir[1] += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.25*((-0.9682458365518543*v1[7])-0.5590169943749475*v1[5]+1.118033988749895*v1[4]+0.8660254037844386*v1[1]+0.5*v1[0])*dfac1; 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  positivityWeightByDir[1] += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.25*(1.936491673103709*(v1[7]+v1[6])+1.118033988749895*(v1[5]+v1[4])+1.5*v1[3]+0.8660254037844386*(v1[2]+v1[1])+0.5*v1[0])*dfac1; 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  positivityWeightByDir[1] += 0.5*(alphaR+std::abs(alphaR)); 
#endif 

  positivityWeightByDir[2] = 0.; 
#if cflType == SURFAVG 
  // evaluate surface-averaged alpha on left 
  alphaL = 0.25*(2.23606797749979*v2[5]-1.732050807568877*v2[2]+v2[0])*dfac2; 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  positivityWeightByDir[2] += -0.5*(alphaL-std::abs(alphaL)); 
  // evaluate surface-averaged alpha on right 
  alphaR = 0.25*(2.23606797749979*v2[5]+1.732050807568877*v2[2]+v2[0])*dfac2; 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  positivityWeightByDir[2] += 0.5*(alphaR+std::abs(alphaR)); 
#elif cflType == QUAD 
  // evaluate alpha at left surface quadrature points 
  alphaL = 0.25*((-1.936491673103709*(v2[7]+v2[6]))+1.118033988749895*(v2[5]+v2[4])+1.5*v2[3]-0.8660254037844386*(v2[2]+v2[1])+0.5*v2[0])*dfac2; 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  positivityWeightByDir[2] += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.25*(0.9682458365518543*v2[6]+1.118033988749895*v2[5]-0.5590169943749475*v2[4]-0.8660254037844386*v2[2]+0.5*v2[0])*dfac2; 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  positivityWeightByDir[2] += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.25*(1.936491673103709*v2[7]-1.936491673103709*v2[6]+1.118033988749895*(v2[5]+v2[4])-1.5*v2[3]-0.8660254037844386*v2[2]+0.8660254037844386*v2[1]+0.5*v2[0])*dfac2; 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  positivityWeightByDir[2] += -0.5*(alphaL-std::abs(alphaL)); 
  // evaluate alpha at right surface quadrature points 
  alphaR = 0.25*((-1.936491673103709*v2[7])+1.936491673103709*v2[6]+1.118033988749895*(v2[5]+v2[4])-1.5*v2[3]+0.8660254037844386*v2[2]-0.8660254037844386*v2[1]+0.5*v2[0])*dfac2; 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  positivityWeightByDir[2] += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.25*((-0.9682458365518543*v2[6])+1.118033988749895*v2[5]-0.5590169943749475*v2[4]+0.8660254037844386*v2[2]+0.5*v2[0])*dfac2; 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  positivityWeightByDir[2] += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.25*(1.936491673103709*(v2[7]+v2[6])+1.118033988749895*(v2[5]+v2[4])+1.5*v2[3]+0.8660254037844386*(v2[2]+v2[1])+0.5*v2[0])*dfac2; 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  positivityWeightByDir[2] += 0.5*(alphaR+std::abs(alphaR)); 
#endif 

  out[1] += 0.8660254037844386*(f[7]*v1[7]+f[6]*v1[6]+f[5]*v1[5]+f[4]*v1[4]+f[3]*v1[3]+f[2]*v1[2]+f[1]*v1[1]+f[0]*v1[0])*dfac1; 
  out[2] += 0.8660254037844386*(f[7]*v2[7]+f[6]*v2[6]+f[5]*v2[5]+f[4]*v2[4]+f[3]*v2[3]+f[2]*v2[2]+f[1]*v2[1]+f[0]*v2[0])*dfac2; 
  out[3] += 0.1*((8.660254037844387*(f[5]*v2[7]+v2[5]*f[7])+7.745966692414834*(f[3]*v2[6]+v2[3]*f[6]+f[1]*v2[4]+v2[1]*f[4])+8.660254037844386*(f[2]*v2[3]+v2[2]*f[3]+f[0]*v2[1]+v2[0]*f[1]))*dfac2+(7.745966692414834*(f[3]*v1[7]+v1[3]*f[7])+8.660254037844387*(f[4]*v1[6]+v1[4]*f[6])+7.745966692414834*(f[2]*v1[5]+v1[2]*f[5])+8.660254037844386*(f[1]*v1[3]+v1[1]*f[3]+f[0]*v1[2]+v1[0]*f[2]))*dfac1); 
  out[4] += 0.1*(19.36491673103708*(f[5]*v1[7]+v1[5]*f[7])+17.32050807568877*(f[3]*v1[6]+v1[3]*f[6])+17.32050807568877*(f[1]*v1[4]+v1[1]*f[4])+19.36491673103709*(f[2]*v1[3]+v1[2]*f[3]+f[0]*v1[1]+v1[0]*f[1]))*dfac1; 
  out[5] += 0.1*(17.32050807568877*(f[3]*v2[7]+v2[3]*f[7])+19.36491673103708*(f[4]*v2[6]+v2[4]*f[6])+17.32050807568877*(f[2]*v2[5]+v2[2]*f[5])+19.36491673103709*(f[1]*v2[3]+v2[1]*f[3]+f[0]*v2[2]+v2[0]*f[2]))*dfac2; 
  out[6] += 0.01428571428571429*((54.22176684690384*f[7]*v2[7]+(38.72983346207417*f[6]+60.6217782649107*f[2])*v2[6]+60.6217782649107*v2[2]*f[6]+(38.72983346207417*f[4]+60.62177826491071*f[0])*v2[4]+60.62177826491071*v2[0]*f[4]+54.22176684690384*(f[3]*v2[3]+f[1]*v2[1]))*dfac2+((108.4435336938077*f[6]+121.2435565298214*f[2])*v1[7]+(108.4435336938077*v1[6]+121.2435565298214*v1[2])*f[7]+121.2435565298214*(f[1]*v1[6]+v1[1]*f[6])+121.2435565298214*(f[3]*v1[5]+v1[3]*f[5]+f[3]*v1[4]+v1[3]*f[4])+135.5544171172596*(f[0]*v1[3]+v1[0]*f[3]+f[1]*v1[2]+v1[1]*f[2]))*dfac1); 
  out[7] += 0.01428571428571429*(((108.4435336938077*f[6]+121.2435565298214*f[2])*v2[7]+(108.4435336938077*v2[6]+121.2435565298214*v2[2])*f[7]+121.2435565298214*(f[1]*v2[6]+v2[1]*f[6])+121.2435565298214*(f[3]*v2[5]+v2[3]*f[5]+f[3]*v2[4]+v2[3]*f[4])+135.5544171172596*(f[0]*v2[3]+v2[0]*f[3]+f[1]*v2[2]+v2[1]*f[2]))*dfac2+((38.72983346207417*f[7]+60.6217782649107*f[1])*v1[7]+60.6217782649107*v1[1]*f[7]+54.22176684690384*f[6]*v1[6]+(38.72983346207417*f[5]+60.62177826491071*f[0])*v1[5]+60.62177826491071*v1[0]*f[5]+54.22176684690384*(f[3]*v1[3]+f[2]*v1[2]))*dfac1); 
  return cflRate; 
} 