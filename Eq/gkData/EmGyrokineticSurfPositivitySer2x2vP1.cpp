#include <GyrokineticModDecl.h> 
double EmGyrokineticSurfPositivity2x2vSer_X_P1_Bvars_0(const double q_, const double m_, const double cflL, const double cflR, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *Apar, const double *dApardt, const double *fl, const double *fr, double *outl, double *outr) 
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
  double incr[16]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = (0.125*(2.0*BdriftX[0]*m_*wv2+q_*((BmagInv[0]*(1.732050807568877*Apar[2]-3.0*Apar[3])*dfac_y+BdriftX[0]*(Apar[0]-1.732050807568877*Apar[1]))*wv+BmagInv[0]*(3.0*Phi[3]-1.732050807568877*Phi[2])*dfac_y)))/q_; 

  double alpha[8]; 
  alpha[0] = (0.7071067811865475*(2.0*BdriftX[0]*m_*wv2+q_*((BmagInv[0]*(1.732050807568877*Apar[2]-3.0*Apar[3])*dfac_y+BdriftX[0]*(Apar[0]-1.732050807568877*Apar[1]))*wv+BmagInv[0]*(3.0*Phi[3]-1.732050807568877*Phi[2])*dfac_y)))/q_; 
  alpha[1] = -0.408248290463863*BdriftX[0]*(3.0*Apar[3]-1.732050807568877*Apar[2])*wv; 
  alpha[2] = (0.8164965809277261*BdriftX[0]*m_*wv)/(dfac_v*q_); 
  double rVal;  // rVal=f1/f0 at each control node in dimensions other than x 
  double fqVal[8];  // fqVal = anti-limited f evaluated at each control node on x surface 
  // determine upwinding at each surface control node 
  if(0.3535533905932737*alpha[0]-0.2041241452319315*(alpha[2]+alpha[1]) > 0) {
  rVal = -(1.0*(3.0*fl[15]-5.196152422706631*(fl[13]+fl[12]+fl[11])+9.0*(fl[8]+fl[6]+fl[5])-15.58845726811989*fl[1]))/(36.0*EPSILON-1.732050807568877*fl[14]+3.0*(fl[10]+fl[9]+fl[7])-5.196152422706631*(fl[4]+fl[3]+fl[2])+9.0*fl[0]); 
  fqVal[0] = -0.02777777777777778*(1.732050807568877*fl[14]-3.0*(fl[10]+fl[9]+fl[7])+5.196152422706631*(fl[4]+fl[3]+fl[2])-9.0*fl[0])*limTheta(rVal,1.0,cflL); 
  } else {
  rVal = -(1.0*(3.0*fr[15]-5.196152422706631*(fr[13]+fr[12]+fr[11])+9.0*(fr[8]+fr[6]+fr[5])-15.58845726811989*fr[1]))/(36.0*EPSILON-1.732050807568877*fr[14]+3.0*(fr[10]+fr[9]+fr[7])-5.196152422706631*(fr[4]+fr[3]+fr[2])+9.0*fr[0]); 
  fqVal[0] = -0.02777777777777778*(1.732050807568877*fr[14]-3.0*(fr[10]+fr[9]+fr[7])+5.196152422706631*(fr[4]+fr[3]+fr[2])-9.0*fr[0])*limTheta(rVal,-1.0,cflR); 
  }
  if((-0.2041241452319315*alpha[2])+0.2041241452319315*alpha[1]+0.3535533905932737*alpha[0] > 0) {
  rVal = (3.0*fl[15]+5.196152422706631*fl[13]-5.196152422706631*(fl[12]+fl[11])-9.0*(fl[8]+fl[6])+9.0*fl[5]+15.58845726811989*fl[1])/(36.0*EPSILON+1.732050807568877*fl[14]+3.0*fl[10]-3.0*(fl[9]+fl[7])-5.196152422706631*(fl[4]+fl[3])+5.196152422706631*fl[2]+9.0*fl[0]); 
  fqVal[1] = 0.02777777777777778*(1.732050807568877*fl[14]+3.0*fl[10]-3.0*(fl[9]+fl[7])-5.196152422706631*(fl[4]+fl[3])+5.196152422706631*fl[2]+9.0*fl[0])*limTheta(rVal,1.0,cflL); 
  } else {
  rVal = (3.0*fr[15]+5.196152422706631*fr[13]-5.196152422706631*(fr[12]+fr[11])-9.0*(fr[8]+fr[6])+9.0*fr[5]+15.58845726811989*fr[1])/(36.0*EPSILON+1.732050807568877*fr[14]+3.0*fr[10]-3.0*(fr[9]+fr[7])-5.196152422706631*(fr[4]+fr[3])+5.196152422706631*fr[2]+9.0*fr[0]); 
  fqVal[1] = 0.02777777777777778*(1.732050807568877*fr[14]+3.0*fr[10]-3.0*(fr[9]+fr[7])-5.196152422706631*(fr[4]+fr[3])+5.196152422706631*fr[2]+9.0*fr[0])*limTheta(rVal,-1.0,cflR); 
  }
  if(0.2041241452319315*alpha[2]-0.2041241452319315*alpha[1]+0.3535533905932737*alpha[0] > 0) {
  rVal = (3.0*fl[15]-5.196152422706631*fl[13]+5.196152422706631*fl[12]-5.196152422706631*fl[11]-9.0*fl[8]+9.0*fl[6]-9.0*fl[5]+15.58845726811989*fl[1])/(36.0*EPSILON+1.732050807568877*fl[14]-3.0*fl[10]+3.0*fl[9]-3.0*fl[7]-5.196152422706631*fl[4]+5.196152422706631*fl[3]-5.196152422706631*fl[2]+9.0*fl[0]); 
  fqVal[2] = 0.02777777777777778*(1.732050807568877*fl[14]-3.0*fl[10]+3.0*fl[9]-3.0*fl[7]-5.196152422706631*fl[4]+5.196152422706631*fl[3]-5.196152422706631*fl[2]+9.0*fl[0])*limTheta(rVal,1.0,cflL); 
  } else {
  rVal = (3.0*fr[15]-5.196152422706631*fr[13]+5.196152422706631*fr[12]-5.196152422706631*fr[11]-9.0*fr[8]+9.0*fr[6]-9.0*fr[5]+15.58845726811989*fr[1])/(36.0*EPSILON+1.732050807568877*fr[14]-3.0*fr[10]+3.0*fr[9]-3.0*fr[7]-5.196152422706631*fr[4]+5.196152422706631*fr[3]-5.196152422706631*fr[2]+9.0*fr[0]); 
  fqVal[2] = 0.02777777777777778*(1.732050807568877*fr[14]-3.0*fr[10]+3.0*fr[9]-3.0*fr[7]-5.196152422706631*fr[4]+5.196152422706631*fr[3]-5.196152422706631*fr[2]+9.0*fr[0])*limTheta(rVal,-1.0,cflR); 
  }
  if(0.2041241452319315*(alpha[2]+alpha[1])+0.3535533905932737*alpha[0] > 0) {
  rVal = -(1.0*(3.0*fl[15]+5.196152422706631*(fl[13]+fl[12])-5.196152422706631*fl[11]+9.0*fl[8]-9.0*(fl[6]+fl[5])-15.58845726811989*fl[1]))/(36.0*EPSILON-1.732050807568877*fl[14]-3.0*(fl[10]+fl[9])+3.0*fl[7]-5.196152422706631*fl[4]+5.196152422706631*(fl[3]+fl[2])+9.0*fl[0]); 
  fqVal[3] = -0.02777777777777778*(1.732050807568877*fl[14]+3.0*(fl[10]+fl[9])-3.0*fl[7]+5.196152422706631*fl[4]-5.196152422706631*(fl[3]+fl[2])-9.0*fl[0])*limTheta(rVal,1.0,cflL); 
  } else {
  rVal = -(1.0*(3.0*fr[15]+5.196152422706631*(fr[13]+fr[12])-5.196152422706631*fr[11]+9.0*fr[8]-9.0*(fr[6]+fr[5])-15.58845726811989*fr[1]))/(36.0*EPSILON-1.732050807568877*fr[14]-3.0*(fr[10]+fr[9])+3.0*fr[7]-5.196152422706631*fr[4]+5.196152422706631*(fr[3]+fr[2])+9.0*fr[0]); 
  fqVal[3] = -0.02777777777777778*(1.732050807568877*fr[14]+3.0*(fr[10]+fr[9])-3.0*fr[7]+5.196152422706631*fr[4]-5.196152422706631*(fr[3]+fr[2])-9.0*fr[0])*limTheta(rVal,-1.0,cflR); 
  }
  if(0.3535533905932737*alpha[0]-0.2041241452319315*(alpha[2]+alpha[1]) > 0) {
  rVal = (3.0*fl[15]-5.196152422706631*(fl[13]+fl[12])+5.196152422706631*fl[11]+9.0*fl[8]-9.0*(fl[6]+fl[5])+15.58845726811989*fl[1])/(36.0*EPSILON+1.732050807568877*fl[14]-3.0*(fl[10]+fl[9])+3.0*fl[7]+5.196152422706631*fl[4]-5.196152422706631*(fl[3]+fl[2])+9.0*fl[0]); 
  fqVal[4] = 0.02777777777777778*(1.732050807568877*fl[14]-3.0*(fl[10]+fl[9])+3.0*fl[7]+5.196152422706631*fl[4]-5.196152422706631*(fl[3]+fl[2])+9.0*fl[0])*limTheta(rVal,1.0,cflL); 
  } else {
  rVal = (3.0*fr[15]-5.196152422706631*(fr[13]+fr[12])+5.196152422706631*fr[11]+9.0*fr[8]-9.0*(fr[6]+fr[5])+15.58845726811989*fr[1])/(36.0*EPSILON+1.732050807568877*fr[14]-3.0*(fr[10]+fr[9])+3.0*fr[7]+5.196152422706631*fr[4]-5.196152422706631*(fr[3]+fr[2])+9.0*fr[0]); 
  fqVal[4] = 0.02777777777777778*(1.732050807568877*fr[14]-3.0*(fr[10]+fr[9])+3.0*fr[7]+5.196152422706631*fr[4]-5.196152422706631*(fr[3]+fr[2])+9.0*fr[0])*limTheta(rVal,-1.0,cflR); 
  }
  if((-0.2041241452319315*alpha[2])+0.2041241452319315*alpha[1]+0.3535533905932737*alpha[0] > 0) {
  rVal = -(1.0*(3.0*fl[15]+5.196152422706631*fl[13]-5.196152422706631*fl[12]+5.196152422706631*fl[11]-9.0*fl[8]+9.0*fl[6]-9.0*fl[5]-15.58845726811989*fl[1]))/(36.0*EPSILON-1.732050807568877*fl[14]-3.0*fl[10]+3.0*fl[9]-3.0*fl[7]+5.196152422706631*fl[4]-5.196152422706631*fl[3]+5.196152422706631*fl[2]+9.0*fl[0]); 
  fqVal[5] = -0.02777777777777778*(1.732050807568877*fl[14]+3.0*fl[10]-3.0*fl[9]+3.0*fl[7]-5.196152422706631*fl[4]+5.196152422706631*fl[3]-5.196152422706631*fl[2]-9.0*fl[0])*limTheta(rVal,1.0,cflL); 
  } else {
  rVal = -(1.0*(3.0*fr[15]+5.196152422706631*fr[13]-5.196152422706631*fr[12]+5.196152422706631*fr[11]-9.0*fr[8]+9.0*fr[6]-9.0*fr[5]-15.58845726811989*fr[1]))/(36.0*EPSILON-1.732050807568877*fr[14]-3.0*fr[10]+3.0*fr[9]-3.0*fr[7]+5.196152422706631*fr[4]-5.196152422706631*fr[3]+5.196152422706631*fr[2]+9.0*fr[0]); 
  fqVal[5] = -0.02777777777777778*(1.732050807568877*fr[14]+3.0*fr[10]-3.0*fr[9]+3.0*fr[7]-5.196152422706631*fr[4]+5.196152422706631*fr[3]-5.196152422706631*fr[2]-9.0*fr[0])*limTheta(rVal,-1.0,cflR); 
  }
  if(0.2041241452319315*alpha[2]-0.2041241452319315*alpha[1]+0.3535533905932737*alpha[0] > 0) {
  rVal = -(1.0*(3.0*fl[15]-5.196152422706631*fl[13]+5.196152422706631*(fl[12]+fl[11])-9.0*(fl[8]+fl[6])+9.0*fl[5]-15.58845726811989*fl[1]))/(36.0*EPSILON-1.732050807568877*fl[14]+3.0*fl[10]-3.0*(fl[9]+fl[7])+5.196152422706631*(fl[4]+fl[3])-5.196152422706631*fl[2]+9.0*fl[0]); 
  fqVal[6] = -0.02777777777777778*(1.732050807568877*fl[14]-3.0*fl[10]+3.0*(fl[9]+fl[7])-5.196152422706631*(fl[4]+fl[3])+5.196152422706631*fl[2]-9.0*fl[0])*limTheta(rVal,1.0,cflL); 
  } else {
  rVal = -(1.0*(3.0*fr[15]-5.196152422706631*fr[13]+5.196152422706631*(fr[12]+fr[11])-9.0*(fr[8]+fr[6])+9.0*fr[5]-15.58845726811989*fr[1]))/(36.0*EPSILON-1.732050807568877*fr[14]+3.0*fr[10]-3.0*(fr[9]+fr[7])+5.196152422706631*(fr[4]+fr[3])-5.196152422706631*fr[2]+9.0*fr[0]); 
  fqVal[6] = -0.02777777777777778*(1.732050807568877*fr[14]-3.0*fr[10]+3.0*(fr[9]+fr[7])-5.196152422706631*(fr[4]+fr[3])+5.196152422706631*fr[2]-9.0*fr[0])*limTheta(rVal,-1.0,cflR); 
  }
  if(0.2041241452319315*(alpha[2]+alpha[1])+0.3535533905932737*alpha[0] > 0) {
  rVal = (3.0*fl[15]+5.196152422706631*(fl[13]+fl[12]+fl[11])+9.0*(fl[8]+fl[6]+fl[5])+15.58845726811989*fl[1])/(36.0*EPSILON+1.732050807568877*fl[14]+3.0*(fl[10]+fl[9]+fl[7])+5.196152422706631*(fl[4]+fl[3]+fl[2])+9.0*fl[0]); 
  fqVal[7] = 0.02777777777777778*(1.732050807568877*fl[14]+3.0*(fl[10]+fl[9]+fl[7])+5.196152422706631*(fl[4]+fl[3]+fl[2])+9.0*fl[0])*limTheta(rVal,1.0,cflL); 
  } else {
  rVal = (3.0*fr[15]+5.196152422706631*(fr[13]+fr[12]+fr[11])+9.0*(fr[8]+fr[6]+fr[5])+15.58845726811989*fr[1])/(36.0*EPSILON+1.732050807568877*fr[14]+3.0*(fr[10]+fr[9]+fr[7])+5.196152422706631*(fr[4]+fr[3]+fr[2])+9.0*fr[0]); 
  fqVal[7] = 0.02777777777777778*(1.732050807568877*fr[14]+3.0*(fr[10]+fr[9]+fr[7])+5.196152422706631*(fr[4]+fr[3]+fr[2])+9.0*fr[0])*limTheta(rVal,-1.0,cflR); 
  }
  double fhatALVal[8];  // fhatALVal = mode coefficients of anti-limited f 
  fhatALVal[0] = 0.3535533905932737*(fqVal[7]+fqVal[6]+fqVal[5]+fqVal[4]+fqVal[3]+fqVal[2]+fqVal[1]+fqVal[0]); 
  fhatALVal[1] = 0.6123724356957944*(fqVal[7]-1.0*fqVal[6]+fqVal[5]-1.0*fqVal[4]+fqVal[3]-1.0*fqVal[2]+fqVal[1]-1.0*fqVal[0]); 
  fhatALVal[2] = 0.6123724356957944*(fqVal[7]+fqVal[6]-1.0*(fqVal[5]+fqVal[4])+fqVal[3]+fqVal[2]-1.0*(fqVal[1]+fqVal[0])); 
  fhatALVal[3] = 0.6123724356957944*(fqVal[7]+fqVal[6]+fqVal[5]+fqVal[4]-1.0*(fqVal[3]+fqVal[2]+fqVal[1]+fqVal[0])); 
  fhatALVal[4] = 1.060660171779821*(fqVal[7]-1.0*(fqVal[6]+fqVal[5])+fqVal[4]+fqVal[3]-1.0*(fqVal[2]+fqVal[1])+fqVal[0]); 
  fhatALVal[5] = 1.060660171779821*(fqVal[7]-1.0*fqVal[6]+fqVal[5]-1.0*(fqVal[4]+fqVal[3])+fqVal[2]-1.0*fqVal[1]+fqVal[0]); 
  fhatALVal[6] = 1.060660171779821*(fqVal[7]+fqVal[6]-1.0*(fqVal[5]+fqVal[4]+fqVal[3]+fqVal[2])+fqVal[1]+fqVal[0]); 
  fhatALVal[7] = 1.837117307087383*(fqVal[7]-1.0*(fqVal[6]+fqVal[5])+fqVal[4]-1.0*fqVal[3]+fqVal[2]+fqVal[1]-1.0*fqVal[0]); 
  incr[0] = 0.25*(alpha[2]*fhatALVal[2]+alpha[1]*fhatALVal[1]+alpha[0]*fhatALVal[0])*dfac_x; 
  incr[1] = -0.4330127018922193*(alpha[2]*fhatALVal[2]+alpha[1]*fhatALVal[1]+alpha[0]*fhatALVal[0])*dfac_x; 
  incr[2] = 0.25*(alpha[2]*fhatALVal[4]+alpha[0]*fhatALVal[1]+fhatALVal[0]*alpha[1])*dfac_x; 
  incr[3] = 0.25*(alpha[1]*fhatALVal[4]+alpha[0]*fhatALVal[2]+fhatALVal[0]*alpha[2])*dfac_x; 
  incr[4] = 0.25*(alpha[2]*fhatALVal[6]+alpha[1]*fhatALVal[5]+alpha[0]*fhatALVal[3])*dfac_x; 
  incr[5] = -0.4330127018922193*(alpha[2]*fhatALVal[4]+alpha[0]*fhatALVal[1]+fhatALVal[0]*alpha[1])*dfac_x; 
  incr[6] = -0.4330127018922193*(alpha[1]*fhatALVal[4]+alpha[0]*fhatALVal[2]+fhatALVal[0]*alpha[2])*dfac_x; 
  incr[7] = 0.25*(alpha[0]*fhatALVal[4]+alpha[1]*fhatALVal[2]+fhatALVal[1]*alpha[2])*dfac_x; 
  incr[8] = -0.4330127018922193*(alpha[2]*fhatALVal[6]+alpha[1]*fhatALVal[5]+alpha[0]*fhatALVal[3])*dfac_x; 
  incr[9] = 0.25*(alpha[2]*fhatALVal[7]+alpha[0]*fhatALVal[5]+alpha[1]*fhatALVal[3])*dfac_x; 
  incr[10] = 0.25*(alpha[1]*fhatALVal[7]+alpha[0]*fhatALVal[6]+alpha[2]*fhatALVal[3])*dfac_x; 
  incr[11] = -0.4330127018922193*(alpha[0]*fhatALVal[4]+alpha[1]*fhatALVal[2]+fhatALVal[1]*alpha[2])*dfac_x; 
  incr[12] = -0.4330127018922193*(alpha[2]*fhatALVal[7]+alpha[0]*fhatALVal[5]+alpha[1]*fhatALVal[3])*dfac_x; 
  incr[13] = -0.4330127018922193*(alpha[1]*fhatALVal[7]+alpha[0]*fhatALVal[6]+alpha[2]*fhatALVal[3])*dfac_x; 
  incr[14] = 0.25*(alpha[0]*fhatALVal[7]+alpha[1]*fhatALVal[6]+alpha[2]*fhatALVal[5])*dfac_x; 
  incr[15] = -0.4330127018922193*(alpha[0]*fhatALVal[7]+alpha[1]*fhatALVal[6]+alpha[2]*fhatALVal[5])*dfac_x; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 
  outr[4] += incr[4]; 
  outr[5] += incr[5]; 
  outr[6] += incr[6]; 
  outr[7] += incr[7]; 
  outr[8] += incr[8]; 
  outr[9] += incr[9]; 
  outr[10] += incr[10]; 
  outr[11] += incr[11]; 
  outr[12] += incr[12]; 
  outr[13] += incr[13]; 
  outr[14] += incr[14]; 
  outr[15] += incr[15]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += incr[1]; 
  outl[2] += -1.0*incr[2]; 
  outl[3] += -1.0*incr[3]; 
  outl[4] += -1.0*incr[4]; 
  outl[5] += incr[5]; 
  outl[6] += incr[6]; 
  outl[7] += -1.0*incr[7]; 
  outl[8] += incr[8]; 
  outl[9] += -1.0*incr[9]; 
  outl[10] += -1.0*incr[10]; 
  outl[11] += incr[11]; 
  outl[12] += incr[12]; 
  outl[13] += incr[13]; 
  outl[14] += -1.0*incr[14]; 
  outl[15] += incr[15]; 
  return std::abs(alpha0); 
} 
double EmGyrokineticSurfPositivity2x2vSer_Y_P1_Bvars_0(const double q_, const double m_, const double cflL, const double cflR, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *Apar, const double *dApardt, const double *fl, const double *fr, double *outl, double *outr) 
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
  double incr[16]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = (0.125*(2.0*BdriftY[0]*m_*wv2+q_*((BmagInv[0]*(3.0*Apar[3]-1.732050807568877*Apar[1])*dfac_x+BdriftY[0]*(Apar[0]-1.732050807568877*Apar[2]))*wv+BmagInv[0]*(1.732050807568877*Phi[1]-3.0*Phi[3])*dfac_x)))/q_; 

  double alpha[8]; 
  alpha[0] = (0.7071067811865475*(2.0*BdriftY[0]*m_*wv2+q_*((BmagInv[0]*(3.0*Apar[3]-1.732050807568877*Apar[1])*dfac_x+BdriftY[0]*(Apar[0]-1.732050807568877*Apar[2]))*wv+BmagInv[0]*(1.732050807568877*Phi[1]-3.0*Phi[3])*dfac_x)))/q_; 
  alpha[1] = -0.408248290463863*BdriftY[0]*(3.0*Apar[3]-1.732050807568877*Apar[1])*wv; 
  alpha[2] = (0.8164965809277261*BdriftY[0]*m_*wv)/(dfac_v*q_); 
  double rVal;  // rVal=f1/f0 at each control node in dimensions other than y 
  double fqVal[8];  // fqVal = anti-limited f evaluated at each control node on y surface 
  // determine upwinding at each surface control node 
  if(0.3535533905932737*alpha[0]-0.2041241452319315*(alpha[2]+alpha[1]) > 0) {
  rVal = -(1.0*(3.0*fl[15]-5.196152422706631*(fl[14]+fl[12]+fl[11])+9.0*(fl[9]+fl[7]+fl[5])-15.58845726811989*fl[2]))/(36.0*EPSILON-1.732050807568877*fl[13]+3.0*(fl[10]+fl[8]+fl[6])-5.196152422706631*(fl[4]+fl[3]+fl[1])+9.0*fl[0]); 
  fqVal[0] = -0.02777777777777778*(1.732050807568877*fl[13]-3.0*(fl[10]+fl[8]+fl[6])+5.196152422706631*(fl[4]+fl[3]+fl[1])-9.0*fl[0])*limTheta(rVal,1.0,cflL); 
  } else {
  rVal = -(1.0*(3.0*fr[15]-5.196152422706631*(fr[14]+fr[12]+fr[11])+9.0*(fr[9]+fr[7]+fr[5])-15.58845726811989*fr[2]))/(36.0*EPSILON-1.732050807568877*fr[13]+3.0*(fr[10]+fr[8]+fr[6])-5.196152422706631*(fr[4]+fr[3]+fr[1])+9.0*fr[0]); 
  fqVal[0] = -0.02777777777777778*(1.732050807568877*fr[13]-3.0*(fr[10]+fr[8]+fr[6])+5.196152422706631*(fr[4]+fr[3]+fr[1])-9.0*fr[0])*limTheta(rVal,-1.0,cflR); 
  }
  if((-0.2041241452319315*alpha[2])+0.2041241452319315*alpha[1]+0.3535533905932737*alpha[0] > 0) {
  rVal = (3.0*fl[15]+5.196152422706631*fl[14]-5.196152422706631*(fl[12]+fl[11])-9.0*(fl[9]+fl[7])+9.0*fl[5]+15.58845726811989*fl[2])/(36.0*EPSILON+1.732050807568877*fl[13]+3.0*fl[10]-3.0*(fl[8]+fl[6])-5.196152422706631*(fl[4]+fl[3])+5.196152422706631*fl[1]+9.0*fl[0]); 
  fqVal[1] = 0.02777777777777778*(1.732050807568877*fl[13]+3.0*fl[10]-3.0*(fl[8]+fl[6])-5.196152422706631*(fl[4]+fl[3])+5.196152422706631*fl[1]+9.0*fl[0])*limTheta(rVal,1.0,cflL); 
  } else {
  rVal = (3.0*fr[15]+5.196152422706631*fr[14]-5.196152422706631*(fr[12]+fr[11])-9.0*(fr[9]+fr[7])+9.0*fr[5]+15.58845726811989*fr[2])/(36.0*EPSILON+1.732050807568877*fr[13]+3.0*fr[10]-3.0*(fr[8]+fr[6])-5.196152422706631*(fr[4]+fr[3])+5.196152422706631*fr[1]+9.0*fr[0]); 
  fqVal[1] = 0.02777777777777778*(1.732050807568877*fr[13]+3.0*fr[10]-3.0*(fr[8]+fr[6])-5.196152422706631*(fr[4]+fr[3])+5.196152422706631*fr[1]+9.0*fr[0])*limTheta(rVal,-1.0,cflR); 
  }
  if(0.2041241452319315*alpha[2]-0.2041241452319315*alpha[1]+0.3535533905932737*alpha[0] > 0) {
  rVal = (3.0*fl[15]-5.196152422706631*fl[14]+5.196152422706631*fl[12]-5.196152422706631*fl[11]-9.0*fl[9]+9.0*fl[7]-9.0*fl[5]+15.58845726811989*fl[2])/(36.0*EPSILON+1.732050807568877*fl[13]-3.0*fl[10]+3.0*fl[8]-3.0*fl[6]-5.196152422706631*fl[4]+5.196152422706631*fl[3]-5.196152422706631*fl[1]+9.0*fl[0]); 
  fqVal[2] = 0.02777777777777778*(1.732050807568877*fl[13]-3.0*fl[10]+3.0*fl[8]-3.0*fl[6]-5.196152422706631*fl[4]+5.196152422706631*fl[3]-5.196152422706631*fl[1]+9.0*fl[0])*limTheta(rVal,1.0,cflL); 
  } else {
  rVal = (3.0*fr[15]-5.196152422706631*fr[14]+5.196152422706631*fr[12]-5.196152422706631*fr[11]-9.0*fr[9]+9.0*fr[7]-9.0*fr[5]+15.58845726811989*fr[2])/(36.0*EPSILON+1.732050807568877*fr[13]-3.0*fr[10]+3.0*fr[8]-3.0*fr[6]-5.196152422706631*fr[4]+5.196152422706631*fr[3]-5.196152422706631*fr[1]+9.0*fr[0]); 
  fqVal[2] = 0.02777777777777778*(1.732050807568877*fr[13]-3.0*fr[10]+3.0*fr[8]-3.0*fr[6]-5.196152422706631*fr[4]+5.196152422706631*fr[3]-5.196152422706631*fr[1]+9.0*fr[0])*limTheta(rVal,-1.0,cflR); 
  }
  if(0.2041241452319315*(alpha[2]+alpha[1])+0.3535533905932737*alpha[0] > 0) {
  rVal = -(1.0*(3.0*fl[15]+5.196152422706631*(fl[14]+fl[12])-5.196152422706631*fl[11]+9.0*fl[9]-9.0*(fl[7]+fl[5])-15.58845726811989*fl[2]))/(36.0*EPSILON-1.732050807568877*fl[13]-3.0*(fl[10]+fl[8])+3.0*fl[6]-5.196152422706631*fl[4]+5.196152422706631*(fl[3]+fl[1])+9.0*fl[0]); 
  fqVal[3] = -0.02777777777777778*(1.732050807568877*fl[13]+3.0*(fl[10]+fl[8])-3.0*fl[6]+5.196152422706631*fl[4]-5.196152422706631*(fl[3]+fl[1])-9.0*fl[0])*limTheta(rVal,1.0,cflL); 
  } else {
  rVal = -(1.0*(3.0*fr[15]+5.196152422706631*(fr[14]+fr[12])-5.196152422706631*fr[11]+9.0*fr[9]-9.0*(fr[7]+fr[5])-15.58845726811989*fr[2]))/(36.0*EPSILON-1.732050807568877*fr[13]-3.0*(fr[10]+fr[8])+3.0*fr[6]-5.196152422706631*fr[4]+5.196152422706631*(fr[3]+fr[1])+9.0*fr[0]); 
  fqVal[3] = -0.02777777777777778*(1.732050807568877*fr[13]+3.0*(fr[10]+fr[8])-3.0*fr[6]+5.196152422706631*fr[4]-5.196152422706631*(fr[3]+fr[1])-9.0*fr[0])*limTheta(rVal,-1.0,cflR); 
  }
  if(0.3535533905932737*alpha[0]-0.2041241452319315*(alpha[2]+alpha[1]) > 0) {
  rVal = (3.0*fl[15]-5.196152422706631*(fl[14]+fl[12])+5.196152422706631*fl[11]+9.0*fl[9]-9.0*(fl[7]+fl[5])+15.58845726811989*fl[2])/(36.0*EPSILON+1.732050807568877*fl[13]-3.0*(fl[10]+fl[8])+3.0*fl[6]+5.196152422706631*fl[4]-5.196152422706631*(fl[3]+fl[1])+9.0*fl[0]); 
  fqVal[4] = 0.02777777777777778*(1.732050807568877*fl[13]-3.0*(fl[10]+fl[8])+3.0*fl[6]+5.196152422706631*fl[4]-5.196152422706631*(fl[3]+fl[1])+9.0*fl[0])*limTheta(rVal,1.0,cflL); 
  } else {
  rVal = (3.0*fr[15]-5.196152422706631*(fr[14]+fr[12])+5.196152422706631*fr[11]+9.0*fr[9]-9.0*(fr[7]+fr[5])+15.58845726811989*fr[2])/(36.0*EPSILON+1.732050807568877*fr[13]-3.0*(fr[10]+fr[8])+3.0*fr[6]+5.196152422706631*fr[4]-5.196152422706631*(fr[3]+fr[1])+9.0*fr[0]); 
  fqVal[4] = 0.02777777777777778*(1.732050807568877*fr[13]-3.0*(fr[10]+fr[8])+3.0*fr[6]+5.196152422706631*fr[4]-5.196152422706631*(fr[3]+fr[1])+9.0*fr[0])*limTheta(rVal,-1.0,cflR); 
  }
  if((-0.2041241452319315*alpha[2])+0.2041241452319315*alpha[1]+0.3535533905932737*alpha[0] > 0) {
  rVal = -(1.0*(3.0*fl[15]+5.196152422706631*fl[14]-5.196152422706631*fl[12]+5.196152422706631*fl[11]-9.0*fl[9]+9.0*fl[7]-9.0*fl[5]-15.58845726811989*fl[2]))/(36.0*EPSILON-1.732050807568877*fl[13]-3.0*fl[10]+3.0*fl[8]-3.0*fl[6]+5.196152422706631*fl[4]-5.196152422706631*fl[3]+5.196152422706631*fl[1]+9.0*fl[0]); 
  fqVal[5] = -0.02777777777777778*(1.732050807568877*fl[13]+3.0*fl[10]-3.0*fl[8]+3.0*fl[6]-5.196152422706631*fl[4]+5.196152422706631*fl[3]-5.196152422706631*fl[1]-9.0*fl[0])*limTheta(rVal,1.0,cflL); 
  } else {
  rVal = -(1.0*(3.0*fr[15]+5.196152422706631*fr[14]-5.196152422706631*fr[12]+5.196152422706631*fr[11]-9.0*fr[9]+9.0*fr[7]-9.0*fr[5]-15.58845726811989*fr[2]))/(36.0*EPSILON-1.732050807568877*fr[13]-3.0*fr[10]+3.0*fr[8]-3.0*fr[6]+5.196152422706631*fr[4]-5.196152422706631*fr[3]+5.196152422706631*fr[1]+9.0*fr[0]); 
  fqVal[5] = -0.02777777777777778*(1.732050807568877*fr[13]+3.0*fr[10]-3.0*fr[8]+3.0*fr[6]-5.196152422706631*fr[4]+5.196152422706631*fr[3]-5.196152422706631*fr[1]-9.0*fr[0])*limTheta(rVal,-1.0,cflR); 
  }
  if(0.2041241452319315*alpha[2]-0.2041241452319315*alpha[1]+0.3535533905932737*alpha[0] > 0) {
  rVal = -(1.0*(3.0*fl[15]-5.196152422706631*fl[14]+5.196152422706631*(fl[12]+fl[11])-9.0*(fl[9]+fl[7])+9.0*fl[5]-15.58845726811989*fl[2]))/(36.0*EPSILON-1.732050807568877*fl[13]+3.0*fl[10]-3.0*(fl[8]+fl[6])+5.196152422706631*(fl[4]+fl[3])-5.196152422706631*fl[1]+9.0*fl[0]); 
  fqVal[6] = -0.02777777777777778*(1.732050807568877*fl[13]-3.0*fl[10]+3.0*(fl[8]+fl[6])-5.196152422706631*(fl[4]+fl[3])+5.196152422706631*fl[1]-9.0*fl[0])*limTheta(rVal,1.0,cflL); 
  } else {
  rVal = -(1.0*(3.0*fr[15]-5.196152422706631*fr[14]+5.196152422706631*(fr[12]+fr[11])-9.0*(fr[9]+fr[7])+9.0*fr[5]-15.58845726811989*fr[2]))/(36.0*EPSILON-1.732050807568877*fr[13]+3.0*fr[10]-3.0*(fr[8]+fr[6])+5.196152422706631*(fr[4]+fr[3])-5.196152422706631*fr[1]+9.0*fr[0]); 
  fqVal[6] = -0.02777777777777778*(1.732050807568877*fr[13]-3.0*fr[10]+3.0*(fr[8]+fr[6])-5.196152422706631*(fr[4]+fr[3])+5.196152422706631*fr[1]-9.0*fr[0])*limTheta(rVal,-1.0,cflR); 
  }
  if(0.2041241452319315*(alpha[2]+alpha[1])+0.3535533905932737*alpha[0] > 0) {
  rVal = (3.0*fl[15]+5.196152422706631*(fl[14]+fl[12]+fl[11])+9.0*(fl[9]+fl[7]+fl[5])+15.58845726811989*fl[2])/(36.0*EPSILON+1.732050807568877*fl[13]+3.0*(fl[10]+fl[8]+fl[6])+5.196152422706631*(fl[4]+fl[3]+fl[1])+9.0*fl[0]); 
  fqVal[7] = 0.02777777777777778*(1.732050807568877*fl[13]+3.0*(fl[10]+fl[8]+fl[6])+5.196152422706631*(fl[4]+fl[3]+fl[1])+9.0*fl[0])*limTheta(rVal,1.0,cflL); 
  } else {
  rVal = (3.0*fr[15]+5.196152422706631*(fr[14]+fr[12]+fr[11])+9.0*(fr[9]+fr[7]+fr[5])+15.58845726811989*fr[2])/(36.0*EPSILON+1.732050807568877*fr[13]+3.0*(fr[10]+fr[8]+fr[6])+5.196152422706631*(fr[4]+fr[3]+fr[1])+9.0*fr[0]); 
  fqVal[7] = 0.02777777777777778*(1.732050807568877*fr[13]+3.0*(fr[10]+fr[8]+fr[6])+5.196152422706631*(fr[4]+fr[3]+fr[1])+9.0*fr[0])*limTheta(rVal,-1.0,cflR); 
  }
  double fhatALVal[8];  // fhatALVal = mode coefficients of anti-limited f 
  fhatALVal[0] = 0.3535533905932737*(fqVal[7]+fqVal[6]+fqVal[5]+fqVal[4]+fqVal[3]+fqVal[2]+fqVal[1]+fqVal[0]); 
  fhatALVal[1] = 0.6123724356957944*(fqVal[7]-1.0*fqVal[6]+fqVal[5]-1.0*fqVal[4]+fqVal[3]-1.0*fqVal[2]+fqVal[1]-1.0*fqVal[0]); 
  fhatALVal[2] = 0.6123724356957944*(fqVal[7]+fqVal[6]-1.0*(fqVal[5]+fqVal[4])+fqVal[3]+fqVal[2]-1.0*(fqVal[1]+fqVal[0])); 
  fhatALVal[3] = 0.6123724356957944*(fqVal[7]+fqVal[6]+fqVal[5]+fqVal[4]-1.0*(fqVal[3]+fqVal[2]+fqVal[1]+fqVal[0])); 
  fhatALVal[4] = 1.060660171779821*(fqVal[7]-1.0*(fqVal[6]+fqVal[5])+fqVal[4]+fqVal[3]-1.0*(fqVal[2]+fqVal[1])+fqVal[0]); 
  fhatALVal[5] = 1.060660171779821*(fqVal[7]-1.0*fqVal[6]+fqVal[5]-1.0*(fqVal[4]+fqVal[3])+fqVal[2]-1.0*fqVal[1]+fqVal[0]); 
  fhatALVal[6] = 1.060660171779821*(fqVal[7]+fqVal[6]-1.0*(fqVal[5]+fqVal[4]+fqVal[3]+fqVal[2])+fqVal[1]+fqVal[0]); 
  fhatALVal[7] = 1.837117307087383*(fqVal[7]-1.0*(fqVal[6]+fqVal[5])+fqVal[4]-1.0*fqVal[3]+fqVal[2]+fqVal[1]-1.0*fqVal[0]); 
  incr[0] = 0.25*(alpha[2]*fhatALVal[2]+alpha[1]*fhatALVal[1]+alpha[0]*fhatALVal[0])*dfac_y; 
  incr[1] = 0.25*(alpha[2]*fhatALVal[4]+alpha[0]*fhatALVal[1]+fhatALVal[0]*alpha[1])*dfac_y; 
  incr[2] = -0.4330127018922193*(alpha[2]*fhatALVal[2]+alpha[1]*fhatALVal[1]+alpha[0]*fhatALVal[0])*dfac_y; 
  incr[3] = 0.25*(alpha[1]*fhatALVal[4]+alpha[0]*fhatALVal[2]+fhatALVal[0]*alpha[2])*dfac_y; 
  incr[4] = 0.25*(alpha[2]*fhatALVal[6]+alpha[1]*fhatALVal[5]+alpha[0]*fhatALVal[3])*dfac_y; 
  incr[5] = -0.4330127018922193*(alpha[2]*fhatALVal[4]+alpha[0]*fhatALVal[1]+fhatALVal[0]*alpha[1])*dfac_y; 
  incr[6] = 0.25*(alpha[0]*fhatALVal[4]+alpha[1]*fhatALVal[2]+fhatALVal[1]*alpha[2])*dfac_y; 
  incr[7] = -0.4330127018922193*(alpha[1]*fhatALVal[4]+alpha[0]*fhatALVal[2]+fhatALVal[0]*alpha[2])*dfac_y; 
  incr[8] = 0.25*(alpha[2]*fhatALVal[7]+alpha[0]*fhatALVal[5]+alpha[1]*fhatALVal[3])*dfac_y; 
  incr[9] = -0.4330127018922193*(alpha[2]*fhatALVal[6]+alpha[1]*fhatALVal[5]+alpha[0]*fhatALVal[3])*dfac_y; 
  incr[10] = 0.25*(alpha[1]*fhatALVal[7]+alpha[0]*fhatALVal[6]+alpha[2]*fhatALVal[3])*dfac_y; 
  incr[11] = -0.4330127018922193*(alpha[0]*fhatALVal[4]+alpha[1]*fhatALVal[2]+fhatALVal[1]*alpha[2])*dfac_y; 
  incr[12] = -0.4330127018922193*(alpha[2]*fhatALVal[7]+alpha[0]*fhatALVal[5]+alpha[1]*fhatALVal[3])*dfac_y; 
  incr[13] = 0.25*(alpha[0]*fhatALVal[7]+alpha[1]*fhatALVal[6]+alpha[2]*fhatALVal[5])*dfac_y; 
  incr[14] = -0.4330127018922193*(alpha[1]*fhatALVal[7]+alpha[0]*fhatALVal[6]+alpha[2]*fhatALVal[3])*dfac_y; 
  incr[15] = -0.4330127018922193*(alpha[0]*fhatALVal[7]+alpha[1]*fhatALVal[6]+alpha[2]*fhatALVal[5])*dfac_y; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 
  outr[4] += incr[4]; 
  outr[5] += incr[5]; 
  outr[6] += incr[6]; 
  outr[7] += incr[7]; 
  outr[8] += incr[8]; 
  outr[9] += incr[9]; 
  outr[10] += incr[10]; 
  outr[11] += incr[11]; 
  outr[12] += incr[12]; 
  outr[13] += incr[13]; 
  outr[14] += incr[14]; 
  outr[15] += incr[15]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += -1.0*incr[1]; 
  outl[2] += incr[2]; 
  outl[3] += -1.0*incr[3]; 
  outl[4] += -1.0*incr[4]; 
  outl[5] += incr[5]; 
  outl[6] += -1.0*incr[6]; 
  outl[7] += incr[7]; 
  outl[8] += -1.0*incr[8]; 
  outl[9] += incr[9]; 
  outl[10] += -1.0*incr[10]; 
  outl[11] += incr[11]; 
  outl[12] += incr[12]; 
  outl[13] += -1.0*incr[13]; 
  outl[14] += incr[14]; 
  outl[15] += incr[15]; 
  return std::abs(alpha0); 
} 
double EmGyrokineticSurfPositivity2x2vSer_Vpar_P1_Bvars_0(const double q_, const double m_, const double cflL, const double cflR, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *Apar, const double *dApardt, const double *fl, const double *fr, double *outl, double *outr) 
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
  double incr[16]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = -(0.0625*(dfac_v*(3.464101615137754*(BdriftY[0]*Phi[2]*dfac_y+BdriftX[0]*Phi[1]*dfac_x)*m_*wv+((3.0*BmagInv[0]*(Phi[1]*Apar[2]-1.0*Apar[1]*Phi[2])*dfac_x+1.732050807568877*BdriftY[0]*(Apar[1]*Phi[3]+Apar[0]*Phi[2]))*dfac_y+1.732050807568877*BdriftX[0]*(Apar[2]*Phi[3]+Apar[0]*Phi[1])*dfac_x+4.0*dApardt[0])*q_)-3.464101615137754*(BdriftY[0]*Phi[2]*dfac_y+BdriftX[0]*Phi[1]*dfac_x)*m_))/(dfac_v*m_); 

  double alpha[8]; 
  alpha[0] = -(0.3535533905932737*(dfac_v*(3.464101615137754*(BdriftY[0]*Phi[2]*dfac_y+BdriftX[0]*Phi[1]*dfac_x)*m_*wv+((BmagInv[0]*(3.0*Phi[1]*Apar[2]-3.0*Apar[1]*Phi[2])*dfac_x+1.732050807568877*BdriftY[0]*(Apar[1]*Phi[3]+Apar[0]*Phi[2]))*dfac_y+1.732050807568877*BdriftX[0]*(Apar[2]*Phi[3]+Apar[0]*Phi[1])*dfac_x+4.0*dApardt[0])*q_)-3.464101615137754*(BdriftY[0]*Phi[2]*dfac_y+BdriftX[0]*Phi[1]*dfac_x)*m_))/(dfac_v*m_); 
  alpha[1] = -(0.3535533905932737*(dfac_v*(3.464101615137754*BdriftY[0]*Phi[3]*dfac_y*m_*wv+((BmagInv[0]*(3.0*Phi[1]*Apar[3]-3.0*Apar[1]*Phi[3])*dfac_x+1.732050807568877*BdriftY[0]*(Apar[0]*Phi[3]+Apar[1]*Phi[2]))*dfac_y+1.732050807568877*BdriftX[0]*(Apar[3]*Phi[3]+Apar[1]*Phi[1])*dfac_x+4.0*dApardt[1])*q_)-3.464101615137754*BdriftY[0]*Phi[3]*dfac_y*m_))/(dfac_v*m_); 
  alpha[2] = -(0.3535533905932737*(dfac_v*(3.464101615137754*BdriftX[0]*Phi[3]*dfac_x*m_*wv+((BmagInv[0]*(3.0*Apar[2]*Phi[3]-3.0*Phi[2]*Apar[3])*dfac_x+1.732050807568877*BdriftY[0]*(Apar[3]*Phi[3]+Apar[2]*Phi[2]))*dfac_y+1.732050807568877*BdriftX[0]*(Apar[0]*Phi[3]+Phi[1]*Apar[2])*dfac_x+4.0*dApardt[2])*q_)-3.464101615137754*BdriftX[0]*Phi[3]*dfac_x*m_))/(dfac_v*m_); 
  alpha[4] = -(0.3535533905932737*(1.732050807568877*(BdriftY[0]*(Apar[2]*Phi[3]+Phi[2]*Apar[3])*dfac_y+BdriftX[0]*(Apar[1]*Phi[3]+Phi[1]*Apar[3])*dfac_x)+4.0*dApardt[3])*q_)/m_; 
  double rVal;  // rVal=f1/f0 at each control node in dimensions other than vx 
  double fqVal[8];  // fqVal = anti-limited f evaluated at each control node on vx surface 
  // determine upwinding at each surface control node 
  if(0.1178511301977579*alpha[4]-0.2041241452319315*(alpha[2]+alpha[1])+0.3535533905932737*alpha[0] > 0) {
  rVal = -(1.0*(3.0*fl[15]-5.196152422706631*(fl[14]+fl[13]+fl[11])+9.0*(fl[10]+fl[7]+fl[6])-15.58845726811989*fl[3]))/(36.0*EPSILON-1.732050807568877*fl[12]+3.0*(fl[9]+fl[8]+fl[5])-5.196152422706631*(fl[4]+fl[2]+fl[1])+9.0*fl[0]); 
  fqVal[0] = -0.02777777777777778*(1.732050807568877*fl[12]-3.0*(fl[9]+fl[8]+fl[5])+5.196152422706631*(fl[4]+fl[2]+fl[1])-9.0*fl[0])*limTheta(rVal,1.0,cflL); 
  } else {
  rVal = -(1.0*(3.0*fr[15]-5.196152422706631*(fr[14]+fr[13]+fr[11])+9.0*(fr[10]+fr[7]+fr[6])-15.58845726811989*fr[3]))/(36.0*EPSILON-1.732050807568877*fr[12]+3.0*(fr[9]+fr[8]+fr[5])-5.196152422706631*(fr[4]+fr[2]+fr[1])+9.0*fr[0]); 
  fqVal[0] = -0.02777777777777778*(1.732050807568877*fr[12]-3.0*(fr[9]+fr[8]+fr[5])+5.196152422706631*(fr[4]+fr[2]+fr[1])-9.0*fr[0])*limTheta(rVal,-1.0,cflR); 
  }
  if((-0.1178511301977579*alpha[4])-0.2041241452319315*alpha[2]+0.2041241452319315*alpha[1]+0.3535533905932737*alpha[0] > 0) {
  rVal = (3.0*fl[15]+5.196152422706631*fl[14]-5.196152422706631*(fl[13]+fl[11])-9.0*(fl[10]+fl[7])+9.0*fl[6]+15.58845726811989*fl[3])/(36.0*EPSILON+1.732050807568877*fl[12]+3.0*fl[9]-3.0*(fl[8]+fl[5])-5.196152422706631*(fl[4]+fl[2])+5.196152422706631*fl[1]+9.0*fl[0]); 
  fqVal[1] = 0.02777777777777778*(1.732050807568877*fl[12]+3.0*fl[9]-3.0*(fl[8]+fl[5])-5.196152422706631*(fl[4]+fl[2])+5.196152422706631*fl[1]+9.0*fl[0])*limTheta(rVal,1.0,cflL); 
  } else {
  rVal = (3.0*fr[15]+5.196152422706631*fr[14]-5.196152422706631*(fr[13]+fr[11])-9.0*(fr[10]+fr[7])+9.0*fr[6]+15.58845726811989*fr[3])/(36.0*EPSILON+1.732050807568877*fr[12]+3.0*fr[9]-3.0*(fr[8]+fr[5])-5.196152422706631*(fr[4]+fr[2])+5.196152422706631*fr[1]+9.0*fr[0]); 
  fqVal[1] = 0.02777777777777778*(1.732050807568877*fr[12]+3.0*fr[9]-3.0*(fr[8]+fr[5])-5.196152422706631*(fr[4]+fr[2])+5.196152422706631*fr[1]+9.0*fr[0])*limTheta(rVal,-1.0,cflR); 
  }
  if((-0.1178511301977579*alpha[4])+0.2041241452319315*alpha[2]-0.2041241452319315*alpha[1]+0.3535533905932737*alpha[0] > 0) {
  rVal = (3.0*fl[15]-5.196152422706631*fl[14]+5.196152422706631*fl[13]-5.196152422706631*fl[11]-9.0*fl[10]+9.0*fl[7]-9.0*fl[6]+15.58845726811989*fl[3])/(36.0*EPSILON+1.732050807568877*fl[12]-3.0*fl[9]+3.0*fl[8]-3.0*fl[5]-5.196152422706631*fl[4]+5.196152422706631*fl[2]-5.196152422706631*fl[1]+9.0*fl[0]); 
  fqVal[2] = 0.02777777777777778*(1.732050807568877*fl[12]-3.0*fl[9]+3.0*fl[8]-3.0*fl[5]-5.196152422706631*fl[4]+5.196152422706631*fl[2]-5.196152422706631*fl[1]+9.0*fl[0])*limTheta(rVal,1.0,cflL); 
  } else {
  rVal = (3.0*fr[15]-5.196152422706631*fr[14]+5.196152422706631*fr[13]-5.196152422706631*fr[11]-9.0*fr[10]+9.0*fr[7]-9.0*fr[6]+15.58845726811989*fr[3])/(36.0*EPSILON+1.732050807568877*fr[12]-3.0*fr[9]+3.0*fr[8]-3.0*fr[5]-5.196152422706631*fr[4]+5.196152422706631*fr[2]-5.196152422706631*fr[1]+9.0*fr[0]); 
  fqVal[2] = 0.02777777777777778*(1.732050807568877*fr[12]-3.0*fr[9]+3.0*fr[8]-3.0*fr[5]-5.196152422706631*fr[4]+5.196152422706631*fr[2]-5.196152422706631*fr[1]+9.0*fr[0])*limTheta(rVal,-1.0,cflR); 
  }
  if(0.1178511301977579*alpha[4]+0.2041241452319315*(alpha[2]+alpha[1])+0.3535533905932737*alpha[0] > 0) {
  rVal = -(1.0*(3.0*fl[15]+5.196152422706631*(fl[14]+fl[13])-5.196152422706631*fl[11]+9.0*fl[10]-9.0*(fl[7]+fl[6])-15.58845726811989*fl[3]))/(36.0*EPSILON-1.732050807568877*fl[12]-3.0*(fl[9]+fl[8])+3.0*fl[5]-5.196152422706631*fl[4]+5.196152422706631*(fl[2]+fl[1])+9.0*fl[0]); 
  fqVal[3] = -0.02777777777777778*(1.732050807568877*fl[12]+3.0*(fl[9]+fl[8])-3.0*fl[5]+5.196152422706631*fl[4]-5.196152422706631*(fl[2]+fl[1])-9.0*fl[0])*limTheta(rVal,1.0,cflL); 
  } else {
  rVal = -(1.0*(3.0*fr[15]+5.196152422706631*(fr[14]+fr[13])-5.196152422706631*fr[11]+9.0*fr[10]-9.0*(fr[7]+fr[6])-15.58845726811989*fr[3]))/(36.0*EPSILON-1.732050807568877*fr[12]-3.0*(fr[9]+fr[8])+3.0*fr[5]-5.196152422706631*fr[4]+5.196152422706631*(fr[2]+fr[1])+9.0*fr[0]); 
  fqVal[3] = -0.02777777777777778*(1.732050807568877*fr[12]+3.0*(fr[9]+fr[8])-3.0*fr[5]+5.196152422706631*fr[4]-5.196152422706631*(fr[2]+fr[1])-9.0*fr[0])*limTheta(rVal,-1.0,cflR); 
  }
  if(0.1178511301977579*alpha[4]-0.2041241452319315*(alpha[2]+alpha[1])+0.3535533905932737*alpha[0] > 0) {
  rVal = (3.0*fl[15]-5.196152422706631*(fl[14]+fl[13])+5.196152422706631*fl[11]+9.0*fl[10]-9.0*(fl[7]+fl[6])+15.58845726811989*fl[3])/(36.0*EPSILON+1.732050807568877*fl[12]-3.0*(fl[9]+fl[8])+3.0*fl[5]+5.196152422706631*fl[4]-5.196152422706631*(fl[2]+fl[1])+9.0*fl[0]); 
  fqVal[4] = 0.02777777777777778*(1.732050807568877*fl[12]-3.0*(fl[9]+fl[8])+3.0*fl[5]+5.196152422706631*fl[4]-5.196152422706631*(fl[2]+fl[1])+9.0*fl[0])*limTheta(rVal,1.0,cflL); 
  } else {
  rVal = (3.0*fr[15]-5.196152422706631*(fr[14]+fr[13])+5.196152422706631*fr[11]+9.0*fr[10]-9.0*(fr[7]+fr[6])+15.58845726811989*fr[3])/(36.0*EPSILON+1.732050807568877*fr[12]-3.0*(fr[9]+fr[8])+3.0*fr[5]+5.196152422706631*fr[4]-5.196152422706631*(fr[2]+fr[1])+9.0*fr[0]); 
  fqVal[4] = 0.02777777777777778*(1.732050807568877*fr[12]-3.0*(fr[9]+fr[8])+3.0*fr[5]+5.196152422706631*fr[4]-5.196152422706631*(fr[2]+fr[1])+9.0*fr[0])*limTheta(rVal,-1.0,cflR); 
  }
  if((-0.1178511301977579*alpha[4])-0.2041241452319315*alpha[2]+0.2041241452319315*alpha[1]+0.3535533905932737*alpha[0] > 0) {
  rVal = -(1.0*(3.0*fl[15]+5.196152422706631*fl[14]-5.196152422706631*fl[13]+5.196152422706631*fl[11]-9.0*fl[10]+9.0*fl[7]-9.0*fl[6]-15.58845726811989*fl[3]))/(36.0*EPSILON-1.732050807568877*fl[12]-3.0*fl[9]+3.0*fl[8]-3.0*fl[5]+5.196152422706631*fl[4]-5.196152422706631*fl[2]+5.196152422706631*fl[1]+9.0*fl[0]); 
  fqVal[5] = -0.02777777777777778*(1.732050807568877*fl[12]+3.0*fl[9]-3.0*fl[8]+3.0*fl[5]-5.196152422706631*fl[4]+5.196152422706631*fl[2]-5.196152422706631*fl[1]-9.0*fl[0])*limTheta(rVal,1.0,cflL); 
  } else {
  rVal = -(1.0*(3.0*fr[15]+5.196152422706631*fr[14]-5.196152422706631*fr[13]+5.196152422706631*fr[11]-9.0*fr[10]+9.0*fr[7]-9.0*fr[6]-15.58845726811989*fr[3]))/(36.0*EPSILON-1.732050807568877*fr[12]-3.0*fr[9]+3.0*fr[8]-3.0*fr[5]+5.196152422706631*fr[4]-5.196152422706631*fr[2]+5.196152422706631*fr[1]+9.0*fr[0]); 
  fqVal[5] = -0.02777777777777778*(1.732050807568877*fr[12]+3.0*fr[9]-3.0*fr[8]+3.0*fr[5]-5.196152422706631*fr[4]+5.196152422706631*fr[2]-5.196152422706631*fr[1]-9.0*fr[0])*limTheta(rVal,-1.0,cflR); 
  }
  if((-0.1178511301977579*alpha[4])+0.2041241452319315*alpha[2]-0.2041241452319315*alpha[1]+0.3535533905932737*alpha[0] > 0) {
  rVal = -(1.0*(3.0*fl[15]-5.196152422706631*fl[14]+5.196152422706631*(fl[13]+fl[11])-9.0*(fl[10]+fl[7])+9.0*fl[6]-15.58845726811989*fl[3]))/(36.0*EPSILON-1.732050807568877*fl[12]+3.0*fl[9]-3.0*(fl[8]+fl[5])+5.196152422706631*(fl[4]+fl[2])-5.196152422706631*fl[1]+9.0*fl[0]); 
  fqVal[6] = -0.02777777777777778*(1.732050807568877*fl[12]-3.0*fl[9]+3.0*(fl[8]+fl[5])-5.196152422706631*(fl[4]+fl[2])+5.196152422706631*fl[1]-9.0*fl[0])*limTheta(rVal,1.0,cflL); 
  } else {
  rVal = -(1.0*(3.0*fr[15]-5.196152422706631*fr[14]+5.196152422706631*(fr[13]+fr[11])-9.0*(fr[10]+fr[7])+9.0*fr[6]-15.58845726811989*fr[3]))/(36.0*EPSILON-1.732050807568877*fr[12]+3.0*fr[9]-3.0*(fr[8]+fr[5])+5.196152422706631*(fr[4]+fr[2])-5.196152422706631*fr[1]+9.0*fr[0]); 
  fqVal[6] = -0.02777777777777778*(1.732050807568877*fr[12]-3.0*fr[9]+3.0*(fr[8]+fr[5])-5.196152422706631*(fr[4]+fr[2])+5.196152422706631*fr[1]-9.0*fr[0])*limTheta(rVal,-1.0,cflR); 
  }
  if(0.1178511301977579*alpha[4]+0.2041241452319315*(alpha[2]+alpha[1])+0.3535533905932737*alpha[0] > 0) {
  rVal = (3.0*fl[15]+5.196152422706631*(fl[14]+fl[13]+fl[11])+9.0*(fl[10]+fl[7]+fl[6])+15.58845726811989*fl[3])/(36.0*EPSILON+1.732050807568877*fl[12]+3.0*(fl[9]+fl[8]+fl[5])+5.196152422706631*(fl[4]+fl[2]+fl[1])+9.0*fl[0]); 
  fqVal[7] = 0.02777777777777778*(1.732050807568877*fl[12]+3.0*(fl[9]+fl[8]+fl[5])+5.196152422706631*(fl[4]+fl[2]+fl[1])+9.0*fl[0])*limTheta(rVal,1.0,cflL); 
  } else {
  rVal = (3.0*fr[15]+5.196152422706631*(fr[14]+fr[13]+fr[11])+9.0*(fr[10]+fr[7]+fr[6])+15.58845726811989*fr[3])/(36.0*EPSILON+1.732050807568877*fr[12]+3.0*(fr[9]+fr[8]+fr[5])+5.196152422706631*(fr[4]+fr[2]+fr[1])+9.0*fr[0]); 
  fqVal[7] = 0.02777777777777778*(1.732050807568877*fr[12]+3.0*(fr[9]+fr[8]+fr[5])+5.196152422706631*(fr[4]+fr[2]+fr[1])+9.0*fr[0])*limTheta(rVal,-1.0,cflR); 
  }
  double fhatALVal[8];  // fhatALVal = mode coefficients of anti-limited f 
  fhatALVal[0] = 0.3535533905932737*(fqVal[7]+fqVal[6]+fqVal[5]+fqVal[4]+fqVal[3]+fqVal[2]+fqVal[1]+fqVal[0]); 
  fhatALVal[1] = 0.6123724356957944*(fqVal[7]-1.0*fqVal[6]+fqVal[5]-1.0*fqVal[4]+fqVal[3]-1.0*fqVal[2]+fqVal[1]-1.0*fqVal[0]); 
  fhatALVal[2] = 0.6123724356957944*(fqVal[7]+fqVal[6]-1.0*(fqVal[5]+fqVal[4])+fqVal[3]+fqVal[2]-1.0*(fqVal[1]+fqVal[0])); 
  fhatALVal[3] = 0.6123724356957944*(fqVal[7]+fqVal[6]+fqVal[5]+fqVal[4]-1.0*(fqVal[3]+fqVal[2]+fqVal[1]+fqVal[0])); 
  fhatALVal[4] = 1.060660171779821*(fqVal[7]-1.0*(fqVal[6]+fqVal[5])+fqVal[4]+fqVal[3]-1.0*(fqVal[2]+fqVal[1])+fqVal[0]); 
  fhatALVal[5] = 1.060660171779821*(fqVal[7]-1.0*fqVal[6]+fqVal[5]-1.0*(fqVal[4]+fqVal[3])+fqVal[2]-1.0*fqVal[1]+fqVal[0]); 
  fhatALVal[6] = 1.060660171779821*(fqVal[7]+fqVal[6]-1.0*(fqVal[5]+fqVal[4]+fqVal[3]+fqVal[2])+fqVal[1]+fqVal[0]); 
  fhatALVal[7] = 1.837117307087383*(fqVal[7]-1.0*(fqVal[6]+fqVal[5])+fqVal[4]-1.0*fqVal[3]+fqVal[2]+fqVal[1]-1.0*fqVal[0]); 
  incr[0] = 0.25*(alpha[4]*fhatALVal[4]+alpha[2]*fhatALVal[2]+alpha[1]*fhatALVal[1]+alpha[0]*fhatALVal[0])*dfac_v; 
  incr[1] = 0.25*(alpha[2]*fhatALVal[4]+fhatALVal[2]*alpha[4]+alpha[0]*fhatALVal[1]+fhatALVal[0]*alpha[1])*dfac_v; 
  incr[2] = 0.25*(alpha[1]*fhatALVal[4]+fhatALVal[1]*alpha[4]+alpha[0]*fhatALVal[2]+fhatALVal[0]*alpha[2])*dfac_v; 
  incr[3] = -0.4330127018922193*(alpha[4]*fhatALVal[4]+alpha[2]*fhatALVal[2]+alpha[1]*fhatALVal[1]+alpha[0]*fhatALVal[0])*dfac_v; 
  incr[4] = 0.25*(alpha[4]*fhatALVal[7]+alpha[2]*fhatALVal[6]+alpha[1]*fhatALVal[5]+alpha[0]*fhatALVal[3])*dfac_v; 
  incr[5] = 0.25*(alpha[0]*fhatALVal[4]+fhatALVal[0]*alpha[4]+alpha[1]*fhatALVal[2]+fhatALVal[1]*alpha[2])*dfac_v; 
  incr[6] = -0.4330127018922193*(alpha[2]*fhatALVal[4]+fhatALVal[2]*alpha[4]+alpha[0]*fhatALVal[1]+fhatALVal[0]*alpha[1])*dfac_v; 
  incr[7] = -0.4330127018922193*(alpha[1]*fhatALVal[4]+fhatALVal[1]*alpha[4]+alpha[0]*fhatALVal[2]+fhatALVal[0]*alpha[2])*dfac_v; 
  incr[8] = 0.25*(alpha[2]*fhatALVal[7]+alpha[4]*fhatALVal[6]+alpha[0]*fhatALVal[5]+alpha[1]*fhatALVal[3])*dfac_v; 
  incr[9] = 0.25*(alpha[1]*fhatALVal[7]+alpha[0]*fhatALVal[6]+alpha[4]*fhatALVal[5]+alpha[2]*fhatALVal[3])*dfac_v; 
  incr[10] = -0.4330127018922193*(alpha[4]*fhatALVal[7]+alpha[2]*fhatALVal[6]+alpha[1]*fhatALVal[5]+alpha[0]*fhatALVal[3])*dfac_v; 
  incr[11] = -0.4330127018922193*(alpha[0]*fhatALVal[4]+fhatALVal[0]*alpha[4]+alpha[1]*fhatALVal[2]+fhatALVal[1]*alpha[2])*dfac_v; 
  incr[12] = 0.25*(alpha[0]*fhatALVal[7]+alpha[1]*fhatALVal[6]+alpha[2]*fhatALVal[5]+fhatALVal[3]*alpha[4])*dfac_v; 
  incr[13] = -0.4330127018922193*(alpha[2]*fhatALVal[7]+alpha[4]*fhatALVal[6]+alpha[0]*fhatALVal[5]+alpha[1]*fhatALVal[3])*dfac_v; 
  incr[14] = -0.4330127018922193*(alpha[1]*fhatALVal[7]+alpha[0]*fhatALVal[6]+alpha[4]*fhatALVal[5]+alpha[2]*fhatALVal[3])*dfac_v; 
  incr[15] = -0.4330127018922193*(alpha[0]*fhatALVal[7]+alpha[1]*fhatALVal[6]+alpha[2]*fhatALVal[5]+fhatALVal[3]*alpha[4])*dfac_v; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 
  outr[4] += incr[4]; 
  outr[5] += incr[5]; 
  outr[6] += incr[6]; 
  outr[7] += incr[7]; 
  outr[8] += incr[8]; 
  outr[9] += incr[9]; 
  outr[10] += incr[10]; 
  outr[11] += incr[11]; 
  outr[12] += incr[12]; 
  outr[13] += incr[13]; 
  outr[14] += incr[14]; 
  outr[15] += incr[15]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += -1.0*incr[1]; 
  outl[2] += -1.0*incr[2]; 
  outl[3] += incr[3]; 
  outl[4] += -1.0*incr[4]; 
  outl[5] += -1.0*incr[5]; 
  outl[6] += incr[6]; 
  outl[7] += incr[7]; 
  outl[8] += -1.0*incr[8]; 
  outl[9] += -1.0*incr[9]; 
  outl[10] += incr[10]; 
  outl[11] += incr[11]; 
  outl[12] += -1.0*incr[12]; 
  outl[13] += incr[13]; 
  outl[14] += incr[14]; 
  outl[15] += incr[15]; 
  return std::abs(alpha0); 
} 
double EmGyrokineticSurfPositivity2x2vSer_X_P1_Bvars_1(const double q_, const double m_, const double cflL, const double cflR, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *Apar, const double *dApardt, const double *fl, const double *fr, double *outl, double *outr) 
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
  double incr[16]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = -(0.125*(2.0*(1.732050807568877*BdriftX[1]-1.0*BdriftX[0])*m_*wv2+q_*((((3.0*BmagInv[0]-5.196152422706631*BmagInv[1])*Apar[3]+(3.0*BmagInv[1]-1.732050807568877*BmagInv[0])*Apar[2])*dfac_y+(1.732050807568877*Apar[0]-3.0*Apar[1])*BdriftX[1]+BdriftX[0]*(1.732050807568877*Apar[1]-1.0*Apar[0]))*wv+((5.196152422706631*BmagInv[1]-3.0*BmagInv[0])*Phi[3]+(1.732050807568877*BmagInv[0]-3.0*BmagInv[1])*Phi[2])*dfac_y)))/q_; 

  double alpha[8]; 
  alpha[0] = -(0.7071067811865475*((3.464101615137754*BdriftX[1]-2.0*BdriftX[0])*m_*wv2+q_*((((3.0*BmagInv[0]-5.196152422706631*BmagInv[1])*Apar[3]+(3.0*BmagInv[1]-1.732050807568877*BmagInv[0])*Apar[2])*dfac_y+(1.732050807568877*Apar[0]-3.0*Apar[1])*BdriftX[1]+BdriftX[0]*(1.732050807568877*Apar[1]-1.0*Apar[0]))*wv+((5.196152422706631*BmagInv[1]-3.0*BmagInv[0])*Phi[3]+(1.732050807568877*BmagInv[0]-3.0*BmagInv[1])*Phi[2])*dfac_y)))/q_; 
  alpha[1] = 0.408248290463863*((5.196152422706631*BdriftX[1]-3.0*BdriftX[0])*Apar[3]+(1.732050807568877*BdriftX[0]-3.0*BdriftX[1])*Apar[2])*wv; 
  alpha[2] = -(0.2357022603955158*(6.0*BdriftX[1]-3.464101615137754*BdriftX[0])*m_*wv)/(dfac_v*q_); 
  double rVal;  // rVal=f1/f0 at each control node in dimensions other than x 
  double fqVal[8];  // fqVal = anti-limited f evaluated at each control node on x surface 
  // determine upwinding at each surface control node 
  if(0.3535533905932737*alpha[0]-0.2041241452319315*(alpha[2]+alpha[1]) > 0) {
  rVal = -(1.0*(3.0*fl[15]-5.196152422706631*(fl[13]+fl[12]+fl[11])+9.0*(fl[8]+fl[6]+fl[5])-15.58845726811989*fl[1]))/(36.0*EPSILON-1.732050807568877*fl[14]+3.0*(fl[10]+fl[9]+fl[7])-5.196152422706631*(fl[4]+fl[3]+fl[2])+9.0*fl[0]); 
  fqVal[0] = -0.02777777777777778*(1.732050807568877*fl[14]-3.0*(fl[10]+fl[9]+fl[7])+5.196152422706631*(fl[4]+fl[3]+fl[2])-9.0*fl[0])*limTheta(rVal,1.0,cflL); 
  } else {
  rVal = -(1.0*(3.0*fr[15]-5.196152422706631*(fr[13]+fr[12]+fr[11])+9.0*(fr[8]+fr[6]+fr[5])-15.58845726811989*fr[1]))/(36.0*EPSILON-1.732050807568877*fr[14]+3.0*(fr[10]+fr[9]+fr[7])-5.196152422706631*(fr[4]+fr[3]+fr[2])+9.0*fr[0]); 
  fqVal[0] = -0.02777777777777778*(1.732050807568877*fr[14]-3.0*(fr[10]+fr[9]+fr[7])+5.196152422706631*(fr[4]+fr[3]+fr[2])-9.0*fr[0])*limTheta(rVal,-1.0,cflR); 
  }
  if((-0.2041241452319315*alpha[2])+0.2041241452319315*alpha[1]+0.3535533905932737*alpha[0] > 0) {
  rVal = (3.0*fl[15]+5.196152422706631*fl[13]-5.196152422706631*(fl[12]+fl[11])-9.0*(fl[8]+fl[6])+9.0*fl[5]+15.58845726811989*fl[1])/(36.0*EPSILON+1.732050807568877*fl[14]+3.0*fl[10]-3.0*(fl[9]+fl[7])-5.196152422706631*(fl[4]+fl[3])+5.196152422706631*fl[2]+9.0*fl[0]); 
  fqVal[1] = 0.02777777777777778*(1.732050807568877*fl[14]+3.0*fl[10]-3.0*(fl[9]+fl[7])-5.196152422706631*(fl[4]+fl[3])+5.196152422706631*fl[2]+9.0*fl[0])*limTheta(rVal,1.0,cflL); 
  } else {
  rVal = (3.0*fr[15]+5.196152422706631*fr[13]-5.196152422706631*(fr[12]+fr[11])-9.0*(fr[8]+fr[6])+9.0*fr[5]+15.58845726811989*fr[1])/(36.0*EPSILON+1.732050807568877*fr[14]+3.0*fr[10]-3.0*(fr[9]+fr[7])-5.196152422706631*(fr[4]+fr[3])+5.196152422706631*fr[2]+9.0*fr[0]); 
  fqVal[1] = 0.02777777777777778*(1.732050807568877*fr[14]+3.0*fr[10]-3.0*(fr[9]+fr[7])-5.196152422706631*(fr[4]+fr[3])+5.196152422706631*fr[2]+9.0*fr[0])*limTheta(rVal,-1.0,cflR); 
  }
  if(0.2041241452319315*alpha[2]-0.2041241452319315*alpha[1]+0.3535533905932737*alpha[0] > 0) {
  rVal = (3.0*fl[15]-5.196152422706631*fl[13]+5.196152422706631*fl[12]-5.196152422706631*fl[11]-9.0*fl[8]+9.0*fl[6]-9.0*fl[5]+15.58845726811989*fl[1])/(36.0*EPSILON+1.732050807568877*fl[14]-3.0*fl[10]+3.0*fl[9]-3.0*fl[7]-5.196152422706631*fl[4]+5.196152422706631*fl[3]-5.196152422706631*fl[2]+9.0*fl[0]); 
  fqVal[2] = 0.02777777777777778*(1.732050807568877*fl[14]-3.0*fl[10]+3.0*fl[9]-3.0*fl[7]-5.196152422706631*fl[4]+5.196152422706631*fl[3]-5.196152422706631*fl[2]+9.0*fl[0])*limTheta(rVal,1.0,cflL); 
  } else {
  rVal = (3.0*fr[15]-5.196152422706631*fr[13]+5.196152422706631*fr[12]-5.196152422706631*fr[11]-9.0*fr[8]+9.0*fr[6]-9.0*fr[5]+15.58845726811989*fr[1])/(36.0*EPSILON+1.732050807568877*fr[14]-3.0*fr[10]+3.0*fr[9]-3.0*fr[7]-5.196152422706631*fr[4]+5.196152422706631*fr[3]-5.196152422706631*fr[2]+9.0*fr[0]); 
  fqVal[2] = 0.02777777777777778*(1.732050807568877*fr[14]-3.0*fr[10]+3.0*fr[9]-3.0*fr[7]-5.196152422706631*fr[4]+5.196152422706631*fr[3]-5.196152422706631*fr[2]+9.0*fr[0])*limTheta(rVal,-1.0,cflR); 
  }
  if(0.2041241452319315*(alpha[2]+alpha[1])+0.3535533905932737*alpha[0] > 0) {
  rVal = -(1.0*(3.0*fl[15]+5.196152422706631*(fl[13]+fl[12])-5.196152422706631*fl[11]+9.0*fl[8]-9.0*(fl[6]+fl[5])-15.58845726811989*fl[1]))/(36.0*EPSILON-1.732050807568877*fl[14]-3.0*(fl[10]+fl[9])+3.0*fl[7]-5.196152422706631*fl[4]+5.196152422706631*(fl[3]+fl[2])+9.0*fl[0]); 
  fqVal[3] = -0.02777777777777778*(1.732050807568877*fl[14]+3.0*(fl[10]+fl[9])-3.0*fl[7]+5.196152422706631*fl[4]-5.196152422706631*(fl[3]+fl[2])-9.0*fl[0])*limTheta(rVal,1.0,cflL); 
  } else {
  rVal = -(1.0*(3.0*fr[15]+5.196152422706631*(fr[13]+fr[12])-5.196152422706631*fr[11]+9.0*fr[8]-9.0*(fr[6]+fr[5])-15.58845726811989*fr[1]))/(36.0*EPSILON-1.732050807568877*fr[14]-3.0*(fr[10]+fr[9])+3.0*fr[7]-5.196152422706631*fr[4]+5.196152422706631*(fr[3]+fr[2])+9.0*fr[0]); 
  fqVal[3] = -0.02777777777777778*(1.732050807568877*fr[14]+3.0*(fr[10]+fr[9])-3.0*fr[7]+5.196152422706631*fr[4]-5.196152422706631*(fr[3]+fr[2])-9.0*fr[0])*limTheta(rVal,-1.0,cflR); 
  }
  if(0.3535533905932737*alpha[0]-0.2041241452319315*(alpha[2]+alpha[1]) > 0) {
  rVal = (3.0*fl[15]-5.196152422706631*(fl[13]+fl[12])+5.196152422706631*fl[11]+9.0*fl[8]-9.0*(fl[6]+fl[5])+15.58845726811989*fl[1])/(36.0*EPSILON+1.732050807568877*fl[14]-3.0*(fl[10]+fl[9])+3.0*fl[7]+5.196152422706631*fl[4]-5.196152422706631*(fl[3]+fl[2])+9.0*fl[0]); 
  fqVal[4] = 0.02777777777777778*(1.732050807568877*fl[14]-3.0*(fl[10]+fl[9])+3.0*fl[7]+5.196152422706631*fl[4]-5.196152422706631*(fl[3]+fl[2])+9.0*fl[0])*limTheta(rVal,1.0,cflL); 
  } else {
  rVal = (3.0*fr[15]-5.196152422706631*(fr[13]+fr[12])+5.196152422706631*fr[11]+9.0*fr[8]-9.0*(fr[6]+fr[5])+15.58845726811989*fr[1])/(36.0*EPSILON+1.732050807568877*fr[14]-3.0*(fr[10]+fr[9])+3.0*fr[7]+5.196152422706631*fr[4]-5.196152422706631*(fr[3]+fr[2])+9.0*fr[0]); 
  fqVal[4] = 0.02777777777777778*(1.732050807568877*fr[14]-3.0*(fr[10]+fr[9])+3.0*fr[7]+5.196152422706631*fr[4]-5.196152422706631*(fr[3]+fr[2])+9.0*fr[0])*limTheta(rVal,-1.0,cflR); 
  }
  if((-0.2041241452319315*alpha[2])+0.2041241452319315*alpha[1]+0.3535533905932737*alpha[0] > 0) {
  rVal = -(1.0*(3.0*fl[15]+5.196152422706631*fl[13]-5.196152422706631*fl[12]+5.196152422706631*fl[11]-9.0*fl[8]+9.0*fl[6]-9.0*fl[5]-15.58845726811989*fl[1]))/(36.0*EPSILON-1.732050807568877*fl[14]-3.0*fl[10]+3.0*fl[9]-3.0*fl[7]+5.196152422706631*fl[4]-5.196152422706631*fl[3]+5.196152422706631*fl[2]+9.0*fl[0]); 
  fqVal[5] = -0.02777777777777778*(1.732050807568877*fl[14]+3.0*fl[10]-3.0*fl[9]+3.0*fl[7]-5.196152422706631*fl[4]+5.196152422706631*fl[3]-5.196152422706631*fl[2]-9.0*fl[0])*limTheta(rVal,1.0,cflL); 
  } else {
  rVal = -(1.0*(3.0*fr[15]+5.196152422706631*fr[13]-5.196152422706631*fr[12]+5.196152422706631*fr[11]-9.0*fr[8]+9.0*fr[6]-9.0*fr[5]-15.58845726811989*fr[1]))/(36.0*EPSILON-1.732050807568877*fr[14]-3.0*fr[10]+3.0*fr[9]-3.0*fr[7]+5.196152422706631*fr[4]-5.196152422706631*fr[3]+5.196152422706631*fr[2]+9.0*fr[0]); 
  fqVal[5] = -0.02777777777777778*(1.732050807568877*fr[14]+3.0*fr[10]-3.0*fr[9]+3.0*fr[7]-5.196152422706631*fr[4]+5.196152422706631*fr[3]-5.196152422706631*fr[2]-9.0*fr[0])*limTheta(rVal,-1.0,cflR); 
  }
  if(0.2041241452319315*alpha[2]-0.2041241452319315*alpha[1]+0.3535533905932737*alpha[0] > 0) {
  rVal = -(1.0*(3.0*fl[15]-5.196152422706631*fl[13]+5.196152422706631*(fl[12]+fl[11])-9.0*(fl[8]+fl[6])+9.0*fl[5]-15.58845726811989*fl[1]))/(36.0*EPSILON-1.732050807568877*fl[14]+3.0*fl[10]-3.0*(fl[9]+fl[7])+5.196152422706631*(fl[4]+fl[3])-5.196152422706631*fl[2]+9.0*fl[0]); 
  fqVal[6] = -0.02777777777777778*(1.732050807568877*fl[14]-3.0*fl[10]+3.0*(fl[9]+fl[7])-5.196152422706631*(fl[4]+fl[3])+5.196152422706631*fl[2]-9.0*fl[0])*limTheta(rVal,1.0,cflL); 
  } else {
  rVal = -(1.0*(3.0*fr[15]-5.196152422706631*fr[13]+5.196152422706631*(fr[12]+fr[11])-9.0*(fr[8]+fr[6])+9.0*fr[5]-15.58845726811989*fr[1]))/(36.0*EPSILON-1.732050807568877*fr[14]+3.0*fr[10]-3.0*(fr[9]+fr[7])+5.196152422706631*(fr[4]+fr[3])-5.196152422706631*fr[2]+9.0*fr[0]); 
  fqVal[6] = -0.02777777777777778*(1.732050807568877*fr[14]-3.0*fr[10]+3.0*(fr[9]+fr[7])-5.196152422706631*(fr[4]+fr[3])+5.196152422706631*fr[2]-9.0*fr[0])*limTheta(rVal,-1.0,cflR); 
  }
  if(0.2041241452319315*(alpha[2]+alpha[1])+0.3535533905932737*alpha[0] > 0) {
  rVal = (3.0*fl[15]+5.196152422706631*(fl[13]+fl[12]+fl[11])+9.0*(fl[8]+fl[6]+fl[5])+15.58845726811989*fl[1])/(36.0*EPSILON+1.732050807568877*fl[14]+3.0*(fl[10]+fl[9]+fl[7])+5.196152422706631*(fl[4]+fl[3]+fl[2])+9.0*fl[0]); 
  fqVal[7] = 0.02777777777777778*(1.732050807568877*fl[14]+3.0*(fl[10]+fl[9]+fl[7])+5.196152422706631*(fl[4]+fl[3]+fl[2])+9.0*fl[0])*limTheta(rVal,1.0,cflL); 
  } else {
  rVal = (3.0*fr[15]+5.196152422706631*(fr[13]+fr[12]+fr[11])+9.0*(fr[8]+fr[6]+fr[5])+15.58845726811989*fr[1])/(36.0*EPSILON+1.732050807568877*fr[14]+3.0*(fr[10]+fr[9]+fr[7])+5.196152422706631*(fr[4]+fr[3]+fr[2])+9.0*fr[0]); 
  fqVal[7] = 0.02777777777777778*(1.732050807568877*fr[14]+3.0*(fr[10]+fr[9]+fr[7])+5.196152422706631*(fr[4]+fr[3]+fr[2])+9.0*fr[0])*limTheta(rVal,-1.0,cflR); 
  }
  double fhatALVal[8];  // fhatALVal = mode coefficients of anti-limited f 
  fhatALVal[0] = 0.3535533905932737*(fqVal[7]+fqVal[6]+fqVal[5]+fqVal[4]+fqVal[3]+fqVal[2]+fqVal[1]+fqVal[0]); 
  fhatALVal[1] = 0.6123724356957944*(fqVal[7]-1.0*fqVal[6]+fqVal[5]-1.0*fqVal[4]+fqVal[3]-1.0*fqVal[2]+fqVal[1]-1.0*fqVal[0]); 
  fhatALVal[2] = 0.6123724356957944*(fqVal[7]+fqVal[6]-1.0*(fqVal[5]+fqVal[4])+fqVal[3]+fqVal[2]-1.0*(fqVal[1]+fqVal[0])); 
  fhatALVal[3] = 0.6123724356957944*(fqVal[7]+fqVal[6]+fqVal[5]+fqVal[4]-1.0*(fqVal[3]+fqVal[2]+fqVal[1]+fqVal[0])); 
  fhatALVal[4] = 1.060660171779821*(fqVal[7]-1.0*(fqVal[6]+fqVal[5])+fqVal[4]+fqVal[3]-1.0*(fqVal[2]+fqVal[1])+fqVal[0]); 
  fhatALVal[5] = 1.060660171779821*(fqVal[7]-1.0*fqVal[6]+fqVal[5]-1.0*(fqVal[4]+fqVal[3])+fqVal[2]-1.0*fqVal[1]+fqVal[0]); 
  fhatALVal[6] = 1.060660171779821*(fqVal[7]+fqVal[6]-1.0*(fqVal[5]+fqVal[4]+fqVal[3]+fqVal[2])+fqVal[1]+fqVal[0]); 
  fhatALVal[7] = 1.837117307087383*(fqVal[7]-1.0*(fqVal[6]+fqVal[5])+fqVal[4]-1.0*fqVal[3]+fqVal[2]+fqVal[1]-1.0*fqVal[0]); 
  incr[0] = 0.25*(alpha[2]*fhatALVal[2]+alpha[1]*fhatALVal[1]+alpha[0]*fhatALVal[0])*dfac_x; 
  incr[1] = -0.4330127018922193*(alpha[2]*fhatALVal[2]+alpha[1]*fhatALVal[1]+alpha[0]*fhatALVal[0])*dfac_x; 
  incr[2] = 0.25*(alpha[2]*fhatALVal[4]+alpha[0]*fhatALVal[1]+fhatALVal[0]*alpha[1])*dfac_x; 
  incr[3] = 0.25*(alpha[1]*fhatALVal[4]+alpha[0]*fhatALVal[2]+fhatALVal[0]*alpha[2])*dfac_x; 
  incr[4] = 0.25*(alpha[2]*fhatALVal[6]+alpha[1]*fhatALVal[5]+alpha[0]*fhatALVal[3])*dfac_x; 
  incr[5] = -0.4330127018922193*(alpha[2]*fhatALVal[4]+alpha[0]*fhatALVal[1]+fhatALVal[0]*alpha[1])*dfac_x; 
  incr[6] = -0.4330127018922193*(alpha[1]*fhatALVal[4]+alpha[0]*fhatALVal[2]+fhatALVal[0]*alpha[2])*dfac_x; 
  incr[7] = 0.25*(alpha[0]*fhatALVal[4]+alpha[1]*fhatALVal[2]+fhatALVal[1]*alpha[2])*dfac_x; 
  incr[8] = -0.4330127018922193*(alpha[2]*fhatALVal[6]+alpha[1]*fhatALVal[5]+alpha[0]*fhatALVal[3])*dfac_x; 
  incr[9] = 0.25*(alpha[2]*fhatALVal[7]+alpha[0]*fhatALVal[5]+alpha[1]*fhatALVal[3])*dfac_x; 
  incr[10] = 0.25*(alpha[1]*fhatALVal[7]+alpha[0]*fhatALVal[6]+alpha[2]*fhatALVal[3])*dfac_x; 
  incr[11] = -0.4330127018922193*(alpha[0]*fhatALVal[4]+alpha[1]*fhatALVal[2]+fhatALVal[1]*alpha[2])*dfac_x; 
  incr[12] = -0.4330127018922193*(alpha[2]*fhatALVal[7]+alpha[0]*fhatALVal[5]+alpha[1]*fhatALVal[3])*dfac_x; 
  incr[13] = -0.4330127018922193*(alpha[1]*fhatALVal[7]+alpha[0]*fhatALVal[6]+alpha[2]*fhatALVal[3])*dfac_x; 
  incr[14] = 0.25*(alpha[0]*fhatALVal[7]+alpha[1]*fhatALVal[6]+alpha[2]*fhatALVal[5])*dfac_x; 
  incr[15] = -0.4330127018922193*(alpha[0]*fhatALVal[7]+alpha[1]*fhatALVal[6]+alpha[2]*fhatALVal[5])*dfac_x; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 
  outr[4] += incr[4]; 
  outr[5] += incr[5]; 
  outr[6] += incr[6]; 
  outr[7] += incr[7]; 
  outr[8] += incr[8]; 
  outr[9] += incr[9]; 
  outr[10] += incr[10]; 
  outr[11] += incr[11]; 
  outr[12] += incr[12]; 
  outr[13] += incr[13]; 
  outr[14] += incr[14]; 
  outr[15] += incr[15]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += incr[1]; 
  outl[2] += -1.0*incr[2]; 
  outl[3] += -1.0*incr[3]; 
  outl[4] += -1.0*incr[4]; 
  outl[5] += incr[5]; 
  outl[6] += incr[6]; 
  outl[7] += -1.0*incr[7]; 
  outl[8] += incr[8]; 
  outl[9] += -1.0*incr[9]; 
  outl[10] += -1.0*incr[10]; 
  outl[11] += incr[11]; 
  outl[12] += incr[12]; 
  outl[13] += incr[13]; 
  outl[14] += -1.0*incr[14]; 
  outl[15] += incr[15]; 
  return std::abs(alpha0); 
} 
double EmGyrokineticSurfPositivity2x2vSer_Y_P1_Bvars_1(const double q_, const double m_, const double cflL, const double cflR, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *Apar, const double *dApardt, const double *fl, const double *fr, double *outl, double *outr) 
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
  double incr[16]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = (0.125*(wv*(2.0*BdriftY[0]*m_*wv+(BmagInv[0]*(3.0*Apar[3]-1.732050807568877*Apar[1])*dfac_x-1.732050807568877*(BdriftY[1]*Apar[3]+BdriftY[0]*Apar[2])+Apar[1]*BdriftY[1]+Apar[0]*BdriftY[0])*q_)+BmagInv[0]*dfac_x*(1.732050807568877*Bmag[1]*wm+(1.732050807568877*Phi[1]-3.0*Phi[3])*q_)))/q_; 

  double alpha[8]; 
  alpha[0] = (0.7071067811865475*(2.0*BdriftY[0]*m_*wv2+(BmagInv[0]*(3.0*Apar[3]-1.732050807568877*Apar[1])*dfac_x-1.732050807568877*(BdriftY[1]*Apar[3]+BdriftY[0]*Apar[2])+Apar[1]*BdriftY[1]+Apar[0]*BdriftY[0])*q_*wv+BmagInv[0]*dfac_x*(1.732050807568877*Bmag[1]*wm+(1.732050807568877*Phi[1]-3.0*Phi[3])*q_)))/q_; 
  alpha[1] = (0.7071067811865475*(2.0*BdriftY[1]*m_*wv2+(BmagInv[1]*(3.0*Apar[3]-1.732050807568877*Apar[1])*dfac_x-1.732050807568877*BdriftY[0]*Apar[3]+BdriftY[1]*(Apar[0]-1.732050807568877*Apar[2])+BdriftY[0]*Apar[1])*q_*wv+BmagInv[1]*dfac_x*(1.732050807568877*Bmag[1]*wm+(1.732050807568877*Phi[1]-3.0*Phi[3])*q_)))/q_; 
  alpha[2] = (0.8164965809277261*BdriftY[0]*m_*wv)/(dfac_v*q_); 
  alpha[3] = (0.7071067811865475*BmagInv[0]*Bmag[1]*dfac_x)/(dfac_m*q_); 
  alpha[4] = (0.8164965809277261*BdriftY[1]*m_*wv)/(dfac_v*q_); 
  alpha[5] = (0.7071067811865475*Bmag[1]*BmagInv[1]*dfac_x)/(dfac_m*q_); 
  double rVal;  // rVal=f1/f0 at each control node in dimensions other than y 
  double fqVal[8];  // fqVal = anti-limited f evaluated at each control node on y surface 
  // determine upwinding at each surface control node 
  if(0.1178511301977579*(alpha[5]+alpha[4])-0.2041241452319315*(alpha[3]+alpha[2]+alpha[1])+0.3535533905932737*alpha[0] > 0) {
  rVal = -(1.0*(3.0*fl[15]-5.196152422706631*(fl[14]+fl[12]+fl[11])+9.0*(fl[9]+fl[7]+fl[5])-15.58845726811989*fl[2]))/(36.0*EPSILON-1.732050807568877*fl[13]+3.0*(fl[10]+fl[8]+fl[6])-5.196152422706631*(fl[4]+fl[3]+fl[1])+9.0*fl[0]); 
  fqVal[0] = -0.02777777777777778*(1.732050807568877*fl[13]-3.0*(fl[10]+fl[8]+fl[6])+5.196152422706631*(fl[4]+fl[3]+fl[1])-9.0*fl[0])*limTheta(rVal,1.0,cflL); 
  } else {
  rVal = -(1.0*(3.0*fr[15]-5.196152422706631*(fr[14]+fr[12]+fr[11])+9.0*(fr[9]+fr[7]+fr[5])-15.58845726811989*fr[2]))/(36.0*EPSILON-1.732050807568877*fr[13]+3.0*(fr[10]+fr[8]+fr[6])-5.196152422706631*(fr[4]+fr[3]+fr[1])+9.0*fr[0]); 
  fqVal[0] = -0.02777777777777778*(1.732050807568877*fr[13]-3.0*(fr[10]+fr[8]+fr[6])+5.196152422706631*(fr[4]+fr[3]+fr[1])-9.0*fr[0])*limTheta(rVal,-1.0,cflR); 
  }
  if((-0.1178511301977579*(alpha[5]+alpha[4]))-0.2041241452319315*(alpha[3]+alpha[2])+0.2041241452319315*alpha[1]+0.3535533905932737*alpha[0] > 0) {
  rVal = (3.0*fl[15]+5.196152422706631*fl[14]-5.196152422706631*(fl[12]+fl[11])-9.0*(fl[9]+fl[7])+9.0*fl[5]+15.58845726811989*fl[2])/(36.0*EPSILON+1.732050807568877*fl[13]+3.0*fl[10]-3.0*(fl[8]+fl[6])-5.196152422706631*(fl[4]+fl[3])+5.196152422706631*fl[1]+9.0*fl[0]); 
  fqVal[1] = 0.02777777777777778*(1.732050807568877*fl[13]+3.0*fl[10]-3.0*(fl[8]+fl[6])-5.196152422706631*(fl[4]+fl[3])+5.196152422706631*fl[1]+9.0*fl[0])*limTheta(rVal,1.0,cflL); 
  } else {
  rVal = (3.0*fr[15]+5.196152422706631*fr[14]-5.196152422706631*(fr[12]+fr[11])-9.0*(fr[9]+fr[7])+9.0*fr[5]+15.58845726811989*fr[2])/(36.0*EPSILON+1.732050807568877*fr[13]+3.0*fr[10]-3.0*(fr[8]+fr[6])-5.196152422706631*(fr[4]+fr[3])+5.196152422706631*fr[1]+9.0*fr[0]); 
  fqVal[1] = 0.02777777777777778*(1.732050807568877*fr[13]+3.0*fr[10]-3.0*(fr[8]+fr[6])-5.196152422706631*(fr[4]+fr[3])+5.196152422706631*fr[1]+9.0*fr[0])*limTheta(rVal,-1.0,cflR); 
  }
  if(0.1178511301977579*alpha[5]-0.1178511301977579*alpha[4]-0.2041241452319315*alpha[3]+0.2041241452319315*alpha[2]-0.2041241452319315*alpha[1]+0.3535533905932737*alpha[0] > 0) {
  rVal = (3.0*fl[15]-5.196152422706631*fl[14]+5.196152422706631*fl[12]-5.196152422706631*fl[11]-9.0*fl[9]+9.0*fl[7]-9.0*fl[5]+15.58845726811989*fl[2])/(36.0*EPSILON+1.732050807568877*fl[13]-3.0*fl[10]+3.0*fl[8]-3.0*fl[6]-5.196152422706631*fl[4]+5.196152422706631*fl[3]-5.196152422706631*fl[1]+9.0*fl[0]); 
  fqVal[2] = 0.02777777777777778*(1.732050807568877*fl[13]-3.0*fl[10]+3.0*fl[8]-3.0*fl[6]-5.196152422706631*fl[4]+5.196152422706631*fl[3]-5.196152422706631*fl[1]+9.0*fl[0])*limTheta(rVal,1.0,cflL); 
  } else {
  rVal = (3.0*fr[15]-5.196152422706631*fr[14]+5.196152422706631*fr[12]-5.196152422706631*fr[11]-9.0*fr[9]+9.0*fr[7]-9.0*fr[5]+15.58845726811989*fr[2])/(36.0*EPSILON+1.732050807568877*fr[13]-3.0*fr[10]+3.0*fr[8]-3.0*fr[6]-5.196152422706631*fr[4]+5.196152422706631*fr[3]-5.196152422706631*fr[1]+9.0*fr[0]); 
  fqVal[2] = 0.02777777777777778*(1.732050807568877*fr[13]-3.0*fr[10]+3.0*fr[8]-3.0*fr[6]-5.196152422706631*fr[4]+5.196152422706631*fr[3]-5.196152422706631*fr[1]+9.0*fr[0])*limTheta(rVal,-1.0,cflR); 
  }
  if((-0.1178511301977579*alpha[5])+0.1178511301977579*alpha[4]-0.2041241452319315*alpha[3]+0.2041241452319315*(alpha[2]+alpha[1])+0.3535533905932737*alpha[0] > 0) {
  rVal = -(1.0*(3.0*fl[15]+5.196152422706631*(fl[14]+fl[12])-5.196152422706631*fl[11]+9.0*fl[9]-9.0*(fl[7]+fl[5])-15.58845726811989*fl[2]))/(36.0*EPSILON-1.732050807568877*fl[13]-3.0*(fl[10]+fl[8])+3.0*fl[6]-5.196152422706631*fl[4]+5.196152422706631*(fl[3]+fl[1])+9.0*fl[0]); 
  fqVal[3] = -0.02777777777777778*(1.732050807568877*fl[13]+3.0*(fl[10]+fl[8])-3.0*fl[6]+5.196152422706631*fl[4]-5.196152422706631*(fl[3]+fl[1])-9.0*fl[0])*limTheta(rVal,1.0,cflL); 
  } else {
  rVal = -(1.0*(3.0*fr[15]+5.196152422706631*(fr[14]+fr[12])-5.196152422706631*fr[11]+9.0*fr[9]-9.0*(fr[7]+fr[5])-15.58845726811989*fr[2]))/(36.0*EPSILON-1.732050807568877*fr[13]-3.0*(fr[10]+fr[8])+3.0*fr[6]-5.196152422706631*fr[4]+5.196152422706631*(fr[3]+fr[1])+9.0*fr[0]); 
  fqVal[3] = -0.02777777777777778*(1.732050807568877*fr[13]+3.0*(fr[10]+fr[8])-3.0*fr[6]+5.196152422706631*fr[4]-5.196152422706631*(fr[3]+fr[1])-9.0*fr[0])*limTheta(rVal,-1.0,cflR); 
  }
  if((-0.1178511301977579*alpha[5])+0.1178511301977579*alpha[4]+0.2041241452319315*alpha[3]-0.2041241452319315*(alpha[2]+alpha[1])+0.3535533905932737*alpha[0] > 0) {
  rVal = (3.0*fl[15]-5.196152422706631*(fl[14]+fl[12])+5.196152422706631*fl[11]+9.0*fl[9]-9.0*(fl[7]+fl[5])+15.58845726811989*fl[2])/(36.0*EPSILON+1.732050807568877*fl[13]-3.0*(fl[10]+fl[8])+3.0*fl[6]+5.196152422706631*fl[4]-5.196152422706631*(fl[3]+fl[1])+9.0*fl[0]); 
  fqVal[4] = 0.02777777777777778*(1.732050807568877*fl[13]-3.0*(fl[10]+fl[8])+3.0*fl[6]+5.196152422706631*fl[4]-5.196152422706631*(fl[3]+fl[1])+9.0*fl[0])*limTheta(rVal,1.0,cflL); 
  } else {
  rVal = (3.0*fr[15]-5.196152422706631*(fr[14]+fr[12])+5.196152422706631*fr[11]+9.0*fr[9]-9.0*(fr[7]+fr[5])+15.58845726811989*fr[2])/(36.0*EPSILON+1.732050807568877*fr[13]-3.0*(fr[10]+fr[8])+3.0*fr[6]+5.196152422706631*fr[4]-5.196152422706631*(fr[3]+fr[1])+9.0*fr[0]); 
  fqVal[4] = 0.02777777777777778*(1.732050807568877*fr[13]-3.0*(fr[10]+fr[8])+3.0*fr[6]+5.196152422706631*fr[4]-5.196152422706631*(fr[3]+fr[1])+9.0*fr[0])*limTheta(rVal,-1.0,cflR); 
  }
  if(0.1178511301977579*alpha[5]-0.1178511301977579*alpha[4]+0.2041241452319315*alpha[3]-0.2041241452319315*alpha[2]+0.2041241452319315*alpha[1]+0.3535533905932737*alpha[0] > 0) {
  rVal = -(1.0*(3.0*fl[15]+5.196152422706631*fl[14]-5.196152422706631*fl[12]+5.196152422706631*fl[11]-9.0*fl[9]+9.0*fl[7]-9.0*fl[5]-15.58845726811989*fl[2]))/(36.0*EPSILON-1.732050807568877*fl[13]-3.0*fl[10]+3.0*fl[8]-3.0*fl[6]+5.196152422706631*fl[4]-5.196152422706631*fl[3]+5.196152422706631*fl[1]+9.0*fl[0]); 
  fqVal[5] = -0.02777777777777778*(1.732050807568877*fl[13]+3.0*fl[10]-3.0*fl[8]+3.0*fl[6]-5.196152422706631*fl[4]+5.196152422706631*fl[3]-5.196152422706631*fl[1]-9.0*fl[0])*limTheta(rVal,1.0,cflL); 
  } else {
  rVal = -(1.0*(3.0*fr[15]+5.196152422706631*fr[14]-5.196152422706631*fr[12]+5.196152422706631*fr[11]-9.0*fr[9]+9.0*fr[7]-9.0*fr[5]-15.58845726811989*fr[2]))/(36.0*EPSILON-1.732050807568877*fr[13]-3.0*fr[10]+3.0*fr[8]-3.0*fr[6]+5.196152422706631*fr[4]-5.196152422706631*fr[3]+5.196152422706631*fr[1]+9.0*fr[0]); 
  fqVal[5] = -0.02777777777777778*(1.732050807568877*fr[13]+3.0*fr[10]-3.0*fr[8]+3.0*fr[6]-5.196152422706631*fr[4]+5.196152422706631*fr[3]-5.196152422706631*fr[1]-9.0*fr[0])*limTheta(rVal,-1.0,cflR); 
  }
  if((-0.1178511301977579*(alpha[5]+alpha[4]))+0.2041241452319315*(alpha[3]+alpha[2])-0.2041241452319315*alpha[1]+0.3535533905932737*alpha[0] > 0) {
  rVal = -(1.0*(3.0*fl[15]-5.196152422706631*fl[14]+5.196152422706631*(fl[12]+fl[11])-9.0*(fl[9]+fl[7])+9.0*fl[5]-15.58845726811989*fl[2]))/(36.0*EPSILON-1.732050807568877*fl[13]+3.0*fl[10]-3.0*(fl[8]+fl[6])+5.196152422706631*(fl[4]+fl[3])-5.196152422706631*fl[1]+9.0*fl[0]); 
  fqVal[6] = -0.02777777777777778*(1.732050807568877*fl[13]-3.0*fl[10]+3.0*(fl[8]+fl[6])-5.196152422706631*(fl[4]+fl[3])+5.196152422706631*fl[1]-9.0*fl[0])*limTheta(rVal,1.0,cflL); 
  } else {
  rVal = -(1.0*(3.0*fr[15]-5.196152422706631*fr[14]+5.196152422706631*(fr[12]+fr[11])-9.0*(fr[9]+fr[7])+9.0*fr[5]-15.58845726811989*fr[2]))/(36.0*EPSILON-1.732050807568877*fr[13]+3.0*fr[10]-3.0*(fr[8]+fr[6])+5.196152422706631*(fr[4]+fr[3])-5.196152422706631*fr[1]+9.0*fr[0]); 
  fqVal[6] = -0.02777777777777778*(1.732050807568877*fr[13]-3.0*fr[10]+3.0*(fr[8]+fr[6])-5.196152422706631*(fr[4]+fr[3])+5.196152422706631*fr[1]-9.0*fr[0])*limTheta(rVal,-1.0,cflR); 
  }
  if(0.1178511301977579*(alpha[5]+alpha[4])+0.2041241452319315*(alpha[3]+alpha[2]+alpha[1])+0.3535533905932737*alpha[0] > 0) {
  rVal = (3.0*fl[15]+5.196152422706631*(fl[14]+fl[12]+fl[11])+9.0*(fl[9]+fl[7]+fl[5])+15.58845726811989*fl[2])/(36.0*EPSILON+1.732050807568877*fl[13]+3.0*(fl[10]+fl[8]+fl[6])+5.196152422706631*(fl[4]+fl[3]+fl[1])+9.0*fl[0]); 
  fqVal[7] = 0.02777777777777778*(1.732050807568877*fl[13]+3.0*(fl[10]+fl[8]+fl[6])+5.196152422706631*(fl[4]+fl[3]+fl[1])+9.0*fl[0])*limTheta(rVal,1.0,cflL); 
  } else {
  rVal = (3.0*fr[15]+5.196152422706631*(fr[14]+fr[12]+fr[11])+9.0*(fr[9]+fr[7]+fr[5])+15.58845726811989*fr[2])/(36.0*EPSILON+1.732050807568877*fr[13]+3.0*(fr[10]+fr[8]+fr[6])+5.196152422706631*(fr[4]+fr[3]+fr[1])+9.0*fr[0]); 
  fqVal[7] = 0.02777777777777778*(1.732050807568877*fr[13]+3.0*(fr[10]+fr[8]+fr[6])+5.196152422706631*(fr[4]+fr[3]+fr[1])+9.0*fr[0])*limTheta(rVal,-1.0,cflR); 
  }
  double fhatALVal[8];  // fhatALVal = mode coefficients of anti-limited f 
  fhatALVal[0] = 0.3535533905932737*(fqVal[7]+fqVal[6]+fqVal[5]+fqVal[4]+fqVal[3]+fqVal[2]+fqVal[1]+fqVal[0]); 
  fhatALVal[1] = 0.6123724356957944*(fqVal[7]-1.0*fqVal[6]+fqVal[5]-1.0*fqVal[4]+fqVal[3]-1.0*fqVal[2]+fqVal[1]-1.0*fqVal[0]); 
  fhatALVal[2] = 0.6123724356957944*(fqVal[7]+fqVal[6]-1.0*(fqVal[5]+fqVal[4])+fqVal[3]+fqVal[2]-1.0*(fqVal[1]+fqVal[0])); 
  fhatALVal[3] = 0.6123724356957944*(fqVal[7]+fqVal[6]+fqVal[5]+fqVal[4]-1.0*(fqVal[3]+fqVal[2]+fqVal[1]+fqVal[0])); 
  fhatALVal[4] = 1.060660171779821*(fqVal[7]-1.0*(fqVal[6]+fqVal[5])+fqVal[4]+fqVal[3]-1.0*(fqVal[2]+fqVal[1])+fqVal[0]); 
  fhatALVal[5] = 1.060660171779821*(fqVal[7]-1.0*fqVal[6]+fqVal[5]-1.0*(fqVal[4]+fqVal[3])+fqVal[2]-1.0*fqVal[1]+fqVal[0]); 
  fhatALVal[6] = 1.060660171779821*(fqVal[7]+fqVal[6]-1.0*(fqVal[5]+fqVal[4]+fqVal[3]+fqVal[2])+fqVal[1]+fqVal[0]); 
  fhatALVal[7] = 1.837117307087383*(fqVal[7]-1.0*(fqVal[6]+fqVal[5])+fqVal[4]-1.0*fqVal[3]+fqVal[2]+fqVal[1]-1.0*fqVal[0]); 
  incr[0] = 0.25*(alpha[5]*fhatALVal[5]+alpha[4]*fhatALVal[4]+alpha[3]*fhatALVal[3]+alpha[2]*fhatALVal[2]+alpha[1]*fhatALVal[1]+alpha[0]*fhatALVal[0])*dfac_y; 
  incr[1] = 0.25*(alpha[3]*fhatALVal[5]+fhatALVal[3]*alpha[5]+alpha[2]*fhatALVal[4]+fhatALVal[2]*alpha[4]+alpha[0]*fhatALVal[1]+fhatALVal[0]*alpha[1])*dfac_y; 
  incr[2] = -0.4330127018922193*(alpha[5]*fhatALVal[5]+alpha[4]*fhatALVal[4]+alpha[3]*fhatALVal[3]+alpha[2]*fhatALVal[2]+alpha[1]*fhatALVal[1]+alpha[0]*fhatALVal[0])*dfac_y; 
  incr[3] = 0.25*(alpha[5]*fhatALVal[7]+alpha[3]*fhatALVal[6]+alpha[1]*fhatALVal[4]+fhatALVal[1]*alpha[4]+alpha[0]*fhatALVal[2]+fhatALVal[0]*alpha[2])*dfac_y; 
  incr[4] = 0.25*(alpha[4]*fhatALVal[7]+alpha[2]*fhatALVal[6]+alpha[1]*fhatALVal[5]+fhatALVal[1]*alpha[5]+alpha[0]*fhatALVal[3]+fhatALVal[0]*alpha[3])*dfac_y; 
  incr[5] = -0.4330127018922193*(alpha[3]*fhatALVal[5]+fhatALVal[3]*alpha[5]+alpha[2]*fhatALVal[4]+fhatALVal[2]*alpha[4]+alpha[0]*fhatALVal[1]+fhatALVal[0]*alpha[1])*dfac_y; 
  incr[6] = 0.25*(alpha[3]*fhatALVal[7]+alpha[5]*fhatALVal[6]+alpha[0]*fhatALVal[4]+fhatALVal[0]*alpha[4]+alpha[1]*fhatALVal[2]+fhatALVal[1]*alpha[2])*dfac_y; 
  incr[7] = -0.4330127018922193*(alpha[5]*fhatALVal[7]+alpha[3]*fhatALVal[6]+alpha[1]*fhatALVal[4]+fhatALVal[1]*alpha[4]+alpha[0]*fhatALVal[2]+fhatALVal[0]*alpha[2])*dfac_y; 
  incr[8] = 0.25*(alpha[2]*fhatALVal[7]+alpha[4]*fhatALVal[6]+alpha[0]*fhatALVal[5]+fhatALVal[0]*alpha[5]+alpha[1]*fhatALVal[3]+fhatALVal[1]*alpha[3])*dfac_y; 
  incr[9] = -0.4330127018922193*(alpha[4]*fhatALVal[7]+alpha[2]*fhatALVal[6]+alpha[1]*fhatALVal[5]+fhatALVal[1]*alpha[5]+alpha[0]*fhatALVal[3]+fhatALVal[0]*alpha[3])*dfac_y; 
  incr[10] = 0.25*(alpha[1]*fhatALVal[7]+alpha[0]*fhatALVal[6]+alpha[4]*fhatALVal[5]+fhatALVal[4]*alpha[5]+alpha[2]*fhatALVal[3]+fhatALVal[2]*alpha[3])*dfac_y; 
  incr[11] = -0.4330127018922193*(alpha[3]*fhatALVal[7]+alpha[5]*fhatALVal[6]+alpha[0]*fhatALVal[4]+fhatALVal[0]*alpha[4]+alpha[1]*fhatALVal[2]+fhatALVal[1]*alpha[2])*dfac_y; 
  incr[12] = -0.4330127018922193*(alpha[2]*fhatALVal[7]+alpha[4]*fhatALVal[6]+alpha[0]*fhatALVal[5]+fhatALVal[0]*alpha[5]+alpha[1]*fhatALVal[3]+fhatALVal[1]*alpha[3])*dfac_y; 
  incr[13] = 0.25*(alpha[0]*fhatALVal[7]+alpha[1]*fhatALVal[6]+alpha[2]*fhatALVal[5]+fhatALVal[2]*alpha[5]+alpha[3]*fhatALVal[4]+fhatALVal[3]*alpha[4])*dfac_y; 
  incr[14] = -0.4330127018922193*(alpha[1]*fhatALVal[7]+alpha[0]*fhatALVal[6]+alpha[4]*fhatALVal[5]+fhatALVal[4]*alpha[5]+alpha[2]*fhatALVal[3]+fhatALVal[2]*alpha[3])*dfac_y; 
  incr[15] = -0.4330127018922193*(alpha[0]*fhatALVal[7]+alpha[1]*fhatALVal[6]+alpha[2]*fhatALVal[5]+fhatALVal[2]*alpha[5]+alpha[3]*fhatALVal[4]+fhatALVal[3]*alpha[4])*dfac_y; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 
  outr[4] += incr[4]; 
  outr[5] += incr[5]; 
  outr[6] += incr[6]; 
  outr[7] += incr[7]; 
  outr[8] += incr[8]; 
  outr[9] += incr[9]; 
  outr[10] += incr[10]; 
  outr[11] += incr[11]; 
  outr[12] += incr[12]; 
  outr[13] += incr[13]; 
  outr[14] += incr[14]; 
  outr[15] += incr[15]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += -1.0*incr[1]; 
  outl[2] += incr[2]; 
  outl[3] += -1.0*incr[3]; 
  outl[4] += -1.0*incr[4]; 
  outl[5] += incr[5]; 
  outl[6] += -1.0*incr[6]; 
  outl[7] += incr[7]; 
  outl[8] += -1.0*incr[8]; 
  outl[9] += incr[9]; 
  outl[10] += -1.0*incr[10]; 
  outl[11] += incr[11]; 
  outl[12] += incr[12]; 
  outl[13] += -1.0*incr[13]; 
  outl[14] += incr[14]; 
  outl[15] += incr[15]; 
  return std::abs(alpha0); 
} 
double EmGyrokineticSurfPositivity2x2vSer_Vpar_P1_Bvars_1(const double q_, const double m_, const double cflL, const double cflR, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *Apar, const double *dApardt, const double *fl, const double *fr, double *outl, double *outr) 
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
  double incr[16]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = -(0.0625*(3.464101615137754*dfac_v*m_*(BdriftX[0]*Bmag[1]*dfac_x*wm+((BdriftY[1]*Phi[3]+BdriftY[0]*Phi[2])*dfac_y+BdriftX[0]*Phi[1]*dfac_x)*q_)*wv+Bmag[1]*dfac_x*(dfac_v*(3.0*(BmagInv[1]*Apar[3]+BmagInv[0]*Apar[2])*dfac_y+1.732050807568877*(Apar[1]*BdriftX[1]+Apar[0]*BdriftX[0]))*q_-3.464101615137754*BdriftX[0]*m_)*wm+q_*(dfac_v*((1.732050807568877*((Apar[0]*BdriftY[1]+BdriftY[0]*Apar[1])*Phi[3]+(Apar[1]*BdriftY[1]+Apar[0]*BdriftY[0])*Phi[2])-3.0*(BmagInv[1]*(Apar[1]*Phi[3]-1.0*Phi[1]*Apar[3])+BmagInv[0]*(Apar[1]*Phi[2]-1.0*Phi[1]*Apar[2]))*dfac_x)*dfac_y+1.732050807568877*((BdriftX[1]*Apar[3]+BdriftX[0]*Apar[2])*Phi[3]+(Apar[1]*BdriftX[1]+Apar[0]*BdriftX[0])*Phi[1])*dfac_x+4.0*dApardt[0])*q_-3.464101615137754*((BdriftY[1]*Phi[3]+BdriftY[0]*Phi[2])*dfac_y+BdriftX[0]*Phi[1]*dfac_x)*m_)))/(dfac_v*m_*q_); 

  double alpha[8]; 
  alpha[0] = -(0.3535533905932737*(3.464101615137754*dfac_v*m_*(BdriftX[0]*Bmag[1]*dfac_x*wm+((BdriftY[1]*Phi[3]+BdriftY[0]*Phi[2])*dfac_y+BdriftX[0]*Phi[1]*dfac_x)*q_)*wv+Bmag[1]*dfac_x*(dfac_v*(3.0*(BmagInv[1]*Apar[3]+BmagInv[0]*Apar[2])*dfac_y+1.732050807568877*(Apar[1]*BdriftX[1]+Apar[0]*BdriftX[0]))*q_-3.464101615137754*BdriftX[0]*m_)*wm-3.464101615137754*((BdriftY[1]*Phi[3]+BdriftY[0]*Phi[2])*dfac_y+BdriftX[0]*Phi[1]*dfac_x)*m_*q_+dfac_v*(((BmagInv[1]*(3.0*Phi[1]*Apar[3]-3.0*Apar[1]*Phi[3])+BmagInv[0]*(3.0*Phi[1]*Apar[2]-3.0*Apar[1]*Phi[2]))*dfac_x+1.732050807568877*((Apar[0]*BdriftY[1]+BdriftY[0]*Apar[1])*Phi[3]+(Apar[1]*BdriftY[1]+Apar[0]*BdriftY[0])*Phi[2]))*dfac_y+1.732050807568877*((BdriftX[1]*Apar[3]+BdriftX[0]*Apar[2])*Phi[3]+(Apar[1]*BdriftX[1]+Apar[0]*BdriftX[0])*Phi[1])*dfac_x+4.0*dApardt[0])*q2))/(dfac_v*m_*q_); 
  alpha[1] = -(0.07071067811865474*(17.32050807568877*dfac_v*m_*(BdriftX[1]*Bmag[1]*dfac_x*wm+((BdriftY[0]*Phi[3]+BdriftY[1]*Phi[2])*dfac_y+BdriftX[1]*Phi[1]*dfac_x)*q_)*wv+Bmag[1]*dfac_x*(dfac_v*(15.0*(BmagInv[0]*Apar[3]+BmagInv[1]*Apar[2])*dfac_y+8.660254037844386*(Apar[0]*BdriftX[1]+BdriftX[0]*Apar[1]))*q_-17.32050807568877*BdriftX[1]*m_)*wm-17.32050807568877*((BdriftY[0]*Phi[3]+BdriftY[1]*Phi[2])*dfac_y+BdriftX[1]*Phi[1]*dfac_x)*m_*q_+dfac_v*(((BmagInv[0]*(15.0*Phi[1]*Apar[3]-15.0*Apar[1]*Phi[3])+BmagInv[1]*(15.0*Phi[1]*Apar[2]-15.0*Apar[1]*Phi[2]))*dfac_x+(15.58845726811989*Apar[1]*BdriftY[1]+8.660254037844386*Apar[0]*BdriftY[0])*Phi[3]+8.660254037844386*(Apar[0]*BdriftY[1]+BdriftY[0]*Apar[1])*Phi[2])*dfac_y+8.660254037844386*((BdriftX[0]*Apar[3]+BdriftX[1]*Apar[2])*Phi[3]+(Apar[0]*BdriftX[1]+BdriftX[0]*Apar[1])*Phi[1])*dfac_x+20.0*dApardt[1])*q2))/(dfac_v*m_*q_); 
  alpha[2] = -(0.3535533905932737*(dfac_v*(dfac_x*(3.464101615137754*BdriftX[0]*Phi[3]*m_*wv+1.732050807568877*Bmag[1]*(BdriftX[1]*Apar[3]+BdriftX[0]*Apar[2])*wm)+((BmagInv[0]*(3.0*Apar[2]*Phi[3]-3.0*Phi[2]*Apar[3])*dfac_x+1.732050807568877*((BdriftY[0]*Apar[3]+BdriftY[1]*Apar[2])*Phi[3]+Phi[2]*(BdriftY[1]*Apar[3]+BdriftY[0]*Apar[2])))*dfac_y+1.732050807568877*((Apar[1]*BdriftX[1]+Apar[0]*BdriftX[0])*Phi[3]+Phi[1]*(BdriftX[1]*Apar[3]+BdriftX[0]*Apar[2]))*dfac_x+4.0*dApardt[2])*q_)-3.464101615137754*BdriftX[0]*Phi[3]*dfac_x*m_))/(dfac_v*m_); 
  alpha[3] = -(0.3535533905932737*Bmag[1]*dfac_x*(dfac_v*(2.0*BdriftX[0]*m_*wv+(1.732050807568877*(BmagInv[1]*Apar[3]+BmagInv[0]*Apar[2])*dfac_y+Apar[1]*BdriftX[1]+Apar[0]*BdriftX[0])*q_)-2.0*BdriftX[0]*m_))/(dfac_m*dfac_v*m_*q_); 
  alpha[4] = -(0.07071067811865474*(dfac_v*(dfac_x*(17.32050807568877*BdriftX[1]*Phi[3]*m_*wv+8.660254037844386*Bmag[1]*(BdriftX[0]*Apar[3]+BdriftX[1]*Apar[2])*wm)+((BmagInv[1]*(15.0*Apar[2]*Phi[3]-15.0*Phi[2]*Apar[3])*dfac_x+(15.58845726811989*BdriftY[1]*Apar[3]+8.660254037844386*BdriftY[0]*Apar[2])*Phi[3]+8.660254037844386*Phi[2]*(BdriftY[0]*Apar[3]+BdriftY[1]*Apar[2]))*dfac_y+8.660254037844386*((Apar[0]*BdriftX[1]+BdriftX[0]*Apar[1])*Phi[3]+Phi[1]*(BdriftX[0]*Apar[3]+BdriftX[1]*Apar[2]))*dfac_x+20.0*dApardt[3])*q_)-17.32050807568877*BdriftX[1]*Phi[3]*dfac_x*m_))/(dfac_v*m_); 
  alpha[5] = -(0.3535533905932737*Bmag[1]*dfac_x*(dfac_v*(2.0*BdriftX[1]*m_*wv+(1.732050807568877*(BmagInv[0]*Apar[3]+BmagInv[1]*Apar[2])*dfac_y+Apar[0]*BdriftX[1]+BdriftX[0]*Apar[1])*q_)-2.0*BdriftX[1]*m_))/(dfac_m*dfac_v*m_*q_); 
  alpha[6] = -(0.3535533905932737*Bmag[1]*(BdriftX[1]*Apar[3]+BdriftX[0]*Apar[2])*dfac_x)/(dfac_m*m_); 
  alpha[7] = -(0.3535533905932737*Bmag[1]*(BdriftX[0]*Apar[3]+BdriftX[1]*Apar[2])*dfac_x)/(dfac_m*m_); 
  double rVal;  // rVal=f1/f0 at each control node in dimensions other than vx 
  double fqVal[8];  // fqVal = anti-limited f evaluated at each control node on vx surface 
  // determine upwinding at each surface control node 
  if((-0.06804138174397717*alpha[7])+0.1178511301977579*(alpha[6]+alpha[5]+alpha[4])-0.2041241452319315*(alpha[3]+alpha[2]+alpha[1])+0.3535533905932737*alpha[0] > 0) {
  rVal = -(1.0*(3.0*fl[15]-5.196152422706631*(fl[14]+fl[13]+fl[11])+9.0*(fl[10]+fl[7]+fl[6])-15.58845726811989*fl[3]))/(36.0*EPSILON-1.732050807568877*fl[12]+3.0*(fl[9]+fl[8]+fl[5])-5.196152422706631*(fl[4]+fl[2]+fl[1])+9.0*fl[0]); 
  fqVal[0] = -0.02777777777777778*(1.732050807568877*fl[12]-3.0*(fl[9]+fl[8]+fl[5])+5.196152422706631*(fl[4]+fl[2]+fl[1])-9.0*fl[0])*limTheta(rVal,1.0,cflL); 
  } else {
  rVal = -(1.0*(3.0*fr[15]-5.196152422706631*(fr[14]+fr[13]+fr[11])+9.0*(fr[10]+fr[7]+fr[6])-15.58845726811989*fr[3]))/(36.0*EPSILON-1.732050807568877*fr[12]+3.0*(fr[9]+fr[8]+fr[5])-5.196152422706631*(fr[4]+fr[2]+fr[1])+9.0*fr[0]); 
  fqVal[0] = -0.02777777777777778*(1.732050807568877*fr[12]-3.0*(fr[9]+fr[8]+fr[5])+5.196152422706631*(fr[4]+fr[2]+fr[1])-9.0*fr[0])*limTheta(rVal,-1.0,cflR); 
  }
  if(0.06804138174397717*alpha[7]+0.1178511301977579*alpha[6]-0.1178511301977579*(alpha[5]+alpha[4])-0.2041241452319315*(alpha[3]+alpha[2])+0.2041241452319315*alpha[1]+0.3535533905932737*alpha[0] > 0) {
  rVal = (3.0*fl[15]+5.196152422706631*fl[14]-5.196152422706631*(fl[13]+fl[11])-9.0*(fl[10]+fl[7])+9.0*fl[6]+15.58845726811989*fl[3])/(36.0*EPSILON+1.732050807568877*fl[12]+3.0*fl[9]-3.0*(fl[8]+fl[5])-5.196152422706631*(fl[4]+fl[2])+5.196152422706631*fl[1]+9.0*fl[0]); 
  fqVal[1] = 0.02777777777777778*(1.732050807568877*fl[12]+3.0*fl[9]-3.0*(fl[8]+fl[5])-5.196152422706631*(fl[4]+fl[2])+5.196152422706631*fl[1]+9.0*fl[0])*limTheta(rVal,1.0,cflL); 
  } else {
  rVal = (3.0*fr[15]+5.196152422706631*fr[14]-5.196152422706631*(fr[13]+fr[11])-9.0*(fr[10]+fr[7])+9.0*fr[6]+15.58845726811989*fr[3])/(36.0*EPSILON+1.732050807568877*fr[12]+3.0*fr[9]-3.0*(fr[8]+fr[5])-5.196152422706631*(fr[4]+fr[2])+5.196152422706631*fr[1]+9.0*fr[0]); 
  fqVal[1] = 0.02777777777777778*(1.732050807568877*fr[12]+3.0*fr[9]-3.0*(fr[8]+fr[5])-5.196152422706631*(fr[4]+fr[2])+5.196152422706631*fr[1]+9.0*fr[0])*limTheta(rVal,-1.0,cflR); 
  }
  if(0.06804138174397717*alpha[7]-0.1178511301977579*alpha[6]+0.1178511301977579*alpha[5]-0.1178511301977579*alpha[4]-0.2041241452319315*alpha[3]+0.2041241452319315*alpha[2]-0.2041241452319315*alpha[1]+0.3535533905932737*alpha[0] > 0) {
  rVal = (3.0*fl[15]-5.196152422706631*fl[14]+5.196152422706631*fl[13]-5.196152422706631*fl[11]-9.0*fl[10]+9.0*fl[7]-9.0*fl[6]+15.58845726811989*fl[3])/(36.0*EPSILON+1.732050807568877*fl[12]-3.0*fl[9]+3.0*fl[8]-3.0*fl[5]-5.196152422706631*fl[4]+5.196152422706631*fl[2]-5.196152422706631*fl[1]+9.0*fl[0]); 
  fqVal[2] = 0.02777777777777778*(1.732050807568877*fl[12]-3.0*fl[9]+3.0*fl[8]-3.0*fl[5]-5.196152422706631*fl[4]+5.196152422706631*fl[2]-5.196152422706631*fl[1]+9.0*fl[0])*limTheta(rVal,1.0,cflL); 
  } else {
  rVal = (3.0*fr[15]-5.196152422706631*fr[14]+5.196152422706631*fr[13]-5.196152422706631*fr[11]-9.0*fr[10]+9.0*fr[7]-9.0*fr[6]+15.58845726811989*fr[3])/(36.0*EPSILON+1.732050807568877*fr[12]-3.0*fr[9]+3.0*fr[8]-3.0*fr[5]-5.196152422706631*fr[4]+5.196152422706631*fr[2]-5.196152422706631*fr[1]+9.0*fr[0]); 
  fqVal[2] = 0.02777777777777778*(1.732050807568877*fr[12]-3.0*fr[9]+3.0*fr[8]-3.0*fr[5]-5.196152422706631*fr[4]+5.196152422706631*fr[2]-5.196152422706631*fr[1]+9.0*fr[0])*limTheta(rVal,-1.0,cflR); 
  }
  if((-0.06804138174397717*alpha[7])-0.1178511301977579*(alpha[6]+alpha[5])+0.1178511301977579*alpha[4]-0.2041241452319315*alpha[3]+0.2041241452319315*(alpha[2]+alpha[1])+0.3535533905932737*alpha[0] > 0) {
  rVal = -(1.0*(3.0*fl[15]+5.196152422706631*(fl[14]+fl[13])-5.196152422706631*fl[11]+9.0*fl[10]-9.0*(fl[7]+fl[6])-15.58845726811989*fl[3]))/(36.0*EPSILON-1.732050807568877*fl[12]-3.0*(fl[9]+fl[8])+3.0*fl[5]-5.196152422706631*fl[4]+5.196152422706631*(fl[2]+fl[1])+9.0*fl[0]); 
  fqVal[3] = -0.02777777777777778*(1.732050807568877*fl[12]+3.0*(fl[9]+fl[8])-3.0*fl[5]+5.196152422706631*fl[4]-5.196152422706631*(fl[2]+fl[1])-9.0*fl[0])*limTheta(rVal,1.0,cflL); 
  } else {
  rVal = -(1.0*(3.0*fr[15]+5.196152422706631*(fr[14]+fr[13])-5.196152422706631*fr[11]+9.0*fr[10]-9.0*(fr[7]+fr[6])-15.58845726811989*fr[3]))/(36.0*EPSILON-1.732050807568877*fr[12]-3.0*(fr[9]+fr[8])+3.0*fr[5]-5.196152422706631*fr[4]+5.196152422706631*(fr[2]+fr[1])+9.0*fr[0]); 
  fqVal[3] = -0.02777777777777778*(1.732050807568877*fr[12]+3.0*(fr[9]+fr[8])-3.0*fr[5]+5.196152422706631*fr[4]-5.196152422706631*(fr[2]+fr[1])-9.0*fr[0])*limTheta(rVal,-1.0,cflR); 
  }
  if(0.06804138174397717*alpha[7]-0.1178511301977579*(alpha[6]+alpha[5])+0.1178511301977579*alpha[4]+0.2041241452319315*alpha[3]-0.2041241452319315*(alpha[2]+alpha[1])+0.3535533905932737*alpha[0] > 0) {
  rVal = (3.0*fl[15]-5.196152422706631*(fl[14]+fl[13])+5.196152422706631*fl[11]+9.0*fl[10]-9.0*(fl[7]+fl[6])+15.58845726811989*fl[3])/(36.0*EPSILON+1.732050807568877*fl[12]-3.0*(fl[9]+fl[8])+3.0*fl[5]+5.196152422706631*fl[4]-5.196152422706631*(fl[2]+fl[1])+9.0*fl[0]); 
  fqVal[4] = 0.02777777777777778*(1.732050807568877*fl[12]-3.0*(fl[9]+fl[8])+3.0*fl[5]+5.196152422706631*fl[4]-5.196152422706631*(fl[2]+fl[1])+9.0*fl[0])*limTheta(rVal,1.0,cflL); 
  } else {
  rVal = (3.0*fr[15]-5.196152422706631*(fr[14]+fr[13])+5.196152422706631*fr[11]+9.0*fr[10]-9.0*(fr[7]+fr[6])+15.58845726811989*fr[3])/(36.0*EPSILON+1.732050807568877*fr[12]-3.0*(fr[9]+fr[8])+3.0*fr[5]+5.196152422706631*fr[4]-5.196152422706631*(fr[2]+fr[1])+9.0*fr[0]); 
  fqVal[4] = 0.02777777777777778*(1.732050807568877*fr[12]-3.0*(fr[9]+fr[8])+3.0*fr[5]+5.196152422706631*fr[4]-5.196152422706631*(fr[2]+fr[1])+9.0*fr[0])*limTheta(rVal,-1.0,cflR); 
  }
  if((-0.06804138174397717*alpha[7])-0.1178511301977579*alpha[6]+0.1178511301977579*alpha[5]-0.1178511301977579*alpha[4]+0.2041241452319315*alpha[3]-0.2041241452319315*alpha[2]+0.2041241452319315*alpha[1]+0.3535533905932737*alpha[0] > 0) {
  rVal = -(1.0*(3.0*fl[15]+5.196152422706631*fl[14]-5.196152422706631*fl[13]+5.196152422706631*fl[11]-9.0*fl[10]+9.0*fl[7]-9.0*fl[6]-15.58845726811989*fl[3]))/(36.0*EPSILON-1.732050807568877*fl[12]-3.0*fl[9]+3.0*fl[8]-3.0*fl[5]+5.196152422706631*fl[4]-5.196152422706631*fl[2]+5.196152422706631*fl[1]+9.0*fl[0]); 
  fqVal[5] = -0.02777777777777778*(1.732050807568877*fl[12]+3.0*fl[9]-3.0*fl[8]+3.0*fl[5]-5.196152422706631*fl[4]+5.196152422706631*fl[2]-5.196152422706631*fl[1]-9.0*fl[0])*limTheta(rVal,1.0,cflL); 
  } else {
  rVal = -(1.0*(3.0*fr[15]+5.196152422706631*fr[14]-5.196152422706631*fr[13]+5.196152422706631*fr[11]-9.0*fr[10]+9.0*fr[7]-9.0*fr[6]-15.58845726811989*fr[3]))/(36.0*EPSILON-1.732050807568877*fr[12]-3.0*fr[9]+3.0*fr[8]-3.0*fr[5]+5.196152422706631*fr[4]-5.196152422706631*fr[2]+5.196152422706631*fr[1]+9.0*fr[0]); 
  fqVal[5] = -0.02777777777777778*(1.732050807568877*fr[12]+3.0*fr[9]-3.0*fr[8]+3.0*fr[5]-5.196152422706631*fr[4]+5.196152422706631*fr[2]-5.196152422706631*fr[1]-9.0*fr[0])*limTheta(rVal,-1.0,cflR); 
  }
  if((-0.06804138174397717*alpha[7])+0.1178511301977579*alpha[6]-0.1178511301977579*(alpha[5]+alpha[4])+0.2041241452319315*(alpha[3]+alpha[2])-0.2041241452319315*alpha[1]+0.3535533905932737*alpha[0] > 0) {
  rVal = -(1.0*(3.0*fl[15]-5.196152422706631*fl[14]+5.196152422706631*(fl[13]+fl[11])-9.0*(fl[10]+fl[7])+9.0*fl[6]-15.58845726811989*fl[3]))/(36.0*EPSILON-1.732050807568877*fl[12]+3.0*fl[9]-3.0*(fl[8]+fl[5])+5.196152422706631*(fl[4]+fl[2])-5.196152422706631*fl[1]+9.0*fl[0]); 
  fqVal[6] = -0.02777777777777778*(1.732050807568877*fl[12]-3.0*fl[9]+3.0*(fl[8]+fl[5])-5.196152422706631*(fl[4]+fl[2])+5.196152422706631*fl[1]-9.0*fl[0])*limTheta(rVal,1.0,cflL); 
  } else {
  rVal = -(1.0*(3.0*fr[15]-5.196152422706631*fr[14]+5.196152422706631*(fr[13]+fr[11])-9.0*(fr[10]+fr[7])+9.0*fr[6]-15.58845726811989*fr[3]))/(36.0*EPSILON-1.732050807568877*fr[12]+3.0*fr[9]-3.0*(fr[8]+fr[5])+5.196152422706631*(fr[4]+fr[2])-5.196152422706631*fr[1]+9.0*fr[0]); 
  fqVal[6] = -0.02777777777777778*(1.732050807568877*fr[12]-3.0*fr[9]+3.0*(fr[8]+fr[5])-5.196152422706631*(fr[4]+fr[2])+5.196152422706631*fr[1]-9.0*fr[0])*limTheta(rVal,-1.0,cflR); 
  }
  if(0.06804138174397717*alpha[7]+0.1178511301977579*(alpha[6]+alpha[5]+alpha[4])+0.2041241452319315*(alpha[3]+alpha[2]+alpha[1])+0.3535533905932737*alpha[0] > 0) {
  rVal = (3.0*fl[15]+5.196152422706631*(fl[14]+fl[13]+fl[11])+9.0*(fl[10]+fl[7]+fl[6])+15.58845726811989*fl[3])/(36.0*EPSILON+1.732050807568877*fl[12]+3.0*(fl[9]+fl[8]+fl[5])+5.196152422706631*(fl[4]+fl[2]+fl[1])+9.0*fl[0]); 
  fqVal[7] = 0.02777777777777778*(1.732050807568877*fl[12]+3.0*(fl[9]+fl[8]+fl[5])+5.196152422706631*(fl[4]+fl[2]+fl[1])+9.0*fl[0])*limTheta(rVal,1.0,cflL); 
  } else {
  rVal = (3.0*fr[15]+5.196152422706631*(fr[14]+fr[13]+fr[11])+9.0*(fr[10]+fr[7]+fr[6])+15.58845726811989*fr[3])/(36.0*EPSILON+1.732050807568877*fr[12]+3.0*(fr[9]+fr[8]+fr[5])+5.196152422706631*(fr[4]+fr[2]+fr[1])+9.0*fr[0]); 
  fqVal[7] = 0.02777777777777778*(1.732050807568877*fr[12]+3.0*(fr[9]+fr[8]+fr[5])+5.196152422706631*(fr[4]+fr[2]+fr[1])+9.0*fr[0])*limTheta(rVal,-1.0,cflR); 
  }
  double fhatALVal[8];  // fhatALVal = mode coefficients of anti-limited f 
  fhatALVal[0] = 0.3535533905932737*(fqVal[7]+fqVal[6]+fqVal[5]+fqVal[4]+fqVal[3]+fqVal[2]+fqVal[1]+fqVal[0]); 
  fhatALVal[1] = 0.6123724356957944*(fqVal[7]-1.0*fqVal[6]+fqVal[5]-1.0*fqVal[4]+fqVal[3]-1.0*fqVal[2]+fqVal[1]-1.0*fqVal[0]); 
  fhatALVal[2] = 0.6123724356957944*(fqVal[7]+fqVal[6]-1.0*(fqVal[5]+fqVal[4])+fqVal[3]+fqVal[2]-1.0*(fqVal[1]+fqVal[0])); 
  fhatALVal[3] = 0.6123724356957944*(fqVal[7]+fqVal[6]+fqVal[5]+fqVal[4]-1.0*(fqVal[3]+fqVal[2]+fqVal[1]+fqVal[0])); 
  fhatALVal[4] = 1.060660171779821*(fqVal[7]-1.0*(fqVal[6]+fqVal[5])+fqVal[4]+fqVal[3]-1.0*(fqVal[2]+fqVal[1])+fqVal[0]); 
  fhatALVal[5] = 1.060660171779821*(fqVal[7]-1.0*fqVal[6]+fqVal[5]-1.0*(fqVal[4]+fqVal[3])+fqVal[2]-1.0*fqVal[1]+fqVal[0]); 
  fhatALVal[6] = 1.060660171779821*(fqVal[7]+fqVal[6]-1.0*(fqVal[5]+fqVal[4]+fqVal[3]+fqVal[2])+fqVal[1]+fqVal[0]); 
  fhatALVal[7] = 1.837117307087383*(fqVal[7]-1.0*(fqVal[6]+fqVal[5])+fqVal[4]-1.0*fqVal[3]+fqVal[2]+fqVal[1]-1.0*fqVal[0]); 
  incr[0] = 0.25*(alpha[7]*fhatALVal[7]+alpha[6]*fhatALVal[6]+alpha[5]*fhatALVal[5]+alpha[4]*fhatALVal[4]+alpha[3]*fhatALVal[3]+alpha[2]*fhatALVal[2]+alpha[1]*fhatALVal[1]+alpha[0]*fhatALVal[0])*dfac_v; 
  incr[1] = 0.25*(alpha[6]*fhatALVal[7]+fhatALVal[6]*alpha[7]+alpha[3]*fhatALVal[5]+fhatALVal[3]*alpha[5]+alpha[2]*fhatALVal[4]+fhatALVal[2]*alpha[4]+alpha[0]*fhatALVal[1]+fhatALVal[0]*alpha[1])*dfac_v; 
  incr[2] = 0.25*(alpha[5]*fhatALVal[7]+fhatALVal[5]*alpha[7]+alpha[3]*fhatALVal[6]+fhatALVal[3]*alpha[6]+alpha[1]*fhatALVal[4]+fhatALVal[1]*alpha[4]+alpha[0]*fhatALVal[2]+fhatALVal[0]*alpha[2])*dfac_v; 
  incr[3] = -0.4330127018922193*(alpha[7]*fhatALVal[7]+alpha[6]*fhatALVal[6]+alpha[5]*fhatALVal[5]+alpha[4]*fhatALVal[4]+alpha[3]*fhatALVal[3]+alpha[2]*fhatALVal[2]+alpha[1]*fhatALVal[1]+alpha[0]*fhatALVal[0])*dfac_v; 
  incr[4] = 0.25*(alpha[4]*fhatALVal[7]+fhatALVal[4]*alpha[7]+alpha[2]*fhatALVal[6]+fhatALVal[2]*alpha[6]+alpha[1]*fhatALVal[5]+fhatALVal[1]*alpha[5]+alpha[0]*fhatALVal[3]+fhatALVal[0]*alpha[3])*dfac_v; 
  incr[5] = 0.25*(alpha[3]*fhatALVal[7]+fhatALVal[3]*alpha[7]+alpha[5]*fhatALVal[6]+fhatALVal[5]*alpha[6]+alpha[0]*fhatALVal[4]+fhatALVal[0]*alpha[4]+alpha[1]*fhatALVal[2]+fhatALVal[1]*alpha[2])*dfac_v; 
  incr[6] = -0.4330127018922193*(alpha[6]*fhatALVal[7]+fhatALVal[6]*alpha[7]+alpha[3]*fhatALVal[5]+fhatALVal[3]*alpha[5]+alpha[2]*fhatALVal[4]+fhatALVal[2]*alpha[4]+alpha[0]*fhatALVal[1]+fhatALVal[0]*alpha[1])*dfac_v; 
  incr[7] = -0.4330127018922193*(alpha[5]*fhatALVal[7]+fhatALVal[5]*alpha[7]+alpha[3]*fhatALVal[6]+fhatALVal[3]*alpha[6]+alpha[1]*fhatALVal[4]+fhatALVal[1]*alpha[4]+alpha[0]*fhatALVal[2]+fhatALVal[0]*alpha[2])*dfac_v; 
  incr[8] = 0.25*(alpha[2]*fhatALVal[7]+fhatALVal[2]*alpha[7]+alpha[4]*fhatALVal[6]+fhatALVal[4]*alpha[6]+alpha[0]*fhatALVal[5]+fhatALVal[0]*alpha[5]+alpha[1]*fhatALVal[3]+fhatALVal[1]*alpha[3])*dfac_v; 
  incr[9] = 0.25*(alpha[1]*fhatALVal[7]+fhatALVal[1]*alpha[7]+alpha[0]*fhatALVal[6]+fhatALVal[0]*alpha[6]+alpha[4]*fhatALVal[5]+fhatALVal[4]*alpha[5]+alpha[2]*fhatALVal[3]+fhatALVal[2]*alpha[3])*dfac_v; 
  incr[10] = -0.4330127018922193*(alpha[4]*fhatALVal[7]+fhatALVal[4]*alpha[7]+alpha[2]*fhatALVal[6]+fhatALVal[2]*alpha[6]+alpha[1]*fhatALVal[5]+fhatALVal[1]*alpha[5]+alpha[0]*fhatALVal[3]+fhatALVal[0]*alpha[3])*dfac_v; 
  incr[11] = -0.4330127018922193*(alpha[3]*fhatALVal[7]+fhatALVal[3]*alpha[7]+alpha[5]*fhatALVal[6]+fhatALVal[5]*alpha[6]+alpha[0]*fhatALVal[4]+fhatALVal[0]*alpha[4]+alpha[1]*fhatALVal[2]+fhatALVal[1]*alpha[2])*dfac_v; 
  incr[12] = 0.25*(alpha[0]*fhatALVal[7]+fhatALVal[0]*alpha[7]+alpha[1]*fhatALVal[6]+fhatALVal[1]*alpha[6]+alpha[2]*fhatALVal[5]+fhatALVal[2]*alpha[5]+alpha[3]*fhatALVal[4]+fhatALVal[3]*alpha[4])*dfac_v; 
  incr[13] = -0.4330127018922193*(alpha[2]*fhatALVal[7]+fhatALVal[2]*alpha[7]+alpha[4]*fhatALVal[6]+fhatALVal[4]*alpha[6]+alpha[0]*fhatALVal[5]+fhatALVal[0]*alpha[5]+alpha[1]*fhatALVal[3]+fhatALVal[1]*alpha[3])*dfac_v; 
  incr[14] = -0.4330127018922193*(alpha[1]*fhatALVal[7]+fhatALVal[1]*alpha[7]+alpha[0]*fhatALVal[6]+fhatALVal[0]*alpha[6]+alpha[4]*fhatALVal[5]+fhatALVal[4]*alpha[5]+alpha[2]*fhatALVal[3]+fhatALVal[2]*alpha[3])*dfac_v; 
  incr[15] = -0.4330127018922193*(alpha[0]*fhatALVal[7]+fhatALVal[0]*alpha[7]+alpha[1]*fhatALVal[6]+fhatALVal[1]*alpha[6]+alpha[2]*fhatALVal[5]+fhatALVal[2]*alpha[5]+alpha[3]*fhatALVal[4]+fhatALVal[3]*alpha[4])*dfac_v; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 
  outr[4] += incr[4]; 
  outr[5] += incr[5]; 
  outr[6] += incr[6]; 
  outr[7] += incr[7]; 
  outr[8] += incr[8]; 
  outr[9] += incr[9]; 
  outr[10] += incr[10]; 
  outr[11] += incr[11]; 
  outr[12] += incr[12]; 
  outr[13] += incr[13]; 
  outr[14] += incr[14]; 
  outr[15] += incr[15]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += -1.0*incr[1]; 
  outl[2] += -1.0*incr[2]; 
  outl[3] += incr[3]; 
  outl[4] += -1.0*incr[4]; 
  outl[5] += -1.0*incr[5]; 
  outl[6] += incr[6]; 
  outl[7] += incr[7]; 
  outl[8] += -1.0*incr[8]; 
  outl[9] += -1.0*incr[9]; 
  outl[10] += incr[10]; 
  outl[11] += incr[11]; 
  outl[12] += -1.0*incr[12]; 
  outl[13] += incr[13]; 
  outl[14] += incr[14]; 
  outl[15] += incr[15]; 
  return std::abs(alpha0); 
} 
