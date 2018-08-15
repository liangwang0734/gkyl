#include <GyrokineticModDecl.h> 
double GyrokineticSurfPositivity1x1vSer_X_P1_Bvars_0(const double q_, const double m_, const double dt, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[4]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = wv; 

  double alpha[4]; 
  alpha[0] = 2.0*wv; 
  if (alpha0>0) { 
  double rVal[4];  // rVal=f1/f0 at each node 
  rVal[0] = -(1.0*(1.732050807568877*fl[3]-3.0*fl[1]))/(3.464101615137754*EPSILON-1.0*fl[2]+1.732050807568877*fl[0]); 
  rVal[1] = -(1.0*(1.732050807568877*fl[3]-3.0*fl[1]))/(3.464101615137754*EPSILON-1.0*fl[2]+1.732050807568877*fl[0]); 
  rVal[2] = (1.732050807568877*fl[3]+3.0*fl[1])/(3.464101615137754*EPSILON+fl[2]+1.732050807568877*fl[0]); 
  rVal[3] = (1.732050807568877*fl[3]+3.0*fl[1])/(3.464101615137754*EPSILON+fl[2]+1.732050807568877*fl[0]); 
  double fqVal[4];  // fqVal = anti-limited f evaluated at each node 
  double cflVal = std::abs((alpha0*dt)/dxv[0]); 
  fqVal[0] = -0.2886751345948129*(fl[2]-1.732050807568877*fl[0])*limTheta(rVal[0],1.0,cflVal); 
  fqVal[1] = -0.2886751345948129*(fl[2]-1.732050807568877*fl[0])*limTheta(rVal[1],1.0,cflVal); 
  fqVal[2] = 0.2886751345948129*(fl[2]+1.732050807568877*fl[0])*limTheta(rVal[2],1.0,cflVal); 
  fqVal[3] = 0.2886751345948129*(fl[2]+1.732050807568877*fl[0])*limTheta(rVal[3],1.0,cflVal); 
  double fhatALVal[4];  // fhatALVal = mode coefficients of anti-limited f 
  fhatALVal[0] = 0.5*(fqVal[3]+fqVal[2]+fqVal[1]+fqVal[0]); 
  fhatALVal[1] = 0.8660254037844386*(fqVal[3]-1.0*fqVal[2]+fqVal[1]-1.0*fqVal[0]); 
  fhatALVal[2] = 0.8660254037844386*(fqVal[3]+fqVal[2]-1.0*(fqVal[1]+fqVal[0])); 
  fhatALVal[3] = 1.5*(fqVal[3]-1.0*(fqVal[2]+fqVal[1])+fqVal[0]); 
  incr[0] = 0.25*alpha[0]*(1.732050807568877*fhatALVal[1]+fhatALVal[0])*dfac_x; 
  incr[1] = -0.25*alpha[0]*(3.0*fhatALVal[1]+1.732050807568877*fhatALVal[0])*dfac_x; 
  incr[2] = 0.25*alpha[0]*(1.732050807568877*fhatALVal[3]+fhatALVal[2])*dfac_x; 
  incr[3] = -0.25*alpha[0]*(3.0*fhatALVal[3]+1.732050807568877*fhatALVal[2])*dfac_x; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += incr[1]; 
  outl[2] += -1.0*incr[2]; 
  outl[3] += incr[3]; 
  } else { 
  double rVal[4];  // rVal=f1/f0 at each node 
  rVal[0] = -(1.0*(1.732050807568877*fr[3]-3.0*fr[1]))/(3.464101615137754*EPSILON-1.0*fr[2]+1.732050807568877*fr[0]); 
  rVal[1] = -(1.0*(1.732050807568877*fr[3]-3.0*fr[1]))/(3.464101615137754*EPSILON-1.0*fr[2]+1.732050807568877*fr[0]); 
  rVal[2] = (1.732050807568877*fr[3]+3.0*fr[1])/(3.464101615137754*EPSILON+fr[2]+1.732050807568877*fr[0]); 
  rVal[3] = (1.732050807568877*fr[3]+3.0*fr[1])/(3.464101615137754*EPSILON+fr[2]+1.732050807568877*fr[0]); 
  double fqVal[4];  // fqVal = anti-limited f evaluated at each node 
  double cflVal = std::abs((alpha0*dt)/dxv[0]); 
  fqVal[0] = -0.2886751345948129*(fr[2]-1.732050807568877*fr[0])*limTheta(rVal[0],-1.0,cflVal); 
  fqVal[1] = -0.2886751345948129*(fr[2]-1.732050807568877*fr[0])*limTheta(rVal[1],-1.0,cflVal); 
  fqVal[2] = 0.2886751345948129*(fr[2]+1.732050807568877*fr[0])*limTheta(rVal[2],-1.0,cflVal); 
  fqVal[3] = 0.2886751345948129*(fr[2]+1.732050807568877*fr[0])*limTheta(rVal[3],-1.0,cflVal); 
  double fhatALVal[4];  // fhatALVal = mode coefficients of anti-limited f 
  fhatALVal[0] = 0.5*(fqVal[3]+fqVal[2]+fqVal[1]+fqVal[0]); 
  fhatALVal[1] = 0.8660254037844386*(fqVal[3]-1.0*fqVal[2]+fqVal[1]-1.0*fqVal[0]); 
  fhatALVal[2] = 0.8660254037844386*(fqVal[3]+fqVal[2]-1.0*(fqVal[1]+fqVal[0])); 
  fhatALVal[3] = 1.5*(fqVal[3]-1.0*(fqVal[2]+fqVal[1])+fqVal[0]); 
  incr[0] = -0.25*alpha[0]*(1.732050807568877*fhatALVal[1]-1.0*fhatALVal[0])*dfac_x; 
  incr[1] = 0.25*alpha[0]*(3.0*fhatALVal[1]-1.732050807568877*fhatALVal[0])*dfac_x; 
  incr[2] = -0.25*alpha[0]*(1.732050807568877*fhatALVal[3]-1.0*fhatALVal[2])*dfac_x; 
  incr[3] = 0.25*alpha[0]*(3.0*fhatALVal[3]-1.732050807568877*fhatALVal[2])*dfac_x; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += incr[1]; 
  outl[2] += -1.0*incr[2]; 
  outl[3] += incr[3]; 
  } 
  return std::abs(alpha0); 
} 
double GyrokineticSurfPositivity1x1vSer_X_P2_Bvars_0(const double q_, const double m_, const double dt, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[8]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = wv; 

  double amax = 0.0; 
  bool upwind = false; 
  if(upwind) 
    amax = std::abs(alpha0); 
  else 
    amax = amax_in; 

  double Ghat[8]; 
  Ghat[0] = 1.118033988749895*fr[4]*wv+1.118033988749895*fl[4]*wv-0.8660254037844386*fr[1]*wv+0.8660254037844386*fl[1]*wv+0.5*fr[0]*wv+0.5*fl[0]*wv+(0.6454972243679028*fr[6])/dfac_v+(0.6454972243679028*fl[6])/dfac_v-(0.5*fr[3])/dfac_v+(0.5*fl[3])/dfac_v+(0.2886751345948129*fr[2])/dfac_v+(0.2886751345948129*fl[2])/dfac_v-1.118033988749895*fr[4]*amax+1.118033988749895*fl[4]*amax+0.8660254037844386*fr[1]*amax+0.8660254037844386*fl[1]*amax-0.5*fr[0]*amax+0.5*fl[0]*amax; 
  Ghat[2] = 1.118033988749895*fr[6]*wv+1.118033988749895*fl[6]*wv-0.8660254037844386*fr[3]*wv+0.8660254037844386*fl[3]*wv+0.5*fr[2]*wv+0.5*fl[2]*wv-(0.447213595499958*fr[7])/dfac_v+(0.447213595499958*fl[7])/dfac_v+(0.2581988897471612*fr[5])/dfac_v+(0.2581988897471612*fl[5])/dfac_v+(0.6454972243679029*fr[4])/dfac_v+(0.6454972243679029*fl[4])/dfac_v-(0.5*fr[1])/dfac_v+(0.5*fl[1])/dfac_v+(0.2886751345948129*fr[0])/dfac_v+(0.2886751345948129*fl[0])/dfac_v-1.118033988749895*fr[6]*amax+1.118033988749895*fl[6]*amax+0.8660254037844386*fr[3]*amax+0.8660254037844386*fl[3]*amax-0.5*fr[2]*amax+0.5*fl[2]*amax; 
  Ghat[5] = (-0.8660254037844387*fr[7]*wv)+0.8660254037844387*fl[7]*wv+0.5*fr[5]*wv+0.5*fl[5]*wv+(0.5773502691896257*fr[6])/dfac_v+(0.5773502691896257*fl[6])/dfac_v-(0.4472135954999579*fr[3])/dfac_v+(0.4472135954999579*fl[3])/dfac_v+(0.2581988897471612*fr[2])/dfac_v+(0.2581988897471612*fl[2])/dfac_v+0.8660254037844387*fr[7]*amax+0.8660254037844387*fl[7]*amax-0.5*fr[5]*amax+0.5*fl[5]*amax; 

  incr[0] = 0.5*Ghat[0]*dfac_x; 
  incr[1] = -0.8660254037844386*Ghat[0]*dfac_x; 
  incr[2] = 0.5*Ghat[2]*dfac_x; 
  incr[3] = -0.8660254037844386*Ghat[2]*dfac_x; 
  incr[4] = 1.118033988749895*Ghat[0]*dfac_x; 
  incr[5] = 0.5*Ghat[5]*dfac_x; 
  incr[6] = 1.118033988749895*Ghat[2]*dfac_x; 
  incr[7] = -0.8660254037844387*Ghat[5]*dfac_x; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 
  outr[4] += incr[4]; 
  outr[5] += incr[5]; 
  outr[6] += incr[6]; 
  outr[7] += incr[7]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += incr[1]; 
  outl[2] += -1.0*incr[2]; 
  outl[3] += incr[3]; 
  outl[4] += -1.0*incr[4]; 
  outl[5] += -1.0*incr[5]; 
  outl[6] += -1.0*incr[6]; 
  outl[7] += incr[7]; 
return std::abs(alpha0); 
} 
double GyrokineticSurfPositivity1x1vSer_Vpar_P1_Bvars_0(const double q_, const double m_, const double dt, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[4]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = -(1.224744871391589*Phi[1]*dfac_x*q_)/m_; 

  double alpha[4]; 
  alpha[0] = -(2.449489742783178*Phi[1]*dfac_x*q_)/m_; 
  if (alpha0>0) { 
  double rVal[4];  // rVal=f1/f0 at each node 
  rVal[0] = -(1.0*(1.732050807568877*fl[3]-3.0*fl[2]))/(3.464101615137754*EPSILON-1.0*fl[1]+1.732050807568877*fl[0]); 
  rVal[1] = (1.732050807568877*fl[3]+3.0*fl[2])/(3.464101615137754*EPSILON+fl[1]+1.732050807568877*fl[0]); 
  rVal[2] = -(1.0*(1.732050807568877*fl[3]-3.0*fl[2]))/(3.464101615137754*EPSILON-1.0*fl[1]+1.732050807568877*fl[0]); 
  rVal[3] = (1.732050807568877*fl[3]+3.0*fl[2])/(3.464101615137754*EPSILON+fl[1]+1.732050807568877*fl[0]); 
  double fqVal[4];  // fqVal = anti-limited f evaluated at each node 
  double cflVal = std::abs((alpha0*dt)/dxv[1]); 
  fqVal[0] = -0.2886751345948129*(fl[1]-1.732050807568877*fl[0])*limTheta(rVal[0],1.0,cflVal); 
  fqVal[1] = 0.2886751345948129*(fl[1]+1.732050807568877*fl[0])*limTheta(rVal[1],1.0,cflVal); 
  fqVal[2] = -0.2886751345948129*(fl[1]-1.732050807568877*fl[0])*limTheta(rVal[2],1.0,cflVal); 
  fqVal[3] = 0.2886751345948129*(fl[1]+1.732050807568877*fl[0])*limTheta(rVal[3],1.0,cflVal); 
  double fhatALVal[4];  // fhatALVal = mode coefficients of anti-limited f 
  fhatALVal[0] = 0.5*(fqVal[3]+fqVal[2]+fqVal[1]+fqVal[0]); 
  fhatALVal[1] = 0.8660254037844386*(fqVal[3]-1.0*fqVal[2]+fqVal[1]-1.0*fqVal[0]); 
  fhatALVal[2] = 0.8660254037844386*(fqVal[3]+fqVal[2]-1.0*(fqVal[1]+fqVal[0])); 
  fhatALVal[3] = 1.5*(fqVal[3]-1.0*(fqVal[2]+fqVal[1])+fqVal[0]); 
  incr[0] = 0.25*alpha[0]*(1.732050807568877*fhatALVal[2]+fhatALVal[0])*dfac_v; 
  incr[1] = 0.25*alpha[0]*(1.732050807568877*fhatALVal[3]+fhatALVal[1])*dfac_v; 
  incr[2] = -0.25*alpha[0]*(3.0*fhatALVal[2]+1.732050807568877*fhatALVal[0])*dfac_v; 
  incr[3] = -0.25*alpha[0]*(3.0*fhatALVal[3]+1.732050807568877*fhatALVal[1])*dfac_v; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += -1.0*incr[1]; 
  outl[2] += incr[2]; 
  outl[3] += incr[3]; 
  } else { 
  double rVal[4];  // rVal=f1/f0 at each node 
  rVal[0] = -(1.0*(1.732050807568877*fr[3]-3.0*fr[2]))/(3.464101615137754*EPSILON-1.0*fr[1]+1.732050807568877*fr[0]); 
  rVal[1] = (1.732050807568877*fr[3]+3.0*fr[2])/(3.464101615137754*EPSILON+fr[1]+1.732050807568877*fr[0]); 
  rVal[2] = -(1.0*(1.732050807568877*fr[3]-3.0*fr[2]))/(3.464101615137754*EPSILON-1.0*fr[1]+1.732050807568877*fr[0]); 
  rVal[3] = (1.732050807568877*fr[3]+3.0*fr[2])/(3.464101615137754*EPSILON+fr[1]+1.732050807568877*fr[0]); 
  double fqVal[4];  // fqVal = anti-limited f evaluated at each node 
  double cflVal = std::abs((alpha0*dt)/dxv[1]); 
  fqVal[0] = -0.2886751345948129*(fr[1]-1.732050807568877*fr[0])*limTheta(rVal[0],-1.0,cflVal); 
  fqVal[1] = 0.2886751345948129*(fr[1]+1.732050807568877*fr[0])*limTheta(rVal[1],-1.0,cflVal); 
  fqVal[2] = -0.2886751345948129*(fr[1]-1.732050807568877*fr[0])*limTheta(rVal[2],-1.0,cflVal); 
  fqVal[3] = 0.2886751345948129*(fr[1]+1.732050807568877*fr[0])*limTheta(rVal[3],-1.0,cflVal); 
  double fhatALVal[4];  // fhatALVal = mode coefficients of anti-limited f 
  fhatALVal[0] = 0.5*(fqVal[3]+fqVal[2]+fqVal[1]+fqVal[0]); 
  fhatALVal[1] = 0.8660254037844386*(fqVal[3]-1.0*fqVal[2]+fqVal[1]-1.0*fqVal[0]); 
  fhatALVal[2] = 0.8660254037844386*(fqVal[3]+fqVal[2]-1.0*(fqVal[1]+fqVal[0])); 
  fhatALVal[3] = 1.5*(fqVal[3]-1.0*(fqVal[2]+fqVal[1])+fqVal[0]); 
  incr[0] = -0.25*alpha[0]*(1.732050807568877*fhatALVal[2]-1.0*fhatALVal[0])*dfac_v; 
  incr[1] = -0.25*alpha[0]*(1.732050807568877*fhatALVal[3]-1.0*fhatALVal[1])*dfac_v; 
  incr[2] = 0.25*alpha[0]*(3.0*fhatALVal[2]-1.732050807568877*fhatALVal[0])*dfac_v; 
  incr[3] = 0.25*alpha[0]*(3.0*fhatALVal[3]-1.732050807568877*fhatALVal[1])*dfac_v; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += -1.0*incr[1]; 
  outl[2] += incr[2]; 
  outl[3] += incr[3]; 
  } 
  return std::abs(alpha0); 
} 
double GyrokineticSurfPositivity1x1vSer_Vpar_P2_Bvars_0(const double q_, const double m_, const double dt, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[8]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = -(1.224744871391589*Phi[1]*dfac_x*q_)/m_; 

  double amax = 0.0; 
  bool upwind = false; 
  if(upwind) 
    amax = std::abs(alpha0); 
  else 
    amax = amax_in; 

  double Ghat[8]; 
  Ghat[0] = (-(3.061862178478972*Phi[2]*fr[7]*dfac_x*q_)/m_)-(3.061862178478972*Phi[2]*fl[7]*dfac_x*q_)/m_-(1.369306393762915*Phi[1]*fr[5]*dfac_x*q_)/m_-(1.369306393762915*Phi[1]*fl[5]*dfac_x*q_)/m_+(2.371708245126284*Phi[2]*fr[3]*dfac_x*q_)/m_-(2.371708245126284*Phi[2]*fl[3]*dfac_x*q_)/m_+(1.060660171779821*Phi[1]*fr[2]*dfac_x*q_)/m_-(1.060660171779821*Phi[1]*fl[2]*dfac_x*q_)/m_-(1.369306393762915*fr[1]*Phi[2]*dfac_x*q_)/m_-(1.369306393762915*fl[1]*Phi[2]*dfac_x*q_)/m_-(0.6123724356957944*fr[0]*Phi[1]*dfac_x*q_)/m_-(0.6123724356957944*fl[0]*Phi[1]*dfac_x*q_)/m_-1.118033988749895*fr[5]*amax+1.118033988749895*fl[5]*amax+0.8660254037844386*fr[2]*amax+0.8660254037844386*fl[2]*amax-0.5*fr[0]*amax+0.5*fl[0]*amax; 
  Ghat[1] = (-(1.369306393762915*Phi[1]*fr[7]*dfac_x*q_)/m_)-(1.369306393762915*Phi[1]*fl[7]*dfac_x*q_)/m_+(2.121320343559642*Phi[2]*fr[6]*dfac_x*q_)/m_-(2.121320343559642*Phi[2]*fl[6]*dfac_x*q_)/m_-(3.061862178478972*Phi[2]*fr[5]*dfac_x*q_)/m_-(3.061862178478972*Phi[2]*fl[5]*dfac_x*q_)/m_-(1.224744871391589*Phi[2]*fr[4]*dfac_x*q_)/m_-(1.224744871391589*Phi[2]*fl[4]*dfac_x*q_)/m_+(1.060660171779821*Phi[1]*fr[3]*dfac_x*q_)/m_-(1.060660171779821*Phi[1]*fl[3]*dfac_x*q_)/m_+(2.371708245126284*Phi[2]*fr[2]*dfac_x*q_)/m_-(2.371708245126284*Phi[2]*fl[2]*dfac_x*q_)/m_-(1.369306393762915*fr[0]*Phi[2]*dfac_x*q_)/m_-(1.369306393762915*fl[0]*Phi[2]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fr[1]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fl[1]*dfac_x*q_)/m_-1.118033988749895*fr[7]*amax+1.118033988749895*fl[7]*amax+0.8660254037844386*fr[3]*amax+0.8660254037844386*fl[3]*amax-0.5*fr[1]*amax+0.5*fl[1]*amax; 
  Ghat[4] = (-(2.738612787525831*Phi[2]*fr[7]*dfac_x*q_)/m_)-(2.738612787525831*Phi[2]*fl[7]*dfac_x*q_)/m_+(1.060660171779821*Phi[1]*fr[6]*dfac_x*q_)/m_-(1.060660171779821*Phi[1]*fl[6]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fr[4]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fl[4]*dfac_x*q_)/m_+(2.121320343559642*Phi[2]*fr[3]*dfac_x*q_)/m_-(2.121320343559642*Phi[2]*fl[3]*dfac_x*q_)/m_-(1.224744871391589*fr[1]*Phi[2]*dfac_x*q_)/m_-(1.224744871391589*fl[1]*Phi[2]*dfac_x*q_)/m_+0.8660254037844387*fr[6]*amax+0.8660254037844387*fl[6]*amax-0.5*fr[4]*amax+0.5*fl[4]*amax; 

  incr[0] = 0.5*Ghat[0]*dfac_v; 
  incr[1] = 0.5*Ghat[1]*dfac_v; 
  incr[2] = -0.8660254037844386*Ghat[0]*dfac_v; 
  incr[3] = -0.8660254037844386*Ghat[1]*dfac_v; 
  incr[4] = 0.5*Ghat[4]*dfac_v; 
  incr[5] = 1.118033988749895*Ghat[0]*dfac_v; 
  incr[6] = -0.8660254037844387*Ghat[4]*dfac_v; 
  incr[7] = 1.118033988749895*Ghat[1]*dfac_v; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 
  outr[4] += incr[4]; 
  outr[5] += incr[5]; 
  outr[6] += incr[6]; 
  outr[7] += incr[7]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += -1.0*incr[1]; 
  outl[2] += incr[2]; 
  outl[3] += incr[3]; 
  outl[4] += -1.0*incr[4]; 
  outl[5] += -1.0*incr[5]; 
  outl[6] += incr[6]; 
  outl[7] += -1.0*incr[7]; 
return std::abs(alpha0); 
} 
double EmGyrokineticSurfPositivity1x1vSer_X_P1_Bvars_0(const double q_, const double m_, const double dt, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *Apar, const double *dApardt, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[4]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = wv; 

  double alpha[4]; 
  alpha[0] = 2.0*wv; 
  if (alpha0>0) { 
  double rVal[4];  // rVal=f1/f0 at each node 
  rVal[0] = -(1.0*(1.732050807568877*fl[3]-3.0*fl[1]))/(3.464101615137754*EPSILON-1.0*fl[2]+1.732050807568877*fl[0]); 
  rVal[1] = -(1.0*(1.732050807568877*fl[3]-3.0*fl[1]))/(3.464101615137754*EPSILON-1.0*fl[2]+1.732050807568877*fl[0]); 
  rVal[2] = (1.732050807568877*fl[3]+3.0*fl[1])/(3.464101615137754*EPSILON+fl[2]+1.732050807568877*fl[0]); 
  rVal[3] = (1.732050807568877*fl[3]+3.0*fl[1])/(3.464101615137754*EPSILON+fl[2]+1.732050807568877*fl[0]); 
  double fqVal[4];  // fqVal = anti-limited f evaluated at each node 
  double cflVal = std::abs((alpha0*dt)/dxv[0]); 
  fqVal[0] = -0.2886751345948129*(fl[2]-1.732050807568877*fl[0])*limTheta(rVal[0],1.0,cflVal); 
  fqVal[1] = -0.2886751345948129*(fl[2]-1.732050807568877*fl[0])*limTheta(rVal[1],1.0,cflVal); 
  fqVal[2] = 0.2886751345948129*(fl[2]+1.732050807568877*fl[0])*limTheta(rVal[2],1.0,cflVal); 
  fqVal[3] = 0.2886751345948129*(fl[2]+1.732050807568877*fl[0])*limTheta(rVal[3],1.0,cflVal); 
  double fhatALVal[4];  // fhatALVal = mode coefficients of anti-limited f 
  fhatALVal[0] = 0.5*(fqVal[3]+fqVal[2]+fqVal[1]+fqVal[0]); 
  fhatALVal[1] = 0.8660254037844386*(fqVal[3]-1.0*fqVal[2]+fqVal[1]-1.0*fqVal[0]); 
  fhatALVal[2] = 0.8660254037844386*(fqVal[3]+fqVal[2]-1.0*(fqVal[1]+fqVal[0])); 
  fhatALVal[3] = 1.5*(fqVal[3]-1.0*(fqVal[2]+fqVal[1])+fqVal[0]); 
  incr[0] = 0.25*alpha[0]*(1.732050807568877*fhatALVal[1]+fhatALVal[0])*dfac_x; 
  incr[1] = -0.25*alpha[0]*(3.0*fhatALVal[1]+1.732050807568877*fhatALVal[0])*dfac_x; 
  incr[2] = 0.25*alpha[0]*(1.732050807568877*fhatALVal[3]+fhatALVal[2])*dfac_x; 
  incr[3] = -0.25*alpha[0]*(3.0*fhatALVal[3]+1.732050807568877*fhatALVal[2])*dfac_x; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += incr[1]; 
  outl[2] += -1.0*incr[2]; 
  outl[3] += incr[3]; 
  } else { 
  double rVal[4];  // rVal=f1/f0 at each node 
  rVal[0] = -(1.0*(1.732050807568877*fr[3]-3.0*fr[1]))/(3.464101615137754*EPSILON-1.0*fr[2]+1.732050807568877*fr[0]); 
  rVal[1] = -(1.0*(1.732050807568877*fr[3]-3.0*fr[1]))/(3.464101615137754*EPSILON-1.0*fr[2]+1.732050807568877*fr[0]); 
  rVal[2] = (1.732050807568877*fr[3]+3.0*fr[1])/(3.464101615137754*EPSILON+fr[2]+1.732050807568877*fr[0]); 
  rVal[3] = (1.732050807568877*fr[3]+3.0*fr[1])/(3.464101615137754*EPSILON+fr[2]+1.732050807568877*fr[0]); 
  double fqVal[4];  // fqVal = anti-limited f evaluated at each node 
  double cflVal = std::abs((alpha0*dt)/dxv[0]); 
  fqVal[0] = -0.2886751345948129*(fr[2]-1.732050807568877*fr[0])*limTheta(rVal[0],-1.0,cflVal); 
  fqVal[1] = -0.2886751345948129*(fr[2]-1.732050807568877*fr[0])*limTheta(rVal[1],-1.0,cflVal); 
  fqVal[2] = 0.2886751345948129*(fr[2]+1.732050807568877*fr[0])*limTheta(rVal[2],-1.0,cflVal); 
  fqVal[3] = 0.2886751345948129*(fr[2]+1.732050807568877*fr[0])*limTheta(rVal[3],-1.0,cflVal); 
  double fhatALVal[4];  // fhatALVal = mode coefficients of anti-limited f 
  fhatALVal[0] = 0.5*(fqVal[3]+fqVal[2]+fqVal[1]+fqVal[0]); 
  fhatALVal[1] = 0.8660254037844386*(fqVal[3]-1.0*fqVal[2]+fqVal[1]-1.0*fqVal[0]); 
  fhatALVal[2] = 0.8660254037844386*(fqVal[3]+fqVal[2]-1.0*(fqVal[1]+fqVal[0])); 
  fhatALVal[3] = 1.5*(fqVal[3]-1.0*(fqVal[2]+fqVal[1])+fqVal[0]); 
  incr[0] = -0.25*alpha[0]*(1.732050807568877*fhatALVal[1]-1.0*fhatALVal[0])*dfac_x; 
  incr[1] = 0.25*alpha[0]*(3.0*fhatALVal[1]-1.732050807568877*fhatALVal[0])*dfac_x; 
  incr[2] = -0.25*alpha[0]*(1.732050807568877*fhatALVal[3]-1.0*fhatALVal[2])*dfac_x; 
  incr[3] = 0.25*alpha[0]*(3.0*fhatALVal[3]-1.732050807568877*fhatALVal[2])*dfac_x; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += incr[1]; 
  outl[2] += -1.0*incr[2]; 
  outl[3] += incr[3]; 
  } 
  return std::abs(alpha0); 
} 
double EmGyrokineticSurfPositivity1x1vSer_X_P2_Bvars_0(const double q_, const double m_, const double dt, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *Apar, const double *dApardt, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[8]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = wv; 

  double amax = 0.0; 
  bool upwind = false; 
  if(upwind) 
    amax = std::abs(alpha0); 
  else 
    amax = amax_in; 

  double Ghat[8]; 
  Ghat[0] = 1.118033988749895*fr[4]*wv+1.118033988749895*fl[4]*wv-0.8660254037844386*fr[1]*wv+0.8660254037844386*fl[1]*wv+0.5*fr[0]*wv+0.5*fl[0]*wv+(0.6454972243679028*fr[6])/dfac_v+(0.6454972243679028*fl[6])/dfac_v-(0.5*fr[3])/dfac_v+(0.5*fl[3])/dfac_v+(0.2886751345948129*fr[2])/dfac_v+(0.2886751345948129*fl[2])/dfac_v-1.118033988749895*fr[4]*amax+1.118033988749895*fl[4]*amax+0.8660254037844386*fr[1]*amax+0.8660254037844386*fl[1]*amax-0.5*fr[0]*amax+0.5*fl[0]*amax; 
  Ghat[2] = 1.118033988749895*fr[6]*wv+1.118033988749895*fl[6]*wv-0.8660254037844386*fr[3]*wv+0.8660254037844386*fl[3]*wv+0.5*fr[2]*wv+0.5*fl[2]*wv-(0.447213595499958*fr[7])/dfac_v+(0.447213595499958*fl[7])/dfac_v+(0.2581988897471612*fr[5])/dfac_v+(0.2581988897471612*fl[5])/dfac_v+(0.6454972243679029*fr[4])/dfac_v+(0.6454972243679029*fl[4])/dfac_v-(0.5*fr[1])/dfac_v+(0.5*fl[1])/dfac_v+(0.2886751345948129*fr[0])/dfac_v+(0.2886751345948129*fl[0])/dfac_v-1.118033988749895*fr[6]*amax+1.118033988749895*fl[6]*amax+0.8660254037844386*fr[3]*amax+0.8660254037844386*fl[3]*amax-0.5*fr[2]*amax+0.5*fl[2]*amax; 
  Ghat[5] = (-0.8660254037844387*fr[7]*wv)+0.8660254037844387*fl[7]*wv+0.5*fr[5]*wv+0.5*fl[5]*wv+(0.5773502691896257*fr[6])/dfac_v+(0.5773502691896257*fl[6])/dfac_v-(0.4472135954999579*fr[3])/dfac_v+(0.4472135954999579*fl[3])/dfac_v+(0.2581988897471612*fr[2])/dfac_v+(0.2581988897471612*fl[2])/dfac_v+0.8660254037844387*fr[7]*amax+0.8660254037844387*fl[7]*amax-0.5*fr[5]*amax+0.5*fl[5]*amax; 

  incr[0] = 0.5*Ghat[0]*dfac_x; 
  incr[1] = -0.8660254037844386*Ghat[0]*dfac_x; 
  incr[2] = 0.5*Ghat[2]*dfac_x; 
  incr[3] = -0.8660254037844386*Ghat[2]*dfac_x; 
  incr[4] = 1.118033988749895*Ghat[0]*dfac_x; 
  incr[5] = 0.5*Ghat[5]*dfac_x; 
  incr[6] = 1.118033988749895*Ghat[2]*dfac_x; 
  incr[7] = -0.8660254037844387*Ghat[5]*dfac_x; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 
  outr[4] += incr[4]; 
  outr[5] += incr[5]; 
  outr[6] += incr[6]; 
  outr[7] += incr[7]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += incr[1]; 
  outl[2] += -1.0*incr[2]; 
  outl[3] += incr[3]; 
  outl[4] += -1.0*incr[4]; 
  outl[5] += -1.0*incr[5]; 
  outl[6] += -1.0*incr[6]; 
  outl[7] += incr[7]; 
return std::abs(alpha0); 
} 
double EmGyrokineticSurfPositivity1x1vSer_Vpar_P1_Bvars_0(const double q_, const double m_, const double dt, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *Apar, const double *dApardt, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[4]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = -(0.5*(2.449489742783178*Phi[1]*dfac_x+1.414213562373095*dApardt[0])*q_)/m_; 

  double alpha[4]; 
  alpha[0] = (-(2.449489742783178*Phi[1]*dfac_x*q_)/m_)-(1.414213562373095*dApardt[0]*q_)/m_; 
  alpha[1] = -(1.414213562373095*dApardt[1]*q_)/m_; 
  if (alpha0>0) { 
  double rVal[4];  // rVal=f1/f0 at each node 
  rVal[0] = -(1.0*(1.732050807568877*fl[3]-3.0*fl[2]))/(3.464101615137754*EPSILON-1.0*fl[1]+1.732050807568877*fl[0]); 
  rVal[1] = (1.732050807568877*fl[3]+3.0*fl[2])/(3.464101615137754*EPSILON+fl[1]+1.732050807568877*fl[0]); 
  rVal[2] = -(1.0*(1.732050807568877*fl[3]-3.0*fl[2]))/(3.464101615137754*EPSILON-1.0*fl[1]+1.732050807568877*fl[0]); 
  rVal[3] = (1.732050807568877*fl[3]+3.0*fl[2])/(3.464101615137754*EPSILON+fl[1]+1.732050807568877*fl[0]); 
  double fqVal[4];  // fqVal = anti-limited f evaluated at each node 
  double cflVal = std::abs((alpha0*dt)/dxv[1]); 
  fqVal[0] = -0.2886751345948129*(fl[1]-1.732050807568877*fl[0])*limTheta(rVal[0],1.0,cflVal); 
  fqVal[1] = 0.2886751345948129*(fl[1]+1.732050807568877*fl[0])*limTheta(rVal[1],1.0,cflVal); 
  fqVal[2] = -0.2886751345948129*(fl[1]-1.732050807568877*fl[0])*limTheta(rVal[2],1.0,cflVal); 
  fqVal[3] = 0.2886751345948129*(fl[1]+1.732050807568877*fl[0])*limTheta(rVal[3],1.0,cflVal); 
  double fhatALVal[4];  // fhatALVal = mode coefficients of anti-limited f 
  fhatALVal[0] = 0.5*(fqVal[3]+fqVal[2]+fqVal[1]+fqVal[0]); 
  fhatALVal[1] = 0.8660254037844386*(fqVal[3]-1.0*fqVal[2]+fqVal[1]-1.0*fqVal[0]); 
  fhatALVal[2] = 0.8660254037844386*(fqVal[3]+fqVal[2]-1.0*(fqVal[1]+fqVal[0])); 
  fhatALVal[3] = 1.5*(fqVal[3]-1.0*(fqVal[2]+fqVal[1])+fqVal[0]); 
  incr[0] = 0.25*(1.732050807568877*alpha[1]*fhatALVal[3]+1.732050807568877*alpha[0]*fhatALVal[2]+alpha[1]*fhatALVal[1]+alpha[0]*fhatALVal[0])*dfac_v; 
  incr[1] = 0.25*(1.732050807568877*alpha[0]*fhatALVal[3]+1.732050807568877*alpha[1]*fhatALVal[2]+alpha[0]*fhatALVal[1]+fhatALVal[0]*alpha[1])*dfac_v; 
  incr[2] = -0.25*(3.0*alpha[1]*fhatALVal[3]+3.0*alpha[0]*fhatALVal[2]+1.732050807568877*alpha[1]*fhatALVal[1]+1.732050807568877*alpha[0]*fhatALVal[0])*dfac_v; 
  incr[3] = -0.25*(3.0*alpha[0]*fhatALVal[3]+3.0*alpha[1]*fhatALVal[2]+1.732050807568877*alpha[0]*fhatALVal[1]+1.732050807568877*fhatALVal[0]*alpha[1])*dfac_v; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += -1.0*incr[1]; 
  outl[2] += incr[2]; 
  outl[3] += incr[3]; 
  } else { 
  double rVal[4];  // rVal=f1/f0 at each node 
  rVal[0] = -(1.0*(1.732050807568877*fr[3]-3.0*fr[2]))/(3.464101615137754*EPSILON-1.0*fr[1]+1.732050807568877*fr[0]); 
  rVal[1] = (1.732050807568877*fr[3]+3.0*fr[2])/(3.464101615137754*EPSILON+fr[1]+1.732050807568877*fr[0]); 
  rVal[2] = -(1.0*(1.732050807568877*fr[3]-3.0*fr[2]))/(3.464101615137754*EPSILON-1.0*fr[1]+1.732050807568877*fr[0]); 
  rVal[3] = (1.732050807568877*fr[3]+3.0*fr[2])/(3.464101615137754*EPSILON+fr[1]+1.732050807568877*fr[0]); 
  double fqVal[4];  // fqVal = anti-limited f evaluated at each node 
  double cflVal = std::abs((alpha0*dt)/dxv[1]); 
  fqVal[0] = -0.2886751345948129*(fr[1]-1.732050807568877*fr[0])*limTheta(rVal[0],-1.0,cflVal); 
  fqVal[1] = 0.2886751345948129*(fr[1]+1.732050807568877*fr[0])*limTheta(rVal[1],-1.0,cflVal); 
  fqVal[2] = -0.2886751345948129*(fr[1]-1.732050807568877*fr[0])*limTheta(rVal[2],-1.0,cflVal); 
  fqVal[3] = 0.2886751345948129*(fr[1]+1.732050807568877*fr[0])*limTheta(rVal[3],-1.0,cflVal); 
  double fhatALVal[4];  // fhatALVal = mode coefficients of anti-limited f 
  fhatALVal[0] = 0.5*(fqVal[3]+fqVal[2]+fqVal[1]+fqVal[0]); 
  fhatALVal[1] = 0.8660254037844386*(fqVal[3]-1.0*fqVal[2]+fqVal[1]-1.0*fqVal[0]); 
  fhatALVal[2] = 0.8660254037844386*(fqVal[3]+fqVal[2]-1.0*(fqVal[1]+fqVal[0])); 
  fhatALVal[3] = 1.5*(fqVal[3]-1.0*(fqVal[2]+fqVal[1])+fqVal[0]); 
  incr[0] = -0.25*(1.732050807568877*alpha[1]*fhatALVal[3]+1.732050807568877*alpha[0]*fhatALVal[2]-1.0*alpha[1]*fhatALVal[1]-1.0*alpha[0]*fhatALVal[0])*dfac_v; 
  incr[1] = -0.25*(1.732050807568877*alpha[0]*fhatALVal[3]+1.732050807568877*alpha[1]*fhatALVal[2]-1.0*alpha[0]*fhatALVal[1]-1.0*fhatALVal[0]*alpha[1])*dfac_v; 
  incr[2] = 0.25*(3.0*alpha[1]*fhatALVal[3]+3.0*alpha[0]*fhatALVal[2]-1.732050807568877*alpha[1]*fhatALVal[1]-1.732050807568877*alpha[0]*fhatALVal[0])*dfac_v; 
  incr[3] = 0.25*(3.0*alpha[0]*fhatALVal[3]+3.0*alpha[1]*fhatALVal[2]-1.732050807568877*alpha[0]*fhatALVal[1]-1.732050807568877*fhatALVal[0]*alpha[1])*dfac_v; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += -1.0*incr[1]; 
  outl[2] += incr[2]; 
  outl[3] += incr[3]; 
  } 
  return std::abs(alpha0); 
} 
double EmGyrokineticSurfPositivity1x1vSer_Vpar_P2_Bvars_0(const double q_, const double m_, const double dt, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *Apar, const double *dApardt, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[8]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = -(0.5*(2.449489742783178*Phi[1]*dfac_x+1.414213562373095*dApardt[0])*q_)/m_; 

  double amax = 0.0; 
  bool upwind = false; 
  if(upwind) 
    amax = std::abs(alpha0); 
  else 
    amax = amax_in; 

  double Ghat[8]; 
  Ghat[0] = (-(3.061862178478972*Phi[2]*fr[7]*dfac_x*q_)/m_)-(3.061862178478972*Phi[2]*fl[7]*dfac_x*q_)/m_-(1.369306393762915*Phi[1]*fr[5]*dfac_x*q_)/m_-(1.369306393762915*Phi[1]*fl[5]*dfac_x*q_)/m_+(2.371708245126284*Phi[2]*fr[3]*dfac_x*q_)/m_-(2.371708245126284*Phi[2]*fl[3]*dfac_x*q_)/m_+(1.060660171779821*Phi[1]*fr[2]*dfac_x*q_)/m_-(1.060660171779821*Phi[1]*fl[2]*dfac_x*q_)/m_-(1.369306393762915*fr[1]*Phi[2]*dfac_x*q_)/m_-(1.369306393762915*fl[1]*Phi[2]*dfac_x*q_)/m_-(0.6123724356957944*fr[0]*Phi[1]*dfac_x*q_)/m_-(0.6123724356957944*fl[0]*Phi[1]*dfac_x*q_)/m_-(0.7905694150420948*dApardt[1]*fr[7]*q_)/m_-(0.7905694150420948*dApardt[1]*fl[7]*q_)/m_+(0.6123724356957944*dApardt[2]*fr[6]*q_)/m_-(0.6123724356957944*dApardt[2]*fl[6]*q_)/m_-(0.7905694150420947*dApardt[0]*fr[5]*q_)/m_-(0.7905694150420947*dApardt[0]*fl[5]*q_)/m_-(0.3535533905932737*dApardt[2]*fr[4]*q_)/m_-(0.3535533905932737*dApardt[2]*fl[4]*q_)/m_+(0.6123724356957944*dApardt[1]*fr[3]*q_)/m_-(0.6123724356957944*dApardt[1]*fl[3]*q_)/m_+(0.6123724356957944*dApardt[0]*fr[2]*q_)/m_-(0.6123724356957944*dApardt[0]*fl[2]*q_)/m_-(0.3535533905932737*dApardt[1]*fr[1]*q_)/m_-(0.3535533905932737*dApardt[1]*fl[1]*q_)/m_-(0.3535533905932737*dApardt[0]*fr[0]*q_)/m_-(0.3535533905932737*dApardt[0]*fl[0]*q_)/m_-1.118033988749895*fr[5]*amax+1.118033988749895*fl[5]*amax+0.8660254037844386*fr[2]*amax+0.8660254037844386*fl[2]*amax-0.5*fr[0]*amax+0.5*fl[0]*amax; 
  Ghat[1] = (-(1.369306393762915*Phi[1]*fr[7]*dfac_x*q_)/m_)-(1.369306393762915*Phi[1]*fl[7]*dfac_x*q_)/m_+(2.121320343559642*Phi[2]*fr[6]*dfac_x*q_)/m_-(2.121320343559642*Phi[2]*fl[6]*dfac_x*q_)/m_-(3.061862178478972*Phi[2]*fr[5]*dfac_x*q_)/m_-(3.061862178478972*Phi[2]*fl[5]*dfac_x*q_)/m_-(1.224744871391589*Phi[2]*fr[4]*dfac_x*q_)/m_-(1.224744871391589*Phi[2]*fl[4]*dfac_x*q_)/m_+(1.060660171779821*Phi[1]*fr[3]*dfac_x*q_)/m_-(1.060660171779821*Phi[1]*fl[3]*dfac_x*q_)/m_+(2.371708245126284*Phi[2]*fr[2]*dfac_x*q_)/m_-(2.371708245126284*Phi[2]*fl[2]*dfac_x*q_)/m_-(1.369306393762915*fr[0]*Phi[2]*dfac_x*q_)/m_-(1.369306393762915*fl[0]*Phi[2]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fr[1]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fl[1]*dfac_x*q_)/m_-(0.7071067811865475*dApardt[2]*fr[7]*q_)/m_-(0.7905694150420948*dApardt[0]*fr[7]*q_)/m_-(0.7071067811865475*dApardt[2]*fl[7]*q_)/m_-(0.7905694150420948*dApardt[0]*fl[7]*q_)/m_+(0.5477225575051661*dApardt[1]*fr[6]*q_)/m_-(0.5477225575051661*dApardt[1]*fl[6]*q_)/m_-(0.7905694150420947*dApardt[1]*fr[5]*q_)/m_-(0.7905694150420947*dApardt[1]*fl[5]*q_)/m_-(0.3162277660168379*dApardt[1]*fr[4]*q_)/m_-(0.3162277660168379*dApardt[1]*fl[4]*q_)/m_+(0.5477225575051661*dApardt[2]*fr[3]*q_)/m_+(0.6123724356957944*dApardt[0]*fr[3]*q_)/m_-(0.5477225575051661*dApardt[2]*fl[3]*q_)/m_-(0.6123724356957944*dApardt[0]*fl[3]*q_)/m_+(0.6123724356957944*dApardt[1]*fr[2]*q_)/m_-(0.6123724356957944*dApardt[1]*fl[2]*q_)/m_-(0.3162277660168379*fr[1]*dApardt[2]*q_)/m_-(0.3162277660168379*fl[1]*dApardt[2]*q_)/m_-(0.3535533905932737*dApardt[0]*fr[1]*q_)/m_-(0.3535533905932737*dApardt[0]*fl[1]*q_)/m_-(0.3535533905932737*fr[0]*dApardt[1]*q_)/m_-(0.3535533905932737*fl[0]*dApardt[1]*q_)/m_-1.118033988749895*fr[7]*amax+1.118033988749895*fl[7]*amax+0.8660254037844386*fr[3]*amax+0.8660254037844386*fl[3]*amax-0.5*fr[1]*amax+0.5*fl[1]*amax; 
  Ghat[4] = (-(2.738612787525831*Phi[2]*fr[7]*dfac_x*q_)/m_)-(2.738612787525831*Phi[2]*fl[7]*dfac_x*q_)/m_+(1.060660171779821*Phi[1]*fr[6]*dfac_x*q_)/m_-(1.060660171779821*Phi[1]*fl[6]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fr[4]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fl[4]*dfac_x*q_)/m_+(2.121320343559642*Phi[2]*fr[3]*dfac_x*q_)/m_-(2.121320343559642*Phi[2]*fl[3]*dfac_x*q_)/m_-(1.224744871391589*fr[1]*Phi[2]*dfac_x*q_)/m_-(1.224744871391589*fl[1]*Phi[2]*dfac_x*q_)/m_-(0.7071067811865475*dApardt[1]*fr[7]*q_)/m_-(0.7071067811865475*dApardt[1]*fl[7]*q_)/m_+(0.3912303982179757*dApardt[2]*fr[6]*q_)/m_+(0.6123724356957944*dApardt[0]*fr[6]*q_)/m_-(0.3912303982179757*dApardt[2]*fl[6]*q_)/m_-(0.6123724356957944*dApardt[0]*fl[6]*q_)/m_-(0.7905694150420947*dApardt[2]*fr[5]*q_)/m_-(0.7905694150420947*dApardt[2]*fl[5]*q_)/m_-(0.2258769757263128*dApardt[2]*fr[4]*q_)/m_-(0.3535533905932737*dApardt[0]*fr[4]*q_)/m_-(0.2258769757263128*dApardt[2]*fl[4]*q_)/m_-(0.3535533905932737*dApardt[0]*fl[4]*q_)/m_+(0.5477225575051661*dApardt[1]*fr[3]*q_)/m_-(0.5477225575051661*dApardt[1]*fl[3]*q_)/m_+(0.6123724356957944*dApardt[2]*fr[2]*q_)/m_-(0.6123724356957944*dApardt[2]*fl[2]*q_)/m_-(0.3535533905932737*fr[0]*dApardt[2]*q_)/m_-(0.3535533905932737*fl[0]*dApardt[2]*q_)/m_-(0.3162277660168379*dApardt[1]*fr[1]*q_)/m_-(0.3162277660168379*dApardt[1]*fl[1]*q_)/m_+0.8660254037844387*fr[6]*amax+0.8660254037844387*fl[6]*amax-0.5*fr[4]*amax+0.5*fl[4]*amax; 

  incr[0] = 0.5*Ghat[0]*dfac_v; 
  incr[1] = 0.5*Ghat[1]*dfac_v; 
  incr[2] = -0.8660254037844386*Ghat[0]*dfac_v; 
  incr[3] = -0.8660254037844386*Ghat[1]*dfac_v; 
  incr[4] = 0.5*Ghat[4]*dfac_v; 
  incr[5] = 1.118033988749895*Ghat[0]*dfac_v; 
  incr[6] = -0.8660254037844387*Ghat[4]*dfac_v; 
  incr[7] = 1.118033988749895*Ghat[1]*dfac_v; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 
  outr[4] += incr[4]; 
  outr[5] += incr[5]; 
  outr[6] += incr[6]; 
  outr[7] += incr[7]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += -1.0*incr[1]; 
  outl[2] += incr[2]; 
  outl[3] += incr[3]; 
  outl[4] += -1.0*incr[4]; 
  outl[5] += -1.0*incr[5]; 
  outl[6] += incr[6]; 
  outl[7] += -1.0*incr[7]; 
return std::abs(alpha0); 
} 
double GyrokineticSurfPositivity1x1vSer_X_P1_Bvars_1(const double q_, const double m_, const double dt, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[4]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = wv; 

  double alpha[4]; 
  alpha[0] = 2.0*wv; 
  if (alpha0>0) { 
  double rVal[4];  // rVal=f1/f0 at each node 
  rVal[0] = -(1.0*(1.732050807568877*fl[3]-3.0*fl[1]))/(3.464101615137754*EPSILON-1.0*fl[2]+1.732050807568877*fl[0]); 
  rVal[1] = -(1.0*(1.732050807568877*fl[3]-3.0*fl[1]))/(3.464101615137754*EPSILON-1.0*fl[2]+1.732050807568877*fl[0]); 
  rVal[2] = (1.732050807568877*fl[3]+3.0*fl[1])/(3.464101615137754*EPSILON+fl[2]+1.732050807568877*fl[0]); 
  rVal[3] = (1.732050807568877*fl[3]+3.0*fl[1])/(3.464101615137754*EPSILON+fl[2]+1.732050807568877*fl[0]); 
  double fqVal[4];  // fqVal = anti-limited f evaluated at each node 
  double cflVal = std::abs((alpha0*dt)/dxv[0]); 
  fqVal[0] = -0.2886751345948129*(fl[2]-1.732050807568877*fl[0])*limTheta(rVal[0],1.0,cflVal); 
  fqVal[1] = -0.2886751345948129*(fl[2]-1.732050807568877*fl[0])*limTheta(rVal[1],1.0,cflVal); 
  fqVal[2] = 0.2886751345948129*(fl[2]+1.732050807568877*fl[0])*limTheta(rVal[2],1.0,cflVal); 
  fqVal[3] = 0.2886751345948129*(fl[2]+1.732050807568877*fl[0])*limTheta(rVal[3],1.0,cflVal); 
  double fhatALVal[4];  // fhatALVal = mode coefficients of anti-limited f 
  fhatALVal[0] = 0.5*(fqVal[3]+fqVal[2]+fqVal[1]+fqVal[0]); 
  fhatALVal[1] = 0.8660254037844386*(fqVal[3]-1.0*fqVal[2]+fqVal[1]-1.0*fqVal[0]); 
  fhatALVal[2] = 0.8660254037844386*(fqVal[3]+fqVal[2]-1.0*(fqVal[1]+fqVal[0])); 
  fhatALVal[3] = 1.5*(fqVal[3]-1.0*(fqVal[2]+fqVal[1])+fqVal[0]); 
  incr[0] = 0.25*alpha[0]*(1.732050807568877*fhatALVal[1]+fhatALVal[0])*dfac_x; 
  incr[1] = -0.25*alpha[0]*(3.0*fhatALVal[1]+1.732050807568877*fhatALVal[0])*dfac_x; 
  incr[2] = 0.25*alpha[0]*(1.732050807568877*fhatALVal[3]+fhatALVal[2])*dfac_x; 
  incr[3] = -0.25*alpha[0]*(3.0*fhatALVal[3]+1.732050807568877*fhatALVal[2])*dfac_x; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += incr[1]; 
  outl[2] += -1.0*incr[2]; 
  outl[3] += incr[3]; 
  } else { 
  double rVal[4];  // rVal=f1/f0 at each node 
  rVal[0] = -(1.0*(1.732050807568877*fr[3]-3.0*fr[1]))/(3.464101615137754*EPSILON-1.0*fr[2]+1.732050807568877*fr[0]); 
  rVal[1] = -(1.0*(1.732050807568877*fr[3]-3.0*fr[1]))/(3.464101615137754*EPSILON-1.0*fr[2]+1.732050807568877*fr[0]); 
  rVal[2] = (1.732050807568877*fr[3]+3.0*fr[1])/(3.464101615137754*EPSILON+fr[2]+1.732050807568877*fr[0]); 
  rVal[3] = (1.732050807568877*fr[3]+3.0*fr[1])/(3.464101615137754*EPSILON+fr[2]+1.732050807568877*fr[0]); 
  double fqVal[4];  // fqVal = anti-limited f evaluated at each node 
  double cflVal = std::abs((alpha0*dt)/dxv[0]); 
  fqVal[0] = -0.2886751345948129*(fr[2]-1.732050807568877*fr[0])*limTheta(rVal[0],-1.0,cflVal); 
  fqVal[1] = -0.2886751345948129*(fr[2]-1.732050807568877*fr[0])*limTheta(rVal[1],-1.0,cflVal); 
  fqVal[2] = 0.2886751345948129*(fr[2]+1.732050807568877*fr[0])*limTheta(rVal[2],-1.0,cflVal); 
  fqVal[3] = 0.2886751345948129*(fr[2]+1.732050807568877*fr[0])*limTheta(rVal[3],-1.0,cflVal); 
  double fhatALVal[4];  // fhatALVal = mode coefficients of anti-limited f 
  fhatALVal[0] = 0.5*(fqVal[3]+fqVal[2]+fqVal[1]+fqVal[0]); 
  fhatALVal[1] = 0.8660254037844386*(fqVal[3]-1.0*fqVal[2]+fqVal[1]-1.0*fqVal[0]); 
  fhatALVal[2] = 0.8660254037844386*(fqVal[3]+fqVal[2]-1.0*(fqVal[1]+fqVal[0])); 
  fhatALVal[3] = 1.5*(fqVal[3]-1.0*(fqVal[2]+fqVal[1])+fqVal[0]); 
  incr[0] = -0.25*alpha[0]*(1.732050807568877*fhatALVal[1]-1.0*fhatALVal[0])*dfac_x; 
  incr[1] = 0.25*alpha[0]*(3.0*fhatALVal[1]-1.732050807568877*fhatALVal[0])*dfac_x; 
  incr[2] = -0.25*alpha[0]*(1.732050807568877*fhatALVal[3]-1.0*fhatALVal[2])*dfac_x; 
  incr[3] = 0.25*alpha[0]*(3.0*fhatALVal[3]-1.732050807568877*fhatALVal[2])*dfac_x; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += incr[1]; 
  outl[2] += -1.0*incr[2]; 
  outl[3] += incr[3]; 
  } 
  return std::abs(alpha0); 
} 
double GyrokineticSurfPositivity1x1vSer_X_P2_Bvars_1(const double q_, const double m_, const double dt, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[8]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = wv; 

  double amax = 0.0; 
  bool upwind = false; 
  if(upwind) 
    amax = std::abs(alpha0); 
  else 
    amax = amax_in; 

  double Ghat[8]; 
  Ghat[0] = 1.118033988749895*fr[4]*wv+1.118033988749895*fl[4]*wv-0.8660254037844386*fr[1]*wv+0.8660254037844386*fl[1]*wv+0.5*fr[0]*wv+0.5*fl[0]*wv+(0.6454972243679028*fr[6])/dfac_v+(0.6454972243679028*fl[6])/dfac_v-(0.5*fr[3])/dfac_v+(0.5*fl[3])/dfac_v+(0.2886751345948129*fr[2])/dfac_v+(0.2886751345948129*fl[2])/dfac_v-1.118033988749895*fr[4]*amax+1.118033988749895*fl[4]*amax+0.8660254037844386*fr[1]*amax+0.8660254037844386*fl[1]*amax-0.5*fr[0]*amax+0.5*fl[0]*amax; 
  Ghat[2] = 1.118033988749895*fr[6]*wv+1.118033988749895*fl[6]*wv-0.8660254037844386*fr[3]*wv+0.8660254037844386*fl[3]*wv+0.5*fr[2]*wv+0.5*fl[2]*wv-(0.447213595499958*fr[7])/dfac_v+(0.447213595499958*fl[7])/dfac_v+(0.2581988897471612*fr[5])/dfac_v+(0.2581988897471612*fl[5])/dfac_v+(0.6454972243679029*fr[4])/dfac_v+(0.6454972243679029*fl[4])/dfac_v-(0.5*fr[1])/dfac_v+(0.5*fl[1])/dfac_v+(0.2886751345948129*fr[0])/dfac_v+(0.2886751345948129*fl[0])/dfac_v-1.118033988749895*fr[6]*amax+1.118033988749895*fl[6]*amax+0.8660254037844386*fr[3]*amax+0.8660254037844386*fl[3]*amax-0.5*fr[2]*amax+0.5*fl[2]*amax; 
  Ghat[5] = (-0.8660254037844387*fr[7]*wv)+0.8660254037844387*fl[7]*wv+0.5*fr[5]*wv+0.5*fl[5]*wv+(0.5773502691896257*fr[6])/dfac_v+(0.5773502691896257*fl[6])/dfac_v-(0.4472135954999579*fr[3])/dfac_v+(0.4472135954999579*fl[3])/dfac_v+(0.2581988897471612*fr[2])/dfac_v+(0.2581988897471612*fl[2])/dfac_v+0.8660254037844387*fr[7]*amax+0.8660254037844387*fl[7]*amax-0.5*fr[5]*amax+0.5*fl[5]*amax; 

  incr[0] = 0.5*Ghat[0]*dfac_x; 
  incr[1] = -0.8660254037844386*Ghat[0]*dfac_x; 
  incr[2] = 0.5*Ghat[2]*dfac_x; 
  incr[3] = -0.8660254037844386*Ghat[2]*dfac_x; 
  incr[4] = 1.118033988749895*Ghat[0]*dfac_x; 
  incr[5] = 0.5*Ghat[5]*dfac_x; 
  incr[6] = 1.118033988749895*Ghat[2]*dfac_x; 
  incr[7] = -0.8660254037844387*Ghat[5]*dfac_x; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 
  outr[4] += incr[4]; 
  outr[5] += incr[5]; 
  outr[6] += incr[6]; 
  outr[7] += incr[7]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += incr[1]; 
  outl[2] += -1.0*incr[2]; 
  outl[3] += incr[3]; 
  outl[4] += -1.0*incr[4]; 
  outl[5] += -1.0*incr[5]; 
  outl[6] += -1.0*incr[6]; 
  outl[7] += incr[7]; 
return std::abs(alpha0); 
} 
double GyrokineticSurfPositivity1x1vSer_Vpar_P1_Bvars_1(const double q_, const double m_, const double dt, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[4]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = -(1.224744871391589*Phi[1]*dfac_x*q_)/m_; 

  double alpha[4]; 
  alpha[0] = -(2.449489742783178*Phi[1]*dfac_x*q_)/m_; 
  if (alpha0>0) { 
  double rVal[4];  // rVal=f1/f0 at each node 
  rVal[0] = -(1.0*(1.732050807568877*fl[3]-3.0*fl[2]))/(3.464101615137754*EPSILON-1.0*fl[1]+1.732050807568877*fl[0]); 
  rVal[1] = (1.732050807568877*fl[3]+3.0*fl[2])/(3.464101615137754*EPSILON+fl[1]+1.732050807568877*fl[0]); 
  rVal[2] = -(1.0*(1.732050807568877*fl[3]-3.0*fl[2]))/(3.464101615137754*EPSILON-1.0*fl[1]+1.732050807568877*fl[0]); 
  rVal[3] = (1.732050807568877*fl[3]+3.0*fl[2])/(3.464101615137754*EPSILON+fl[1]+1.732050807568877*fl[0]); 
  double fqVal[4];  // fqVal = anti-limited f evaluated at each node 
  double cflVal = std::abs((alpha0*dt)/dxv[1]); 
  fqVal[0] = -0.2886751345948129*(fl[1]-1.732050807568877*fl[0])*limTheta(rVal[0],1.0,cflVal); 
  fqVal[1] = 0.2886751345948129*(fl[1]+1.732050807568877*fl[0])*limTheta(rVal[1],1.0,cflVal); 
  fqVal[2] = -0.2886751345948129*(fl[1]-1.732050807568877*fl[0])*limTheta(rVal[2],1.0,cflVal); 
  fqVal[3] = 0.2886751345948129*(fl[1]+1.732050807568877*fl[0])*limTheta(rVal[3],1.0,cflVal); 
  double fhatALVal[4];  // fhatALVal = mode coefficients of anti-limited f 
  fhatALVal[0] = 0.5*(fqVal[3]+fqVal[2]+fqVal[1]+fqVal[0]); 
  fhatALVal[1] = 0.8660254037844386*(fqVal[3]-1.0*fqVal[2]+fqVal[1]-1.0*fqVal[0]); 
  fhatALVal[2] = 0.8660254037844386*(fqVal[3]+fqVal[2]-1.0*(fqVal[1]+fqVal[0])); 
  fhatALVal[3] = 1.5*(fqVal[3]-1.0*(fqVal[2]+fqVal[1])+fqVal[0]); 
  incr[0] = 0.25*alpha[0]*(1.732050807568877*fhatALVal[2]+fhatALVal[0])*dfac_v; 
  incr[1] = 0.25*alpha[0]*(1.732050807568877*fhatALVal[3]+fhatALVal[1])*dfac_v; 
  incr[2] = -0.25*alpha[0]*(3.0*fhatALVal[2]+1.732050807568877*fhatALVal[0])*dfac_v; 
  incr[3] = -0.25*alpha[0]*(3.0*fhatALVal[3]+1.732050807568877*fhatALVal[1])*dfac_v; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += -1.0*incr[1]; 
  outl[2] += incr[2]; 
  outl[3] += incr[3]; 
  } else { 
  double rVal[4];  // rVal=f1/f0 at each node 
  rVal[0] = -(1.0*(1.732050807568877*fr[3]-3.0*fr[2]))/(3.464101615137754*EPSILON-1.0*fr[1]+1.732050807568877*fr[0]); 
  rVal[1] = (1.732050807568877*fr[3]+3.0*fr[2])/(3.464101615137754*EPSILON+fr[1]+1.732050807568877*fr[0]); 
  rVal[2] = -(1.0*(1.732050807568877*fr[3]-3.0*fr[2]))/(3.464101615137754*EPSILON-1.0*fr[1]+1.732050807568877*fr[0]); 
  rVal[3] = (1.732050807568877*fr[3]+3.0*fr[2])/(3.464101615137754*EPSILON+fr[1]+1.732050807568877*fr[0]); 
  double fqVal[4];  // fqVal = anti-limited f evaluated at each node 
  double cflVal = std::abs((alpha0*dt)/dxv[1]); 
  fqVal[0] = -0.2886751345948129*(fr[1]-1.732050807568877*fr[0])*limTheta(rVal[0],-1.0,cflVal); 
  fqVal[1] = 0.2886751345948129*(fr[1]+1.732050807568877*fr[0])*limTheta(rVal[1],-1.0,cflVal); 
  fqVal[2] = -0.2886751345948129*(fr[1]-1.732050807568877*fr[0])*limTheta(rVal[2],-1.0,cflVal); 
  fqVal[3] = 0.2886751345948129*(fr[1]+1.732050807568877*fr[0])*limTheta(rVal[3],-1.0,cflVal); 
  double fhatALVal[4];  // fhatALVal = mode coefficients of anti-limited f 
  fhatALVal[0] = 0.5*(fqVal[3]+fqVal[2]+fqVal[1]+fqVal[0]); 
  fhatALVal[1] = 0.8660254037844386*(fqVal[3]-1.0*fqVal[2]+fqVal[1]-1.0*fqVal[0]); 
  fhatALVal[2] = 0.8660254037844386*(fqVal[3]+fqVal[2]-1.0*(fqVal[1]+fqVal[0])); 
  fhatALVal[3] = 1.5*(fqVal[3]-1.0*(fqVal[2]+fqVal[1])+fqVal[0]); 
  incr[0] = -0.25*alpha[0]*(1.732050807568877*fhatALVal[2]-1.0*fhatALVal[0])*dfac_v; 
  incr[1] = -0.25*alpha[0]*(1.732050807568877*fhatALVal[3]-1.0*fhatALVal[1])*dfac_v; 
  incr[2] = 0.25*alpha[0]*(3.0*fhatALVal[2]-1.732050807568877*fhatALVal[0])*dfac_v; 
  incr[3] = 0.25*alpha[0]*(3.0*fhatALVal[3]-1.732050807568877*fhatALVal[1])*dfac_v; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += -1.0*incr[1]; 
  outl[2] += incr[2]; 
  outl[3] += incr[3]; 
  } 
  return std::abs(alpha0); 
} 
double GyrokineticSurfPositivity1x1vSer_Vpar_P2_Bvars_1(const double q_, const double m_, const double dt, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[8]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = -(1.224744871391589*Phi[1]*dfac_x*q_)/m_; 

  double amax = 0.0; 
  bool upwind = false; 
  if(upwind) 
    amax = std::abs(alpha0); 
  else 
    amax = amax_in; 

  double Ghat[8]; 
  Ghat[0] = (-(3.061862178478972*Phi[2]*fr[7]*dfac_x*q_)/m_)-(3.061862178478972*Phi[2]*fl[7]*dfac_x*q_)/m_-(1.369306393762915*Phi[1]*fr[5]*dfac_x*q_)/m_-(1.369306393762915*Phi[1]*fl[5]*dfac_x*q_)/m_+(2.371708245126284*Phi[2]*fr[3]*dfac_x*q_)/m_-(2.371708245126284*Phi[2]*fl[3]*dfac_x*q_)/m_+(1.060660171779821*Phi[1]*fr[2]*dfac_x*q_)/m_-(1.060660171779821*Phi[1]*fl[2]*dfac_x*q_)/m_-(1.369306393762915*fr[1]*Phi[2]*dfac_x*q_)/m_-(1.369306393762915*fl[1]*Phi[2]*dfac_x*q_)/m_-(0.6123724356957944*fr[0]*Phi[1]*dfac_x*q_)/m_-(0.6123724356957944*fl[0]*Phi[1]*dfac_x*q_)/m_-1.118033988749895*fr[5]*amax+1.118033988749895*fl[5]*amax+0.8660254037844386*fr[2]*amax+0.8660254037844386*fl[2]*amax-0.5*fr[0]*amax+0.5*fl[0]*amax; 
  Ghat[1] = (-(1.369306393762915*Phi[1]*fr[7]*dfac_x*q_)/m_)-(1.369306393762915*Phi[1]*fl[7]*dfac_x*q_)/m_+(2.121320343559642*Phi[2]*fr[6]*dfac_x*q_)/m_-(2.121320343559642*Phi[2]*fl[6]*dfac_x*q_)/m_-(3.061862178478972*Phi[2]*fr[5]*dfac_x*q_)/m_-(3.061862178478972*Phi[2]*fl[5]*dfac_x*q_)/m_-(1.224744871391589*Phi[2]*fr[4]*dfac_x*q_)/m_-(1.224744871391589*Phi[2]*fl[4]*dfac_x*q_)/m_+(1.060660171779821*Phi[1]*fr[3]*dfac_x*q_)/m_-(1.060660171779821*Phi[1]*fl[3]*dfac_x*q_)/m_+(2.371708245126284*Phi[2]*fr[2]*dfac_x*q_)/m_-(2.371708245126284*Phi[2]*fl[2]*dfac_x*q_)/m_-(1.369306393762915*fr[0]*Phi[2]*dfac_x*q_)/m_-(1.369306393762915*fl[0]*Phi[2]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fr[1]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fl[1]*dfac_x*q_)/m_-1.118033988749895*fr[7]*amax+1.118033988749895*fl[7]*amax+0.8660254037844386*fr[3]*amax+0.8660254037844386*fl[3]*amax-0.5*fr[1]*amax+0.5*fl[1]*amax; 
  Ghat[4] = (-(2.738612787525831*Phi[2]*fr[7]*dfac_x*q_)/m_)-(2.738612787525831*Phi[2]*fl[7]*dfac_x*q_)/m_+(1.060660171779821*Phi[1]*fr[6]*dfac_x*q_)/m_-(1.060660171779821*Phi[1]*fl[6]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fr[4]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fl[4]*dfac_x*q_)/m_+(2.121320343559642*Phi[2]*fr[3]*dfac_x*q_)/m_-(2.121320343559642*Phi[2]*fl[3]*dfac_x*q_)/m_-(1.224744871391589*fr[1]*Phi[2]*dfac_x*q_)/m_-(1.224744871391589*fl[1]*Phi[2]*dfac_x*q_)/m_+0.8660254037844387*fr[6]*amax+0.8660254037844387*fl[6]*amax-0.5*fr[4]*amax+0.5*fl[4]*amax; 

  incr[0] = 0.5*Ghat[0]*dfac_v; 
  incr[1] = 0.5*Ghat[1]*dfac_v; 
  incr[2] = -0.8660254037844386*Ghat[0]*dfac_v; 
  incr[3] = -0.8660254037844386*Ghat[1]*dfac_v; 
  incr[4] = 0.5*Ghat[4]*dfac_v; 
  incr[5] = 1.118033988749895*Ghat[0]*dfac_v; 
  incr[6] = -0.8660254037844387*Ghat[4]*dfac_v; 
  incr[7] = 1.118033988749895*Ghat[1]*dfac_v; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 
  outr[4] += incr[4]; 
  outr[5] += incr[5]; 
  outr[6] += incr[6]; 
  outr[7] += incr[7]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += -1.0*incr[1]; 
  outl[2] += incr[2]; 
  outl[3] += incr[3]; 
  outl[4] += -1.0*incr[4]; 
  outl[5] += -1.0*incr[5]; 
  outl[6] += incr[6]; 
  outl[7] += -1.0*incr[7]; 
return std::abs(alpha0); 
} 
double EmGyrokineticSurfPositivity1x1vSer_X_P1_Bvars_1(const double q_, const double m_, const double dt, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *Apar, const double *dApardt, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[4]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = wv; 

  double alpha[4]; 
  alpha[0] = 2.0*wv; 
  if (alpha0>0) { 
  double rVal[4];  // rVal=f1/f0 at each node 
  rVal[0] = -(1.0*(1.732050807568877*fl[3]-3.0*fl[1]))/(3.464101615137754*EPSILON-1.0*fl[2]+1.732050807568877*fl[0]); 
  rVal[1] = -(1.0*(1.732050807568877*fl[3]-3.0*fl[1]))/(3.464101615137754*EPSILON-1.0*fl[2]+1.732050807568877*fl[0]); 
  rVal[2] = (1.732050807568877*fl[3]+3.0*fl[1])/(3.464101615137754*EPSILON+fl[2]+1.732050807568877*fl[0]); 
  rVal[3] = (1.732050807568877*fl[3]+3.0*fl[1])/(3.464101615137754*EPSILON+fl[2]+1.732050807568877*fl[0]); 
  double fqVal[4];  // fqVal = anti-limited f evaluated at each node 
  double cflVal = std::abs((alpha0*dt)/dxv[0]); 
  fqVal[0] = -0.2886751345948129*(fl[2]-1.732050807568877*fl[0])*limTheta(rVal[0],1.0,cflVal); 
  fqVal[1] = -0.2886751345948129*(fl[2]-1.732050807568877*fl[0])*limTheta(rVal[1],1.0,cflVal); 
  fqVal[2] = 0.2886751345948129*(fl[2]+1.732050807568877*fl[0])*limTheta(rVal[2],1.0,cflVal); 
  fqVal[3] = 0.2886751345948129*(fl[2]+1.732050807568877*fl[0])*limTheta(rVal[3],1.0,cflVal); 
  double fhatALVal[4];  // fhatALVal = mode coefficients of anti-limited f 
  fhatALVal[0] = 0.5*(fqVal[3]+fqVal[2]+fqVal[1]+fqVal[0]); 
  fhatALVal[1] = 0.8660254037844386*(fqVal[3]-1.0*fqVal[2]+fqVal[1]-1.0*fqVal[0]); 
  fhatALVal[2] = 0.8660254037844386*(fqVal[3]+fqVal[2]-1.0*(fqVal[1]+fqVal[0])); 
  fhatALVal[3] = 1.5*(fqVal[3]-1.0*(fqVal[2]+fqVal[1])+fqVal[0]); 
  incr[0] = 0.25*alpha[0]*(1.732050807568877*fhatALVal[1]+fhatALVal[0])*dfac_x; 
  incr[1] = -0.25*alpha[0]*(3.0*fhatALVal[1]+1.732050807568877*fhatALVal[0])*dfac_x; 
  incr[2] = 0.25*alpha[0]*(1.732050807568877*fhatALVal[3]+fhatALVal[2])*dfac_x; 
  incr[3] = -0.25*alpha[0]*(3.0*fhatALVal[3]+1.732050807568877*fhatALVal[2])*dfac_x; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += incr[1]; 
  outl[2] += -1.0*incr[2]; 
  outl[3] += incr[3]; 
  } else { 
  double rVal[4];  // rVal=f1/f0 at each node 
  rVal[0] = -(1.0*(1.732050807568877*fr[3]-3.0*fr[1]))/(3.464101615137754*EPSILON-1.0*fr[2]+1.732050807568877*fr[0]); 
  rVal[1] = -(1.0*(1.732050807568877*fr[3]-3.0*fr[1]))/(3.464101615137754*EPSILON-1.0*fr[2]+1.732050807568877*fr[0]); 
  rVal[2] = (1.732050807568877*fr[3]+3.0*fr[1])/(3.464101615137754*EPSILON+fr[2]+1.732050807568877*fr[0]); 
  rVal[3] = (1.732050807568877*fr[3]+3.0*fr[1])/(3.464101615137754*EPSILON+fr[2]+1.732050807568877*fr[0]); 
  double fqVal[4];  // fqVal = anti-limited f evaluated at each node 
  double cflVal = std::abs((alpha0*dt)/dxv[0]); 
  fqVal[0] = -0.2886751345948129*(fr[2]-1.732050807568877*fr[0])*limTheta(rVal[0],-1.0,cflVal); 
  fqVal[1] = -0.2886751345948129*(fr[2]-1.732050807568877*fr[0])*limTheta(rVal[1],-1.0,cflVal); 
  fqVal[2] = 0.2886751345948129*(fr[2]+1.732050807568877*fr[0])*limTheta(rVal[2],-1.0,cflVal); 
  fqVal[3] = 0.2886751345948129*(fr[2]+1.732050807568877*fr[0])*limTheta(rVal[3],-1.0,cflVal); 
  double fhatALVal[4];  // fhatALVal = mode coefficients of anti-limited f 
  fhatALVal[0] = 0.5*(fqVal[3]+fqVal[2]+fqVal[1]+fqVal[0]); 
  fhatALVal[1] = 0.8660254037844386*(fqVal[3]-1.0*fqVal[2]+fqVal[1]-1.0*fqVal[0]); 
  fhatALVal[2] = 0.8660254037844386*(fqVal[3]+fqVal[2]-1.0*(fqVal[1]+fqVal[0])); 
  fhatALVal[3] = 1.5*(fqVal[3]-1.0*(fqVal[2]+fqVal[1])+fqVal[0]); 
  incr[0] = -0.25*alpha[0]*(1.732050807568877*fhatALVal[1]-1.0*fhatALVal[0])*dfac_x; 
  incr[1] = 0.25*alpha[0]*(3.0*fhatALVal[1]-1.732050807568877*fhatALVal[0])*dfac_x; 
  incr[2] = -0.25*alpha[0]*(1.732050807568877*fhatALVal[3]-1.0*fhatALVal[2])*dfac_x; 
  incr[3] = 0.25*alpha[0]*(3.0*fhatALVal[3]-1.732050807568877*fhatALVal[2])*dfac_x; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += incr[1]; 
  outl[2] += -1.0*incr[2]; 
  outl[3] += incr[3]; 
  } 
  return std::abs(alpha0); 
} 
double EmGyrokineticSurfPositivity1x1vSer_X_P2_Bvars_1(const double q_, const double m_, const double dt, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *Apar, const double *dApardt, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[8]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = wv; 

  double amax = 0.0; 
  bool upwind = false; 
  if(upwind) 
    amax = std::abs(alpha0); 
  else 
    amax = amax_in; 

  double Ghat[8]; 
  Ghat[0] = 1.118033988749895*fr[4]*wv+1.118033988749895*fl[4]*wv-0.8660254037844386*fr[1]*wv+0.8660254037844386*fl[1]*wv+0.5*fr[0]*wv+0.5*fl[0]*wv+(0.6454972243679028*fr[6])/dfac_v+(0.6454972243679028*fl[6])/dfac_v-(0.5*fr[3])/dfac_v+(0.5*fl[3])/dfac_v+(0.2886751345948129*fr[2])/dfac_v+(0.2886751345948129*fl[2])/dfac_v-1.118033988749895*fr[4]*amax+1.118033988749895*fl[4]*amax+0.8660254037844386*fr[1]*amax+0.8660254037844386*fl[1]*amax-0.5*fr[0]*amax+0.5*fl[0]*amax; 
  Ghat[2] = 1.118033988749895*fr[6]*wv+1.118033988749895*fl[6]*wv-0.8660254037844386*fr[3]*wv+0.8660254037844386*fl[3]*wv+0.5*fr[2]*wv+0.5*fl[2]*wv-(0.447213595499958*fr[7])/dfac_v+(0.447213595499958*fl[7])/dfac_v+(0.2581988897471612*fr[5])/dfac_v+(0.2581988897471612*fl[5])/dfac_v+(0.6454972243679029*fr[4])/dfac_v+(0.6454972243679029*fl[4])/dfac_v-(0.5*fr[1])/dfac_v+(0.5*fl[1])/dfac_v+(0.2886751345948129*fr[0])/dfac_v+(0.2886751345948129*fl[0])/dfac_v-1.118033988749895*fr[6]*amax+1.118033988749895*fl[6]*amax+0.8660254037844386*fr[3]*amax+0.8660254037844386*fl[3]*amax-0.5*fr[2]*amax+0.5*fl[2]*amax; 
  Ghat[5] = (-0.8660254037844387*fr[7]*wv)+0.8660254037844387*fl[7]*wv+0.5*fr[5]*wv+0.5*fl[5]*wv+(0.5773502691896257*fr[6])/dfac_v+(0.5773502691896257*fl[6])/dfac_v-(0.4472135954999579*fr[3])/dfac_v+(0.4472135954999579*fl[3])/dfac_v+(0.2581988897471612*fr[2])/dfac_v+(0.2581988897471612*fl[2])/dfac_v+0.8660254037844387*fr[7]*amax+0.8660254037844387*fl[7]*amax-0.5*fr[5]*amax+0.5*fl[5]*amax; 

  incr[0] = 0.5*Ghat[0]*dfac_x; 
  incr[1] = -0.8660254037844386*Ghat[0]*dfac_x; 
  incr[2] = 0.5*Ghat[2]*dfac_x; 
  incr[3] = -0.8660254037844386*Ghat[2]*dfac_x; 
  incr[4] = 1.118033988749895*Ghat[0]*dfac_x; 
  incr[5] = 0.5*Ghat[5]*dfac_x; 
  incr[6] = 1.118033988749895*Ghat[2]*dfac_x; 
  incr[7] = -0.8660254037844387*Ghat[5]*dfac_x; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 
  outr[4] += incr[4]; 
  outr[5] += incr[5]; 
  outr[6] += incr[6]; 
  outr[7] += incr[7]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += incr[1]; 
  outl[2] += -1.0*incr[2]; 
  outl[3] += incr[3]; 
  outl[4] += -1.0*incr[4]; 
  outl[5] += -1.0*incr[5]; 
  outl[6] += -1.0*incr[6]; 
  outl[7] += incr[7]; 
return std::abs(alpha0); 
} 
double EmGyrokineticSurfPositivity1x1vSer_Vpar_P1_Bvars_1(const double q_, const double m_, const double dt, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *Apar, const double *dApardt, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[4]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = -(0.5*(2.449489742783178*Phi[1]*dfac_x+1.414213562373095*dApardt[0])*q_)/m_; 

  double alpha[4]; 
  alpha[0] = (-(2.449489742783178*Phi[1]*dfac_x*q_)/m_)-(1.414213562373095*dApardt[0]*q_)/m_; 
  alpha[1] = -(1.414213562373095*dApardt[1]*q_)/m_; 
  if (alpha0>0) { 
  double rVal[4];  // rVal=f1/f0 at each node 
  rVal[0] = -(1.0*(1.732050807568877*fl[3]-3.0*fl[2]))/(3.464101615137754*EPSILON-1.0*fl[1]+1.732050807568877*fl[0]); 
  rVal[1] = (1.732050807568877*fl[3]+3.0*fl[2])/(3.464101615137754*EPSILON+fl[1]+1.732050807568877*fl[0]); 
  rVal[2] = -(1.0*(1.732050807568877*fl[3]-3.0*fl[2]))/(3.464101615137754*EPSILON-1.0*fl[1]+1.732050807568877*fl[0]); 
  rVal[3] = (1.732050807568877*fl[3]+3.0*fl[2])/(3.464101615137754*EPSILON+fl[1]+1.732050807568877*fl[0]); 
  double fqVal[4];  // fqVal = anti-limited f evaluated at each node 
  double cflVal = std::abs((alpha0*dt)/dxv[1]); 
  fqVal[0] = -0.2886751345948129*(fl[1]-1.732050807568877*fl[0])*limTheta(rVal[0],1.0,cflVal); 
  fqVal[1] = 0.2886751345948129*(fl[1]+1.732050807568877*fl[0])*limTheta(rVal[1],1.0,cflVal); 
  fqVal[2] = -0.2886751345948129*(fl[1]-1.732050807568877*fl[0])*limTheta(rVal[2],1.0,cflVal); 
  fqVal[3] = 0.2886751345948129*(fl[1]+1.732050807568877*fl[0])*limTheta(rVal[3],1.0,cflVal); 
  double fhatALVal[4];  // fhatALVal = mode coefficients of anti-limited f 
  fhatALVal[0] = 0.5*(fqVal[3]+fqVal[2]+fqVal[1]+fqVal[0]); 
  fhatALVal[1] = 0.8660254037844386*(fqVal[3]-1.0*fqVal[2]+fqVal[1]-1.0*fqVal[0]); 
  fhatALVal[2] = 0.8660254037844386*(fqVal[3]+fqVal[2]-1.0*(fqVal[1]+fqVal[0])); 
  fhatALVal[3] = 1.5*(fqVal[3]-1.0*(fqVal[2]+fqVal[1])+fqVal[0]); 
  incr[0] = 0.25*(1.732050807568877*alpha[1]*fhatALVal[3]+1.732050807568877*alpha[0]*fhatALVal[2]+alpha[1]*fhatALVal[1]+alpha[0]*fhatALVal[0])*dfac_v; 
  incr[1] = 0.25*(1.732050807568877*alpha[0]*fhatALVal[3]+1.732050807568877*alpha[1]*fhatALVal[2]+alpha[0]*fhatALVal[1]+fhatALVal[0]*alpha[1])*dfac_v; 
  incr[2] = -0.25*(3.0*alpha[1]*fhatALVal[3]+3.0*alpha[0]*fhatALVal[2]+1.732050807568877*alpha[1]*fhatALVal[1]+1.732050807568877*alpha[0]*fhatALVal[0])*dfac_v; 
  incr[3] = -0.25*(3.0*alpha[0]*fhatALVal[3]+3.0*alpha[1]*fhatALVal[2]+1.732050807568877*alpha[0]*fhatALVal[1]+1.732050807568877*fhatALVal[0]*alpha[1])*dfac_v; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += -1.0*incr[1]; 
  outl[2] += incr[2]; 
  outl[3] += incr[3]; 
  } else { 
  double rVal[4];  // rVal=f1/f0 at each node 
  rVal[0] = -(1.0*(1.732050807568877*fr[3]-3.0*fr[2]))/(3.464101615137754*EPSILON-1.0*fr[1]+1.732050807568877*fr[0]); 
  rVal[1] = (1.732050807568877*fr[3]+3.0*fr[2])/(3.464101615137754*EPSILON+fr[1]+1.732050807568877*fr[0]); 
  rVal[2] = -(1.0*(1.732050807568877*fr[3]-3.0*fr[2]))/(3.464101615137754*EPSILON-1.0*fr[1]+1.732050807568877*fr[0]); 
  rVal[3] = (1.732050807568877*fr[3]+3.0*fr[2])/(3.464101615137754*EPSILON+fr[1]+1.732050807568877*fr[0]); 
  double fqVal[4];  // fqVal = anti-limited f evaluated at each node 
  double cflVal = std::abs((alpha0*dt)/dxv[1]); 
  fqVal[0] = -0.2886751345948129*(fr[1]-1.732050807568877*fr[0])*limTheta(rVal[0],-1.0,cflVal); 
  fqVal[1] = 0.2886751345948129*(fr[1]+1.732050807568877*fr[0])*limTheta(rVal[1],-1.0,cflVal); 
  fqVal[2] = -0.2886751345948129*(fr[1]-1.732050807568877*fr[0])*limTheta(rVal[2],-1.0,cflVal); 
  fqVal[3] = 0.2886751345948129*(fr[1]+1.732050807568877*fr[0])*limTheta(rVal[3],-1.0,cflVal); 
  double fhatALVal[4];  // fhatALVal = mode coefficients of anti-limited f 
  fhatALVal[0] = 0.5*(fqVal[3]+fqVal[2]+fqVal[1]+fqVal[0]); 
  fhatALVal[1] = 0.8660254037844386*(fqVal[3]-1.0*fqVal[2]+fqVal[1]-1.0*fqVal[0]); 
  fhatALVal[2] = 0.8660254037844386*(fqVal[3]+fqVal[2]-1.0*(fqVal[1]+fqVal[0])); 
  fhatALVal[3] = 1.5*(fqVal[3]-1.0*(fqVal[2]+fqVal[1])+fqVal[0]); 
  incr[0] = -0.25*(1.732050807568877*alpha[1]*fhatALVal[3]+1.732050807568877*alpha[0]*fhatALVal[2]-1.0*alpha[1]*fhatALVal[1]-1.0*alpha[0]*fhatALVal[0])*dfac_v; 
  incr[1] = -0.25*(1.732050807568877*alpha[0]*fhatALVal[3]+1.732050807568877*alpha[1]*fhatALVal[2]-1.0*alpha[0]*fhatALVal[1]-1.0*fhatALVal[0]*alpha[1])*dfac_v; 
  incr[2] = 0.25*(3.0*alpha[1]*fhatALVal[3]+3.0*alpha[0]*fhatALVal[2]-1.732050807568877*alpha[1]*fhatALVal[1]-1.732050807568877*alpha[0]*fhatALVal[0])*dfac_v; 
  incr[3] = 0.25*(3.0*alpha[0]*fhatALVal[3]+3.0*alpha[1]*fhatALVal[2]-1.732050807568877*alpha[0]*fhatALVal[1]-1.732050807568877*fhatALVal[0]*alpha[1])*dfac_v; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += -1.0*incr[1]; 
  outl[2] += incr[2]; 
  outl[3] += incr[3]; 
  } 
  return std::abs(alpha0); 
} 
double EmGyrokineticSurfPositivity1x1vSer_Vpar_P2_Bvars_1(const double q_, const double m_, const double dt, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *Apar, const double *dApardt, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[8]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = -(0.5*(2.449489742783178*Phi[1]*dfac_x+1.414213562373095*dApardt[0])*q_)/m_; 

  double amax = 0.0; 
  bool upwind = false; 
  if(upwind) 
    amax = std::abs(alpha0); 
  else 
    amax = amax_in; 

  double Ghat[8]; 
  Ghat[0] = (-(3.061862178478972*Phi[2]*fr[7]*dfac_x*q_)/m_)-(3.061862178478972*Phi[2]*fl[7]*dfac_x*q_)/m_-(1.369306393762915*Phi[1]*fr[5]*dfac_x*q_)/m_-(1.369306393762915*Phi[1]*fl[5]*dfac_x*q_)/m_+(2.371708245126284*Phi[2]*fr[3]*dfac_x*q_)/m_-(2.371708245126284*Phi[2]*fl[3]*dfac_x*q_)/m_+(1.060660171779821*Phi[1]*fr[2]*dfac_x*q_)/m_-(1.060660171779821*Phi[1]*fl[2]*dfac_x*q_)/m_-(1.369306393762915*fr[1]*Phi[2]*dfac_x*q_)/m_-(1.369306393762915*fl[1]*Phi[2]*dfac_x*q_)/m_-(0.6123724356957944*fr[0]*Phi[1]*dfac_x*q_)/m_-(0.6123724356957944*fl[0]*Phi[1]*dfac_x*q_)/m_-(0.7905694150420948*dApardt[1]*fr[7]*q_)/m_-(0.7905694150420948*dApardt[1]*fl[7]*q_)/m_+(0.6123724356957944*dApardt[2]*fr[6]*q_)/m_-(0.6123724356957944*dApardt[2]*fl[6]*q_)/m_-(0.7905694150420947*dApardt[0]*fr[5]*q_)/m_-(0.7905694150420947*dApardt[0]*fl[5]*q_)/m_-(0.3535533905932737*dApardt[2]*fr[4]*q_)/m_-(0.3535533905932737*dApardt[2]*fl[4]*q_)/m_+(0.6123724356957944*dApardt[1]*fr[3]*q_)/m_-(0.6123724356957944*dApardt[1]*fl[3]*q_)/m_+(0.6123724356957944*dApardt[0]*fr[2]*q_)/m_-(0.6123724356957944*dApardt[0]*fl[2]*q_)/m_-(0.3535533905932737*dApardt[1]*fr[1]*q_)/m_-(0.3535533905932737*dApardt[1]*fl[1]*q_)/m_-(0.3535533905932737*dApardt[0]*fr[0]*q_)/m_-(0.3535533905932737*dApardt[0]*fl[0]*q_)/m_-1.118033988749895*fr[5]*amax+1.118033988749895*fl[5]*amax+0.8660254037844386*fr[2]*amax+0.8660254037844386*fl[2]*amax-0.5*fr[0]*amax+0.5*fl[0]*amax; 
  Ghat[1] = (-(1.369306393762915*Phi[1]*fr[7]*dfac_x*q_)/m_)-(1.369306393762915*Phi[1]*fl[7]*dfac_x*q_)/m_+(2.121320343559642*Phi[2]*fr[6]*dfac_x*q_)/m_-(2.121320343559642*Phi[2]*fl[6]*dfac_x*q_)/m_-(3.061862178478972*Phi[2]*fr[5]*dfac_x*q_)/m_-(3.061862178478972*Phi[2]*fl[5]*dfac_x*q_)/m_-(1.224744871391589*Phi[2]*fr[4]*dfac_x*q_)/m_-(1.224744871391589*Phi[2]*fl[4]*dfac_x*q_)/m_+(1.060660171779821*Phi[1]*fr[3]*dfac_x*q_)/m_-(1.060660171779821*Phi[1]*fl[3]*dfac_x*q_)/m_+(2.371708245126284*Phi[2]*fr[2]*dfac_x*q_)/m_-(2.371708245126284*Phi[2]*fl[2]*dfac_x*q_)/m_-(1.369306393762915*fr[0]*Phi[2]*dfac_x*q_)/m_-(1.369306393762915*fl[0]*Phi[2]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fr[1]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fl[1]*dfac_x*q_)/m_-(0.7071067811865475*dApardt[2]*fr[7]*q_)/m_-(0.7905694150420948*dApardt[0]*fr[7]*q_)/m_-(0.7071067811865475*dApardt[2]*fl[7]*q_)/m_-(0.7905694150420948*dApardt[0]*fl[7]*q_)/m_+(0.5477225575051661*dApardt[1]*fr[6]*q_)/m_-(0.5477225575051661*dApardt[1]*fl[6]*q_)/m_-(0.7905694150420947*dApardt[1]*fr[5]*q_)/m_-(0.7905694150420947*dApardt[1]*fl[5]*q_)/m_-(0.3162277660168379*dApardt[1]*fr[4]*q_)/m_-(0.3162277660168379*dApardt[1]*fl[4]*q_)/m_+(0.5477225575051661*dApardt[2]*fr[3]*q_)/m_+(0.6123724356957944*dApardt[0]*fr[3]*q_)/m_-(0.5477225575051661*dApardt[2]*fl[3]*q_)/m_-(0.6123724356957944*dApardt[0]*fl[3]*q_)/m_+(0.6123724356957944*dApardt[1]*fr[2]*q_)/m_-(0.6123724356957944*dApardt[1]*fl[2]*q_)/m_-(0.3162277660168379*fr[1]*dApardt[2]*q_)/m_-(0.3162277660168379*fl[1]*dApardt[2]*q_)/m_-(0.3535533905932737*dApardt[0]*fr[1]*q_)/m_-(0.3535533905932737*dApardt[0]*fl[1]*q_)/m_-(0.3535533905932737*fr[0]*dApardt[1]*q_)/m_-(0.3535533905932737*fl[0]*dApardt[1]*q_)/m_-1.118033988749895*fr[7]*amax+1.118033988749895*fl[7]*amax+0.8660254037844386*fr[3]*amax+0.8660254037844386*fl[3]*amax-0.5*fr[1]*amax+0.5*fl[1]*amax; 
  Ghat[4] = (-(2.738612787525831*Phi[2]*fr[7]*dfac_x*q_)/m_)-(2.738612787525831*Phi[2]*fl[7]*dfac_x*q_)/m_+(1.060660171779821*Phi[1]*fr[6]*dfac_x*q_)/m_-(1.060660171779821*Phi[1]*fl[6]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fr[4]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fl[4]*dfac_x*q_)/m_+(2.121320343559642*Phi[2]*fr[3]*dfac_x*q_)/m_-(2.121320343559642*Phi[2]*fl[3]*dfac_x*q_)/m_-(1.224744871391589*fr[1]*Phi[2]*dfac_x*q_)/m_-(1.224744871391589*fl[1]*Phi[2]*dfac_x*q_)/m_-(0.7071067811865475*dApardt[1]*fr[7]*q_)/m_-(0.7071067811865475*dApardt[1]*fl[7]*q_)/m_+(0.3912303982179757*dApardt[2]*fr[6]*q_)/m_+(0.6123724356957944*dApardt[0]*fr[6]*q_)/m_-(0.3912303982179757*dApardt[2]*fl[6]*q_)/m_-(0.6123724356957944*dApardt[0]*fl[6]*q_)/m_-(0.7905694150420947*dApardt[2]*fr[5]*q_)/m_-(0.7905694150420947*dApardt[2]*fl[5]*q_)/m_-(0.2258769757263128*dApardt[2]*fr[4]*q_)/m_-(0.3535533905932737*dApardt[0]*fr[4]*q_)/m_-(0.2258769757263128*dApardt[2]*fl[4]*q_)/m_-(0.3535533905932737*dApardt[0]*fl[4]*q_)/m_+(0.5477225575051661*dApardt[1]*fr[3]*q_)/m_-(0.5477225575051661*dApardt[1]*fl[3]*q_)/m_+(0.6123724356957944*dApardt[2]*fr[2]*q_)/m_-(0.6123724356957944*dApardt[2]*fl[2]*q_)/m_-(0.3535533905932737*fr[0]*dApardt[2]*q_)/m_-(0.3535533905932737*fl[0]*dApardt[2]*q_)/m_-(0.3162277660168379*dApardt[1]*fr[1]*q_)/m_-(0.3162277660168379*dApardt[1]*fl[1]*q_)/m_+0.8660254037844387*fr[6]*amax+0.8660254037844387*fl[6]*amax-0.5*fr[4]*amax+0.5*fl[4]*amax; 

  incr[0] = 0.5*Ghat[0]*dfac_v; 
  incr[1] = 0.5*Ghat[1]*dfac_v; 
  incr[2] = -0.8660254037844386*Ghat[0]*dfac_v; 
  incr[3] = -0.8660254037844386*Ghat[1]*dfac_v; 
  incr[4] = 0.5*Ghat[4]*dfac_v; 
  incr[5] = 1.118033988749895*Ghat[0]*dfac_v; 
  incr[6] = -0.8660254037844387*Ghat[4]*dfac_v; 
  incr[7] = 1.118033988749895*Ghat[1]*dfac_v; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 
  outr[4] += incr[4]; 
  outr[5] += incr[5]; 
  outr[6] += incr[6]; 
  outr[7] += incr[7]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += -1.0*incr[1]; 
  outl[2] += incr[2]; 
  outl[3] += incr[3]; 
  outl[4] += -1.0*incr[4]; 
  outl[5] += -1.0*incr[5]; 
  outl[6] += incr[6]; 
  outl[7] += -1.0*incr[7]; 
return std::abs(alpha0); 
} 
