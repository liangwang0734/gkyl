#include <VmLBOModDecl.h> 
double VmLBOconstNuVol1x3vSerP1(const double *w, const double *dxv, const double nu, const double *u, const double *vtSq, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates. 
  // dxv[4]: Cell spacing. 
  // nu:     diffusion coefficient (collisionality). 
  // u:    bulk velocity. 
  // vtSq: thermal speed squared. 
  // f:      Input distribution function.
  // out:    Incremented output 
  const double rdvx2nu = 2.0*nu/dxv[1]; 
  const double rdvxSq4nu = 4.0*nu/(dxv[1]*dxv[1]); 
  const double rdvy2nu = 2.0*nu/dxv[2]; 
  const double rdvySq4nu = 4.0*nu/(dxv[2]*dxv[2]); 
  const double rdvz2nu = 2.0*nu/dxv[3]; 
  const double rdvzSq4nu = 4.0*nu/(dxv[3]*dxv[3]); 

  double alpha_mid = 0.0; 
  double alpha_drag[48]; 
  double alpha_diffusion[16]; 

  // Expand rdv2nu*(vx-ux) in phase basis.
  alpha_drag[0] = (2.828427124746191*u[0]-4.0*w[1])*rdvx2nu; 
  alpha_drag[1] = 2.828427124746191*u[1]*rdvx2nu; 
  alpha_drag[2] = -1.154700538379252*dxv[1]*rdvx2nu; 

  // x-drag contribution to the midpoint value used for CFL.
  alpha_mid += std::abs(0.125*alpha_drag[0]); 

  // Expand rdv2nu*(vy-uy) in phase basis.
  alpha_drag[16] = (2.828427124746191*u[2]-4.0*w[2])*rdvy2nu; 
  alpha_drag[17] = 2.828427124746191*u[3]*rdvy2nu; 
  alpha_drag[19] = -1.154700538379252*dxv[2]*rdvy2nu; 

  // y-drag contribution to the midpoint value used for CFL.
  alpha_mid += std::abs(0.125*alpha_drag[16]); 

  // Expand rdv2nu*(vz-uz) in phase basis.
  alpha_drag[32] = (2.828427124746191*u[4]-4.0*w[3])*rdvz2nu; 
  alpha_drag[33] = 2.828427124746191*u[5]*rdvz2nu; 
  alpha_drag[36] = -1.154700538379252*dxv[3]*rdvz2nu; 

  // z-drag contribution to the midpoint value used for CFL.
  alpha_mid += std::abs(0.125*alpha_drag[32]); 

  // Expand nu*vthSq in phase basis.
  alpha_diffusion[0] = 2.828427124746191*vtSq[0]; 
  alpha_diffusion[1] = 2.828427124746191*vtSq[1]; 

  // Diffusion contribution to the midpoint value used for CFL.
  alpha_mid += std::abs(0.3333333333333333*alpha_diffusion[0]*rdvxSq4nu); 
  alpha_mid += std::abs(0.3333333333333333*alpha_diffusion[0]*rdvySq4nu); 
  alpha_mid += std::abs(0.3333333333333333*alpha_diffusion[0]*rdvzSq4nu); 

  // Put together updates due to drag and diffusion terms.
  out[2] += 0.4330127018922193*(alpha_drag[2]*f[2]+alpha_drag[1]*f[1]+alpha_drag[0]*f[0]); 
  out[3] += 0.4330127018922193*(f[3]*alpha_drag[19]+f[1]*alpha_drag[17]+f[0]*alpha_drag[16]); 
  out[4] += 0.4330127018922193*(f[4]*alpha_drag[36]+f[1]*alpha_drag[33]+f[0]*alpha_drag[32]); 
  out[5] += 0.4330127018922193*(alpha_drag[2]*f[5]+alpha_drag[0]*f[1]+f[0]*alpha_drag[1]); 
  out[6] += 0.4330127018922193*(f[6]*alpha_drag[19]+f[0]*alpha_drag[17]+f[1]*alpha_drag[16]); 
  out[7] += 0.4330127018922193*(f[7]*alpha_drag[19]+f[5]*alpha_drag[17]+f[2]*alpha_drag[16]+alpha_drag[2]*f[7]+alpha_drag[1]*f[6]+alpha_drag[0]*f[3]); 
  out[8] += 0.4330127018922193*(f[8]*alpha_drag[36]+f[0]*alpha_drag[33]+f[1]*alpha_drag[32]); 
  out[9] += 0.4330127018922193*(f[9]*alpha_drag[36]+f[5]*alpha_drag[33]+f[2]*alpha_drag[32]+alpha_drag[2]*f[9]+alpha_drag[1]*f[8]+alpha_drag[0]*f[4]); 
  out[10] += 0.4330127018922193*(f[10]*alpha_drag[36]+f[6]*alpha_drag[33]+f[3]*alpha_drag[32]+f[10]*alpha_drag[19]+f[8]*alpha_drag[17]+f[4]*alpha_drag[16]); 
  out[11] += 0.4330127018922193*(f[11]*alpha_drag[19]+f[2]*alpha_drag[17]+f[5]*alpha_drag[16]+alpha_drag[2]*f[11]+alpha_drag[0]*f[6]+alpha_drag[1]*f[3]); 
  out[12] += 0.4330127018922193*(f[12]*alpha_drag[36]+f[2]*alpha_drag[33]+f[5]*alpha_drag[32]+alpha_drag[2]*f[12]+alpha_drag[0]*f[8]+alpha_drag[1]*f[4]); 
  out[13] += 0.4330127018922193*(f[13]*alpha_drag[36]+f[3]*alpha_drag[33]+f[6]*alpha_drag[32]+f[13]*alpha_drag[19]+f[4]*alpha_drag[17]+f[8]*alpha_drag[16]); 
  out[14] += 0.4330127018922193*(f[14]*alpha_drag[36]+f[11]*alpha_drag[33]+f[7]*alpha_drag[32]+f[14]*alpha_drag[19]+f[12]*alpha_drag[17]+f[9]*alpha_drag[16]+alpha_drag[2]*f[14]+alpha_drag[1]*f[13]+alpha_drag[0]*f[10]); 
  out[15] += 0.4330127018922193*(f[15]*alpha_drag[36]+f[7]*alpha_drag[33]+f[11]*alpha_drag[32]+f[15]*alpha_drag[19]+f[9]*alpha_drag[17]+f[12]*alpha_drag[16]+alpha_drag[2]*f[15]+alpha_drag[0]*f[13]+alpha_drag[1]*f[10]); 

  return alpha_mid; 

} 
double VmLBOconstNuVol1x3vSerP2(const double *w, const double *dxv, const double nu, const double *u, const double *vtSq, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates. 
  // dxv[4]: Cell spacing. 
  // nu:     diffusion coefficient (collisionality). 
  // u:    bulk velocity. 
  // vtSq: thermal speed squared. 
  // f:      Input distribution function.
  // out:    Incremented output 
  const double rdvx2nu = 2.0*nu/dxv[1]; 
  const double rdvxSq4nu = 4.0*nu/(dxv[1]*dxv[1]); 
  const double rdvy2nu = 2.0*nu/dxv[2]; 
  const double rdvySq4nu = 4.0*nu/(dxv[2]*dxv[2]); 
  const double rdvz2nu = 2.0*nu/dxv[3]; 
  const double rdvzSq4nu = 4.0*nu/(dxv[3]*dxv[3]); 

  double alpha_mid = 0.0; 
  double alpha_drag[144]; 
  double alpha_diffusion[48]; 

  // Expand rdv2nu*(vx-ux) in phase basis.
  alpha_drag[0] = (2.828427124746191*u[0]-4.0*w[1])*rdvx2nu; 
  alpha_drag[1] = 2.828427124746191*u[1]*rdvx2nu; 
  alpha_drag[2] = -1.154700538379252*dxv[1]*rdvx2nu; 
  alpha_drag[11] = 2.828427124746191*u[2]*rdvx2nu; 

  // x-drag contribution to the midpoint value used for CFL.
  alpha_mid += std::abs(0.125*alpha_drag[0]-0.1397542485937369*alpha_drag[11]); 

  // Expand rdv2nu*(vy-uy) in phase basis.
  alpha_drag[48] = (2.828427124746191*u[3]-4.0*w[2])*rdvy2nu; 
  alpha_drag[49] = 2.828427124746191*u[4]*rdvy2nu; 
  alpha_drag[51] = -1.154700538379252*dxv[2]*rdvy2nu; 
  alpha_drag[59] = 2.828427124746191*u[5]*rdvy2nu; 

  // y-drag contribution to the midpoint value used for CFL.
  alpha_mid += std::abs(0.125*alpha_drag[48]-0.1397542485937369*alpha_drag[59]); 

  // Expand rdv2nu*(vz-uz) in phase basis.
  alpha_drag[96] = (2.828427124746191*u[6]-4.0*w[3])*rdvz2nu; 
  alpha_drag[97] = 2.828427124746191*u[7]*rdvz2nu; 
  alpha_drag[100] = -1.154700538379252*dxv[3]*rdvz2nu; 
  alpha_drag[107] = 2.828427124746191*u[8]*rdvz2nu; 

  // z-drag contribution to the midpoint value used for CFL.
  alpha_mid += std::abs(0.125*alpha_drag[96]-0.1397542485937369*alpha_drag[107]); 

  // Expand nu*vthSq in phase basis.
  alpha_diffusion[0] = 2.828427124746191*vtSq[0]; 
  alpha_diffusion[1] = 2.828427124746191*vtSq[1]; 
  alpha_diffusion[11] = 2.828427124746191*vtSq[2]; 

  // Diffusion contribution to the midpoint value used for CFL.
  alpha_mid += std::abs((0.45*alpha_diffusion[0]-0.5031152949374527*alpha_diffusion[11])*rdvxSq4nu); 
  alpha_mid += std::abs((0.45*alpha_diffusion[0]-0.5031152949374527*alpha_diffusion[11])*rdvySq4nu); 
  alpha_mid += std::abs((0.45*alpha_diffusion[0]-0.5031152949374527*alpha_diffusion[11])*rdvzSq4nu); 

  // Put together updates due to drag and diffusion terms.
  out[2] += 0.4330127018922193*(alpha_drag[11]*f[11]+alpha_drag[2]*f[2]+alpha_drag[1]*f[1]+alpha_drag[0]*f[0]); 
  out[3] += 0.4330127018922193*(f[11]*alpha_drag[59]+f[3]*alpha_drag[51]+f[1]*alpha_drag[49]+f[0]*alpha_drag[48]); 
  out[4] += 0.4330127018922193*(f[11]*alpha_drag[107]+f[4]*alpha_drag[100]+f[1]*alpha_drag[97]+f[0]*alpha_drag[96]); 
  out[5] += 0.3872983346207416*(alpha_drag[1]*f[11]+f[1]*alpha_drag[11])+0.4330127018922193*(alpha_drag[2]*f[5]+alpha_drag[0]*f[1]+f[0]*alpha_drag[1]); 
  out[6] += 0.3872983346207416*f[1]*alpha_drag[59]+0.4330127018922193*f[6]*alpha_drag[51]+0.3872983346207416*f[11]*alpha_drag[49]+0.4330127018922193*(f[0]*alpha_drag[49]+f[1]*alpha_drag[48]); 
  out[7] += 0.4330127018922193*(f[19]*alpha_drag[59]+f[7]*alpha_drag[51]+f[5]*alpha_drag[49]+f[2]*alpha_drag[48]+alpha_drag[11]*f[21]+alpha_drag[2]*f[7]+alpha_drag[1]*f[6]+alpha_drag[0]*f[3]); 
  out[8] += 0.3872983346207416*f[1]*alpha_drag[107]+0.4330127018922193*f[8]*alpha_drag[100]+0.3872983346207416*f[11]*alpha_drag[97]+0.4330127018922193*(f[0]*alpha_drag[97]+f[1]*alpha_drag[96]); 
  out[9] += 0.4330127018922193*(f[19]*alpha_drag[107]+f[9]*alpha_drag[100]+f[5]*alpha_drag[97]+f[2]*alpha_drag[96]+alpha_drag[11]*f[25]+alpha_drag[2]*f[9]+alpha_drag[1]*f[8]+alpha_drag[0]*f[4]); 
  out[10] += 0.4330127018922193*(f[21]*alpha_drag[107]+f[10]*alpha_drag[100]+f[6]*alpha_drag[97]+f[3]*alpha_drag[96]+f[25]*alpha_drag[59]+f[10]*alpha_drag[51]+f[8]*alpha_drag[49]+f[4]*alpha_drag[48]); 
  out[12] += 1.677050983124842*(alpha_diffusion[11]*f[11]+alpha_diffusion[1]*f[1]+alpha_diffusion[0]*f[0])*rdvxSq4nu+0.9682458365518543*alpha_drag[11]*f[19]+0.8660254037844386*alpha_drag[2]*f[12]+0.9682458365518543*(alpha_drag[1]*f[5]+alpha_drag[0]*f[2]+f[0]*alpha_drag[2]); 
  out[13] += 1.677050983124842*(alpha_diffusion[11]*f[11]+alpha_diffusion[1]*f[1]+alpha_diffusion[0]*f[0])*rdvySq4nu+0.9682458365518543*f[21]*alpha_drag[59]+0.8660254037844386*f[13]*alpha_drag[51]+0.9682458365518543*(f[0]*alpha_drag[51]+f[6]*alpha_drag[49]+f[3]*alpha_drag[48]); 
  out[14] += 1.677050983124842*(alpha_diffusion[11]*f[11]+alpha_diffusion[1]*f[1]+alpha_diffusion[0]*f[0])*rdvzSq4nu+0.9682458365518543*f[25]*alpha_drag[107]+0.8660254037844386*f[14]*alpha_drag[100]+0.9682458365518543*(f[0]*alpha_drag[100]+f[8]*alpha_drag[97]+f[4]*alpha_drag[96]); 
  out[15] += 0.3872983346207416*f[5]*alpha_drag[59]+0.4330127018922193*f[15]*alpha_drag[51]+0.3872983346207416*f[19]*alpha_drag[49]+0.4330127018922193*(f[2]*alpha_drag[49]+f[5]*alpha_drag[48])+0.3872983346207416*alpha_drag[1]*f[21]+0.4330127018922193*alpha_drag[2]*f[15]+0.3872983346207416*f[6]*alpha_drag[11]+0.4330127018922193*(alpha_drag[0]*f[6]+alpha_drag[1]*f[3]); 
  out[16] += 0.3872983346207416*f[5]*alpha_drag[107]+0.4330127018922193*f[16]*alpha_drag[100]+0.3872983346207416*f[19]*alpha_drag[97]+0.4330127018922193*(f[2]*alpha_drag[97]+f[5]*alpha_drag[96])+0.3872983346207416*alpha_drag[1]*f[25]+0.4330127018922193*alpha_drag[2]*f[16]+0.3872983346207416*f[8]*alpha_drag[11]+0.4330127018922193*(alpha_drag[0]*f[8]+alpha_drag[1]*f[4]); 
  out[17] += 0.3872983346207416*f[6]*alpha_drag[107]+0.4330127018922193*f[17]*alpha_drag[100]+0.3872983346207416*f[21]*alpha_drag[97]+0.4330127018922193*(f[3]*alpha_drag[97]+f[6]*alpha_drag[96])+0.3872983346207416*f[8]*alpha_drag[59]+0.4330127018922193*f[17]*alpha_drag[51]+0.3872983346207416*f[25]*alpha_drag[49]+0.4330127018922193*(f[4]*alpha_drag[49]+f[8]*alpha_drag[48]); 
  out[18] += 0.4330127018922193*(f[32]*alpha_drag[107]+f[18]*alpha_drag[100]+f[15]*alpha_drag[97]+f[7]*alpha_drag[96]+f[35]*alpha_drag[59]+f[18]*alpha_drag[51]+f[16]*alpha_drag[49]+f[9]*alpha_drag[48]+alpha_drag[11]*f[37]+alpha_drag[2]*f[18]+alpha_drag[1]*f[17]+alpha_drag[0]*f[10]); 
  out[19] += 0.4330127018922193*alpha_drag[2]*f[19]+0.276641667586244*alpha_drag[11]*f[11]+0.4330127018922193*(alpha_drag[0]*f[11]+f[0]*alpha_drag[11])+0.3872983346207416*alpha_drag[1]*f[1]; 
  out[20] += (1.5*(alpha_diffusion[1]*f[11]+f[1]*alpha_diffusion[11])+1.677050983124842*(alpha_diffusion[0]*f[1]+f[0]*alpha_diffusion[1]))*rdvxSq4nu+0.8660254037844386*(alpha_drag[2]*f[20]+alpha_drag[1]*f[19]+f[5]*alpha_drag[11])+0.9682458365518543*(alpha_drag[0]*f[5]+alpha_drag[1]*f[2]+f[1]*alpha_drag[2]); 
  out[21] += 0.276641667586244*f[11]*alpha_drag[59]+0.4330127018922193*(f[0]*alpha_drag[59]+f[21]*alpha_drag[51])+0.3872983346207416*f[1]*alpha_drag[49]+0.4330127018922193*f[11]*alpha_drag[48]; 
  out[22] += 1.677050983124842*(alpha_diffusion[11]*f[21]+alpha_diffusion[1]*f[6]+alpha_diffusion[0]*f[3])*rdvxSq4nu+0.4330127018922193*(f[22]*alpha_drag[51]+f[20]*alpha_drag[49]+f[12]*alpha_drag[48])+0.9682458365518543*alpha_drag[11]*f[32]+0.8660254037844386*alpha_drag[2]*f[22]+0.9682458365518543*(alpha_drag[1]*f[15]+alpha_drag[0]*f[7]+alpha_drag[2]*f[3]); 
  out[23] += (1.5*(alpha_diffusion[1]*f[11]+f[1]*alpha_diffusion[11])+1.677050983124842*(alpha_diffusion[0]*f[1]+f[0]*alpha_diffusion[1]))*rdvySq4nu+0.8660254037844386*f[6]*alpha_drag[59]+(0.8660254037844386*f[23]+0.9682458365518543*f[1])*alpha_drag[51]+0.8660254037844386*f[21]*alpha_drag[49]+0.9682458365518543*(f[3]*alpha_drag[49]+f[6]*alpha_drag[48]); 
  out[24] += 1.677050983124842*(alpha_diffusion[11]*f[19]+alpha_diffusion[1]*f[5]+alpha_diffusion[0]*f[2])*rdvySq4nu+0.9682458365518543*f[32]*alpha_drag[59]+0.8660254037844386*f[24]*alpha_drag[51]+0.9682458365518543*(f[2]*alpha_drag[51]+f[15]*alpha_drag[49]+f[7]*alpha_drag[48])+0.4330127018922193*(alpha_drag[2]*f[24]+alpha_drag[1]*f[23]+alpha_drag[0]*f[13]); 
  out[25] += 0.276641667586244*f[11]*alpha_drag[107]+0.4330127018922193*(f[0]*alpha_drag[107]+f[25]*alpha_drag[100])+0.3872983346207416*f[1]*alpha_drag[97]+0.4330127018922193*f[11]*alpha_drag[96]; 
  out[26] += 1.677050983124842*(alpha_diffusion[11]*f[25]+alpha_diffusion[1]*f[8]+alpha_diffusion[0]*f[4])*rdvxSq4nu+0.4330127018922193*(f[26]*alpha_drag[100]+f[20]*alpha_drag[97]+f[12]*alpha_drag[96])+0.9682458365518543*alpha_drag[11]*f[35]+0.8660254037844386*alpha_drag[2]*f[26]+0.9682458365518543*(alpha_drag[1]*f[16]+alpha_drag[0]*f[9]+alpha_drag[2]*f[4]); 
  out[27] += 1.677050983124842*(alpha_diffusion[11]*f[25]+alpha_diffusion[1]*f[8]+alpha_diffusion[0]*f[4])*rdvySq4nu+0.4330127018922193*(f[27]*alpha_drag[100]+f[23]*alpha_drag[97]+f[13]*alpha_drag[96])+0.9682458365518543*f[37]*alpha_drag[59]+0.8660254037844386*f[27]*alpha_drag[51]+0.9682458365518543*(f[4]*alpha_drag[51]+f[17]*alpha_drag[49]+f[10]*alpha_drag[48]); 
  out[28] += (1.5*(alpha_diffusion[1]*f[11]+f[1]*alpha_diffusion[11])+1.677050983124842*(alpha_diffusion[0]*f[1]+f[0]*alpha_diffusion[1]))*rdvzSq4nu+0.8660254037844386*f[8]*alpha_drag[107]+(0.8660254037844386*f[28]+0.9682458365518543*f[1])*alpha_drag[100]+0.8660254037844386*f[25]*alpha_drag[97]+0.9682458365518543*(f[4]*alpha_drag[97]+f[8]*alpha_drag[96]); 
  out[29] += 1.677050983124842*(alpha_diffusion[11]*f[19]+alpha_diffusion[1]*f[5]+alpha_diffusion[0]*f[2])*rdvzSq4nu+0.9682458365518543*f[35]*alpha_drag[107]+0.8660254037844386*f[29]*alpha_drag[100]+0.9682458365518543*(f[2]*alpha_drag[100]+f[16]*alpha_drag[97]+f[9]*alpha_drag[96])+0.4330127018922193*(alpha_drag[2]*f[29]+alpha_drag[1]*f[28]+alpha_drag[0]*f[14]); 
  out[30] += 1.677050983124842*(alpha_diffusion[11]*f[21]+alpha_diffusion[1]*f[6]+alpha_diffusion[0]*f[3])*rdvzSq4nu+0.9682458365518543*f[37]*alpha_drag[107]+0.8660254037844386*f[30]*alpha_drag[100]+0.9682458365518543*(f[3]*alpha_drag[100]+f[17]*alpha_drag[97]+f[10]*alpha_drag[96])+0.4330127018922193*(f[30]*alpha_drag[51]+f[28]*alpha_drag[49]+f[14]*alpha_drag[48]); 
  out[31] += 0.3872983346207416*f[15]*alpha_drag[107]+0.4330127018922193*f[31]*alpha_drag[100]+0.3872983346207416*f[32]*alpha_drag[97]+0.4330127018922193*(f[7]*alpha_drag[97]+f[15]*alpha_drag[96])+0.3872983346207416*f[16]*alpha_drag[59]+0.4330127018922193*f[31]*alpha_drag[51]+0.3872983346207416*f[35]*alpha_drag[49]+0.4330127018922193*(f[9]*alpha_drag[49]+f[16]*alpha_drag[48])+0.3872983346207416*alpha_drag[1]*f[37]+0.4330127018922193*alpha_drag[2]*f[31]+0.3872983346207416*alpha_drag[11]*f[17]+0.4330127018922193*(alpha_drag[0]*f[17]+alpha_drag[1]*f[10]); 
  out[32] += 0.276641667586244*f[19]*alpha_drag[59]+0.4330127018922193*(f[2]*alpha_drag[59]+f[32]*alpha_drag[51])+0.3872983346207416*f[5]*alpha_drag[49]+0.4330127018922193*(f[19]*alpha_drag[48]+alpha_drag[2]*f[32])+0.276641667586244*alpha_drag[11]*f[21]+0.4330127018922193*(alpha_drag[0]*f[21]+f[3]*alpha_drag[11])+0.3872983346207416*alpha_drag[1]*f[6]; 
  out[33] += (1.5*(alpha_diffusion[1]*f[21]+f[6]*alpha_diffusion[11])+1.677050983124842*(alpha_diffusion[0]*f[6]+alpha_diffusion[1]*f[3]))*rdvxSq4nu+0.3872983346207416*f[20]*alpha_drag[59]+0.4330127018922193*(f[33]*alpha_drag[51]+f[12]*alpha_drag[49]+f[20]*alpha_drag[48])+0.8660254037844386*(alpha_drag[2]*f[33]+alpha_drag[1]*f[32]+alpha_drag[11]*f[15])+0.9682458365518543*(alpha_drag[0]*f[15]+alpha_drag[1]*f[7]+alpha_drag[2]*f[6]); 
  out[34] += (1.5*(alpha_diffusion[1]*f[19]+f[5]*alpha_diffusion[11])+1.677050983124842*(alpha_diffusion[0]*f[5]+alpha_diffusion[1]*f[2]))*rdvySq4nu+0.8660254037844386*f[15]*alpha_drag[59]+(0.8660254037844386*f[34]+0.9682458365518543*f[5])*alpha_drag[51]+0.8660254037844386*f[32]*alpha_drag[49]+0.9682458365518543*(f[7]*alpha_drag[49]+f[15]*alpha_drag[48])+0.4330127018922193*alpha_drag[2]*f[34]+0.3872983346207416*alpha_drag[11]*f[23]+0.4330127018922193*(alpha_drag[0]*f[23]+alpha_drag[1]*f[13]); 
  out[35] += 0.276641667586244*f[19]*alpha_drag[107]+0.4330127018922193*(f[2]*alpha_drag[107]+f[35]*alpha_drag[100])+0.3872983346207416*f[5]*alpha_drag[97]+0.4330127018922193*(f[19]*alpha_drag[96]+alpha_drag[2]*f[35])+0.276641667586244*alpha_drag[11]*f[25]+0.4330127018922193*(alpha_drag[0]*f[25]+f[4]*alpha_drag[11])+0.3872983346207416*alpha_drag[1]*f[8]; 
  out[36] += (1.5*(alpha_diffusion[1]*f[25]+f[8]*alpha_diffusion[11])+1.677050983124842*(alpha_diffusion[0]*f[8]+alpha_diffusion[1]*f[4]))*rdvxSq4nu+0.3872983346207416*f[20]*alpha_drag[107]+0.4330127018922193*(f[36]*alpha_drag[100]+f[12]*alpha_drag[97]+f[20]*alpha_drag[96])+0.8660254037844386*(alpha_drag[2]*f[36]+alpha_drag[1]*f[35]+alpha_drag[11]*f[16])+0.9682458365518543*(alpha_drag[0]*f[16]+alpha_drag[1]*f[9]+alpha_drag[2]*f[8]); 
  out[37] += 0.276641667586244*f[21]*alpha_drag[107]+0.4330127018922193*(f[3]*alpha_drag[107]+f[37]*alpha_drag[100])+0.3872983346207416*f[6]*alpha_drag[97]+0.4330127018922193*f[21]*alpha_drag[96]+0.276641667586244*f[25]*alpha_drag[59]+0.4330127018922193*(f[4]*alpha_drag[59]+f[37]*alpha_drag[51])+0.3872983346207416*f[8]*alpha_drag[49]+0.4330127018922193*f[25]*alpha_drag[48]; 
  out[38] += 1.677050983124842*(alpha_diffusion[11]*f[37]+alpha_diffusion[1]*f[17]+alpha_diffusion[0]*f[10])*rdvxSq4nu+0.4330127018922193*(f[38]*alpha_drag[100]+f[33]*alpha_drag[97]+f[22]*alpha_drag[96]+f[38]*alpha_drag[51]+f[36]*alpha_drag[49]+f[26]*alpha_drag[48])+0.9682458365518543*alpha_drag[11]*f[44]+0.8660254037844386*alpha_drag[2]*f[38]+0.9682458365518543*(alpha_drag[1]*f[31]+alpha_drag[0]*f[18]+alpha_drag[2]*f[10]); 
  out[39] += (1.5*(alpha_diffusion[1]*f[25]+f[8]*alpha_diffusion[11])+1.677050983124842*(alpha_diffusion[0]*f[8]+alpha_diffusion[1]*f[4]))*rdvySq4nu+0.3872983346207416*f[23]*alpha_drag[107]+0.4330127018922193*(f[39]*alpha_drag[100]+f[13]*alpha_drag[97]+f[23]*alpha_drag[96])+0.8660254037844386*f[17]*alpha_drag[59]+(0.8660254037844386*f[39]+0.9682458365518543*f[8])*alpha_drag[51]+0.8660254037844386*f[37]*alpha_drag[49]+0.9682458365518543*(f[10]*alpha_drag[49]+f[17]*alpha_drag[48]); 
  out[40] += 1.677050983124842*(alpha_diffusion[11]*f[35]+alpha_diffusion[1]*f[16]+alpha_diffusion[0]*f[9])*rdvySq4nu+0.4330127018922193*(f[40]*alpha_drag[100]+f[34]*alpha_drag[97]+f[24]*alpha_drag[96])+0.9682458365518543*f[44]*alpha_drag[59]+0.8660254037844386*f[40]*alpha_drag[51]+0.9682458365518543*(f[9]*alpha_drag[51]+f[31]*alpha_drag[49]+f[18]*alpha_drag[48])+0.4330127018922193*(alpha_drag[2]*f[40]+alpha_drag[1]*f[39]+alpha_drag[0]*f[27]); 
  out[41] += (1.5*(alpha_diffusion[1]*f[19]+f[5]*alpha_diffusion[11])+1.677050983124842*(alpha_diffusion[0]*f[5]+alpha_diffusion[1]*f[2]))*rdvzSq4nu+0.8660254037844386*f[16]*alpha_drag[107]+(0.8660254037844386*f[41]+0.9682458365518543*f[5])*alpha_drag[100]+0.8660254037844386*f[35]*alpha_drag[97]+0.9682458365518543*(f[9]*alpha_drag[97]+f[16]*alpha_drag[96])+0.4330127018922193*alpha_drag[2]*f[41]+0.3872983346207416*alpha_drag[11]*f[28]+0.4330127018922193*(alpha_drag[0]*f[28]+alpha_drag[1]*f[14]); 
  out[42] += (1.5*(alpha_diffusion[1]*f[21]+f[6]*alpha_diffusion[11])+1.677050983124842*(alpha_diffusion[0]*f[6]+alpha_diffusion[1]*f[3]))*rdvzSq4nu+0.8660254037844386*f[17]*alpha_drag[107]+(0.8660254037844386*f[42]+0.9682458365518543*f[6])*alpha_drag[100]+0.8660254037844386*f[37]*alpha_drag[97]+0.9682458365518543*(f[10]*alpha_drag[97]+f[17]*alpha_drag[96])+0.3872983346207416*f[28]*alpha_drag[59]+0.4330127018922193*(f[42]*alpha_drag[51]+f[14]*alpha_drag[49]+f[28]*alpha_drag[48]); 
  out[43] += 1.677050983124842*(alpha_diffusion[11]*f[32]+alpha_diffusion[1]*f[15]+alpha_diffusion[0]*f[7])*rdvzSq4nu+0.9682458365518543*f[44]*alpha_drag[107]+0.8660254037844386*f[43]*alpha_drag[100]+0.9682458365518543*(f[7]*alpha_drag[100]+f[31]*alpha_drag[97]+f[18]*alpha_drag[96])+0.4330127018922193*(f[43]*alpha_drag[51]+f[41]*alpha_drag[49]+f[29]*alpha_drag[48]+alpha_drag[2]*f[43]+alpha_drag[1]*f[42]+alpha_drag[0]*f[30]); 
  out[44] += 0.276641667586244*f[32]*alpha_drag[107]+0.4330127018922193*(f[7]*alpha_drag[107]+f[44]*alpha_drag[100])+0.3872983346207416*f[15]*alpha_drag[97]+0.4330127018922193*f[32]*alpha_drag[96]+0.276641667586244*f[35]*alpha_drag[59]+0.4330127018922193*(f[9]*alpha_drag[59]+f[44]*alpha_drag[51])+0.3872983346207416*f[16]*alpha_drag[49]+0.4330127018922193*(f[35]*alpha_drag[48]+alpha_drag[2]*f[44])+(0.276641667586244*alpha_drag[11]+0.4330127018922193*alpha_drag[0])*f[37]+0.3872983346207416*alpha_drag[1]*f[17]+0.4330127018922193*f[10]*alpha_drag[11]; 
  out[45] += (1.5*alpha_diffusion[1]*f[37]+1.677050983124842*(alpha_diffusion[0]*f[17]+alpha_diffusion[1]*f[10])+1.5*alpha_diffusion[11]*f[17])*rdvxSq4nu+0.3872983346207416*f[33]*alpha_drag[107]+0.4330127018922193*(f[45]*alpha_drag[100]+f[22]*alpha_drag[97]+f[33]*alpha_drag[96])+0.3872983346207416*f[36]*alpha_drag[59]+0.4330127018922193*(f[45]*alpha_drag[51]+f[26]*alpha_drag[49]+f[36]*alpha_drag[48])+0.8660254037844386*(alpha_drag[2]*f[45]+alpha_drag[1]*f[44]+alpha_drag[11]*f[31])+0.9682458365518543*(alpha_drag[0]*f[31]+alpha_drag[1]*f[18]+alpha_drag[2]*f[17]); 
  out[46] += (1.5*alpha_diffusion[1]*f[35]+1.677050983124842*(alpha_diffusion[0]*f[16]+alpha_diffusion[1]*f[9])+1.5*alpha_diffusion[11]*f[16])*rdvySq4nu+0.3872983346207416*f[34]*alpha_drag[107]+0.4330127018922193*(f[46]*alpha_drag[100]+f[24]*alpha_drag[97]+f[34]*alpha_drag[96])+0.8660254037844386*f[31]*alpha_drag[59]+(0.8660254037844386*f[46]+0.9682458365518543*f[16])*alpha_drag[51]+0.8660254037844386*f[44]*alpha_drag[49]+0.9682458365518543*(f[18]*alpha_drag[49]+f[31]*alpha_drag[48])+0.4330127018922193*alpha_drag[2]*f[46]+0.3872983346207416*alpha_drag[11]*f[39]+0.4330127018922193*(alpha_drag[0]*f[39]+alpha_drag[1]*f[27]); 
  out[47] += (1.5*alpha_diffusion[1]*f[32]+1.677050983124842*(alpha_diffusion[0]*f[15]+alpha_diffusion[1]*f[7])+1.5*alpha_diffusion[11]*f[15])*rdvzSq4nu+0.8660254037844386*f[31]*alpha_drag[107]+(0.8660254037844386*f[47]+0.9682458365518543*f[15])*alpha_drag[100]+0.8660254037844386*f[44]*alpha_drag[97]+0.9682458365518543*(f[18]*alpha_drag[97]+f[31]*alpha_drag[96])+0.3872983346207416*f[41]*alpha_drag[59]+0.4330127018922193*(f[47]*alpha_drag[51]+f[29]*alpha_drag[49]+f[41]*alpha_drag[48]+alpha_drag[2]*f[47])+0.3872983346207416*alpha_drag[11]*f[42]+0.4330127018922193*(alpha_drag[0]*f[42]+alpha_drag[1]*f[30]); 

  return alpha_mid; 

} 
double VmLBOconstNuVol1x3vSerP3(const double *w, const double *dxv, const double nu, const double *u, const double *vtSq, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates. 
  // dxv[4]: Cell spacing. 
  // nu:     diffusion coefficient (collisionality). 
  // u:    bulk velocity. 
  // vtSq: thermal speed squared. 
  // f:      Input distribution function.
  // out:    Incremented output 
  const double rdvx2nu = 2.0*nu/dxv[1]; 
  const double rdvxSq4nu = 4.0*nu/(dxv[1]*dxv[1]); 
  const double rdvy2nu = 2.0*nu/dxv[2]; 
  const double rdvySq4nu = 4.0*nu/(dxv[2]*dxv[2]); 
  const double rdvz2nu = 2.0*nu/dxv[3]; 
  const double rdvzSq4nu = 4.0*nu/(dxv[3]*dxv[3]); 

  double alpha_mid = 0.0; 
  double alpha_drag[240]; 
  double alpha_diffusion[80]; 

  // Expand rdv2nu*(vx-ux) in phase basis.
  alpha_drag[0] = (2.828427124746191*u[0]-4.0*w[1])*rdvx2nu; 
  alpha_drag[1] = 2.828427124746191*u[1]*rdvx2nu; 
  alpha_drag[2] = -1.154700538379252*dxv[1]*rdvx2nu; 
  alpha_drag[11] = 2.828427124746191*u[2]*rdvx2nu; 
  alpha_drag[31] = 2.828427124746191*u[3]*rdvx2nu; 

  // x-drag contribution to the midpoint value used for CFL.
  alpha_mid += std::abs(0.125*alpha_drag[0]-0.1397542485937369*alpha_drag[11]); 

  // Expand rdv2nu*(vy-uy) in phase basis.
  alpha_drag[80] = (2.828427124746191*u[4]-4.0*w[2])*rdvy2nu; 
  alpha_drag[81] = 2.828427124746191*u[5]*rdvy2nu; 
  alpha_drag[83] = -1.154700538379252*dxv[2]*rdvy2nu; 
  alpha_drag[91] = 2.828427124746191*u[6]*rdvy2nu; 
  alpha_drag[111] = 2.828427124746191*u[7]*rdvy2nu; 

  // y-drag contribution to the midpoint value used for CFL.
  alpha_mid += std::abs(0.125*alpha_drag[80]-0.1397542485937369*alpha_drag[91]); 

  // Expand rdv2nu*(vz-uz) in phase basis.
  alpha_drag[160] = (2.828427124746191*u[8]-4.0*w[3])*rdvz2nu; 
  alpha_drag[161] = 2.828427124746191*u[9]*rdvz2nu; 
  alpha_drag[164] = -1.154700538379252*dxv[3]*rdvz2nu; 
  alpha_drag[171] = 2.828427124746191*u[10]*rdvz2nu; 
  alpha_drag[191] = 2.828427124746191*u[11]*rdvz2nu; 

  // z-drag contribution to the midpoint value used for CFL.
  alpha_mid += std::abs(0.125*alpha_drag[160]-0.1397542485937369*alpha_drag[171]); 

  // Expand nu*vthSq in phase basis.
  alpha_diffusion[0] = 2.828427124746191*vtSq[0]; 
  alpha_diffusion[1] = 2.828427124746191*vtSq[1]; 
  alpha_diffusion[11] = 2.828427124746191*vtSq[2]; 
  alpha_diffusion[31] = 2.828427124746191*vtSq[3]; 

  // Diffusion contribution to the midpoint value used for CFL.
  alpha_mid += std::abs((0.5714285714285714*alpha_diffusion[0]-0.6388765649999399*alpha_diffusion[11])*rdvxSq4nu); 
  alpha_mid += std::abs((0.5714285714285714*alpha_diffusion[0]-0.6388765649999399*alpha_diffusion[11])*rdvySq4nu); 
  alpha_mid += std::abs((0.5714285714285714*alpha_diffusion[0]-0.6388765649999399*alpha_diffusion[11])*rdvzSq4nu); 

  // Put together updates due to drag and diffusion terms.
  out[2] += 0.4330127018922193*(alpha_drag[31]*f[31]+alpha_drag[11]*f[11]+alpha_drag[2]*f[2]+alpha_drag[1]*f[1]+alpha_drag[0]*f[0]); 
  out[3] += 0.4330127018922193*(f[31]*alpha_drag[111]+f[11]*alpha_drag[91]+f[3]*alpha_drag[83]+f[1]*alpha_drag[81]+f[0]*alpha_drag[80]); 
  out[4] += 0.4330127018922193*(f[31]*alpha_drag[191]+f[11]*alpha_drag[171]+f[4]*alpha_drag[164]+f[1]*alpha_drag[161]+f[0]*alpha_drag[160]); 
  out[5] += 0.3803194146278324*(alpha_drag[11]*f[31]+f[11]*alpha_drag[31])+0.3872983346207416*(alpha_drag[1]*f[11]+f[1]*alpha_drag[11])+0.4330127018922193*(alpha_drag[2]*f[5]+alpha_drag[0]*f[1]+f[0]*alpha_drag[1]); 
  out[6] += 0.3803194146278324*f[11]*alpha_drag[111]+(0.3803194146278324*f[31]+0.3872983346207416*f[1])*alpha_drag[91]+0.4330127018922193*f[6]*alpha_drag[83]+0.3872983346207416*f[11]*alpha_drag[81]+0.4330127018922193*(f[0]*alpha_drag[81]+f[1]*alpha_drag[80]); 
  out[7] += 0.4330127018922193*(f[48]*alpha_drag[111]+f[19]*alpha_drag[91]+f[7]*alpha_drag[83]+f[5]*alpha_drag[81]+f[2]*alpha_drag[80]+alpha_drag[31]*f[50]+alpha_drag[11]*f[21]+alpha_drag[2]*f[7]+alpha_drag[1]*f[6]+alpha_drag[0]*f[3]); 
  out[8] += 0.3803194146278324*f[11]*alpha_drag[191]+(0.3803194146278324*f[31]+0.3872983346207416*f[1])*alpha_drag[171]+0.4330127018922193*f[8]*alpha_drag[164]+0.3872983346207416*f[11]*alpha_drag[161]+0.4330127018922193*(f[0]*alpha_drag[161]+f[1]*alpha_drag[160]); 
  out[9] += 0.4330127018922193*(f[48]*alpha_drag[191]+f[19]*alpha_drag[171]+f[9]*alpha_drag[164]+f[5]*alpha_drag[161]+f[2]*alpha_drag[160]+alpha_drag[31]*f[54]+alpha_drag[11]*f[25]+alpha_drag[2]*f[9]+alpha_drag[1]*f[8]+alpha_drag[0]*f[4]); 
  out[10] += 0.4330127018922193*(f[50]*alpha_drag[191]+f[21]*alpha_drag[171]+f[10]*alpha_drag[164]+f[6]*alpha_drag[161]+f[3]*alpha_drag[160]+f[54]*alpha_drag[111]+f[25]*alpha_drag[91]+f[10]*alpha_drag[83]+f[8]*alpha_drag[81]+f[4]*alpha_drag[80]); 
  out[12] += 1.677050983124842*(alpha_diffusion[31]*f[31]+alpha_diffusion[11]*f[11]+alpha_diffusion[1]*f[1]+alpha_diffusion[0]*f[0])*rdvxSq4nu+0.9682458365518543*(alpha_drag[31]*f[48]+alpha_drag[11]*f[19])+0.8660254037844386*alpha_drag[2]*f[12]+0.9682458365518543*(alpha_drag[1]*f[5]+alpha_drag[0]*f[2]+f[0]*alpha_drag[2]); 
  out[13] += 1.677050983124842*(alpha_diffusion[31]*f[31]+alpha_diffusion[11]*f[11]+alpha_diffusion[1]*f[1]+alpha_diffusion[0]*f[0])*rdvySq4nu+0.9682458365518543*(f[50]*alpha_drag[111]+f[21]*alpha_drag[91])+0.8660254037844386*f[13]*alpha_drag[83]+0.9682458365518543*(f[0]*alpha_drag[83]+f[6]*alpha_drag[81]+f[3]*alpha_drag[80]); 
  out[14] += 1.677050983124842*(alpha_diffusion[31]*f[31]+alpha_diffusion[11]*f[11]+alpha_diffusion[1]*f[1]+alpha_diffusion[0]*f[0])*rdvzSq4nu+0.9682458365518543*(f[54]*alpha_drag[191]+f[25]*alpha_drag[171])+0.8660254037844386*f[14]*alpha_drag[164]+0.9682458365518543*(f[0]*alpha_drag[164]+f[8]*alpha_drag[161]+f[4]*alpha_drag[160]); 
  out[15] += 0.3803194146278324*f[19]*alpha_drag[111]+(0.3803194146278324*f[48]+0.3872983346207416*f[5])*alpha_drag[91]+0.4330127018922193*f[15]*alpha_drag[83]+0.3872983346207416*f[19]*alpha_drag[81]+0.4330127018922193*(f[2]*alpha_drag[81]+f[5]*alpha_drag[80])+0.3803194146278324*alpha_drag[11]*f[50]+f[21]*(0.3803194146278324*alpha_drag[31]+0.3872983346207416*alpha_drag[1])+0.4330127018922193*alpha_drag[2]*f[15]+0.3872983346207416*f[6]*alpha_drag[11]+0.4330127018922193*(alpha_drag[0]*f[6]+alpha_drag[1]*f[3]); 
  out[16] += 0.3803194146278324*f[19]*alpha_drag[191]+(0.3803194146278324*f[48]+0.3872983346207416*f[5])*alpha_drag[171]+0.4330127018922193*f[16]*alpha_drag[164]+0.3872983346207416*f[19]*alpha_drag[161]+0.4330127018922193*(f[2]*alpha_drag[161]+f[5]*alpha_drag[160])+0.3803194146278324*alpha_drag[11]*f[54]+f[25]*(0.3803194146278324*alpha_drag[31]+0.3872983346207416*alpha_drag[1])+0.4330127018922193*alpha_drag[2]*f[16]+0.3872983346207416*f[8]*alpha_drag[11]+0.4330127018922193*(alpha_drag[0]*f[8]+alpha_drag[1]*f[4]); 
  out[17] += 0.3803194146278324*f[21]*alpha_drag[191]+(0.3803194146278324*f[50]+0.3872983346207416*f[6])*alpha_drag[171]+0.4330127018922193*f[17]*alpha_drag[164]+0.3872983346207416*f[21]*alpha_drag[161]+0.4330127018922193*(f[3]*alpha_drag[161]+f[6]*alpha_drag[160])+0.3803194146278324*f[25]*alpha_drag[111]+(0.3803194146278324*f[54]+0.3872983346207416*f[8])*alpha_drag[91]+0.4330127018922193*f[17]*alpha_drag[83]+0.3872983346207416*f[25]*alpha_drag[81]+0.4330127018922193*(f[4]*alpha_drag[81]+f[8]*alpha_drag[80]); 
  out[18] += 0.4330127018922193*(f[64]*alpha_drag[191]+f[36]*alpha_drag[171]+f[18]*alpha_drag[164]+f[15]*alpha_drag[161]+f[7]*alpha_drag[160]+f[67]*alpha_drag[111]+f[39]*alpha_drag[91]+f[18]*alpha_drag[83]+f[16]*alpha_drag[81]+f[9]*alpha_drag[80]+alpha_drag[31]*f[69]+alpha_drag[11]*f[41]+alpha_drag[2]*f[18]+alpha_drag[1]*f[17]+alpha_drag[0]*f[10]); 
  out[19] += 0.2581988897471612*alpha_drag[31]*f[31]+0.3803194146278324*(alpha_drag[1]*f[31]+f[1]*alpha_drag[31])+0.4330127018922193*alpha_drag[2]*f[19]+0.276641667586244*alpha_drag[11]*f[11]+0.4330127018922193*(alpha_drag[0]*f[11]+f[0]*alpha_drag[11])+0.3872983346207416*alpha_drag[1]*f[1]; 
  out[20] += (1.472970759092948*(alpha_diffusion[11]*f[31]+f[11]*alpha_diffusion[31])+1.5*(alpha_diffusion[1]*f[11]+f[1]*alpha_diffusion[11])+1.677050983124842*(alpha_diffusion[0]*f[1]+f[0]*alpha_diffusion[1]))*rdvxSq4nu+0.8504200642707612*(alpha_drag[11]*f[48]+f[19]*alpha_drag[31])+0.8660254037844386*(alpha_drag[2]*f[20]+alpha_drag[1]*f[19]+f[5]*alpha_drag[11])+0.9682458365518543*(alpha_drag[0]*f[5]+alpha_drag[1]*f[2]+f[1]*alpha_drag[2]); 
  out[21] += (0.2581988897471612*f[31]+0.3803194146278324*f[1])*alpha_drag[111]+0.276641667586244*f[11]*alpha_drag[91]+0.4330127018922193*(f[0]*alpha_drag[91]+f[21]*alpha_drag[83])+(0.3803194146278324*f[31]+0.3872983346207416*f[1])*alpha_drag[81]+0.4330127018922193*f[11]*alpha_drag[80]; 
  out[22] += 1.677050983124842*(alpha_diffusion[31]*f[50]+alpha_diffusion[11]*f[21]+alpha_diffusion[1]*f[6]+alpha_diffusion[0]*f[3])*rdvxSq4nu+0.4330127018922193*(f[22]*alpha_drag[83]+f[20]*alpha_drag[81]+f[12]*alpha_drag[80])+0.9682458365518543*(alpha_drag[31]*f[64]+alpha_drag[11]*f[36])+0.8660254037844386*alpha_drag[2]*f[22]+0.9682458365518543*(alpha_drag[1]*f[15]+alpha_drag[0]*f[7]+alpha_drag[2]*f[3]); 
  out[23] += (1.472970759092948*(alpha_diffusion[11]*f[31]+f[11]*alpha_diffusion[31])+1.5*(alpha_diffusion[1]*f[11]+f[1]*alpha_diffusion[11])+1.677050983124842*(alpha_diffusion[0]*f[1]+f[0]*alpha_diffusion[1]))*rdvySq4nu+0.8504200642707612*f[21]*alpha_drag[111]+(0.8504200642707612*f[50]+0.8660254037844386*f[6])*alpha_drag[91]+(0.8660254037844386*f[23]+0.9682458365518543*f[1])*alpha_drag[83]+0.8660254037844386*f[21]*alpha_drag[81]+0.9682458365518543*(f[3]*alpha_drag[81]+f[6]*alpha_drag[80]); 
  out[24] += 1.677050983124842*(alpha_diffusion[31]*f[48]+alpha_diffusion[11]*f[19]+alpha_diffusion[1]*f[5]+alpha_diffusion[0]*f[2])*rdvySq4nu+0.9682458365518543*(f[64]*alpha_drag[111]+f[36]*alpha_drag[91])+0.8660254037844386*f[24]*alpha_drag[83]+0.9682458365518543*(f[2]*alpha_drag[83]+f[15]*alpha_drag[81]+f[7]*alpha_drag[80])+0.4330127018922193*(alpha_drag[2]*f[24]+alpha_drag[1]*f[23]+alpha_drag[0]*f[13]); 
  out[25] += (0.2581988897471612*f[31]+0.3803194146278324*f[1])*alpha_drag[191]+0.276641667586244*f[11]*alpha_drag[171]+0.4330127018922193*(f[0]*alpha_drag[171]+f[25]*alpha_drag[164])+(0.3803194146278324*f[31]+0.3872983346207416*f[1])*alpha_drag[161]+0.4330127018922193*f[11]*alpha_drag[160]; 
  out[26] += 1.677050983124842*(alpha_diffusion[31]*f[54]+alpha_diffusion[11]*f[25]+alpha_diffusion[1]*f[8]+alpha_diffusion[0]*f[4])*rdvxSq4nu+0.4330127018922193*(f[26]*alpha_drag[164]+f[20]*alpha_drag[161]+f[12]*alpha_drag[160])+0.9682458365518543*(alpha_drag[31]*f[67]+alpha_drag[11]*f[39])+0.8660254037844386*alpha_drag[2]*f[26]+0.9682458365518543*(alpha_drag[1]*f[16]+alpha_drag[0]*f[9]+alpha_drag[2]*f[4]); 
  out[27] += 1.677050983124842*(alpha_diffusion[31]*f[54]+alpha_diffusion[11]*f[25]+alpha_diffusion[1]*f[8]+alpha_diffusion[0]*f[4])*rdvySq4nu+0.4330127018922193*(f[27]*alpha_drag[164]+f[23]*alpha_drag[161]+f[13]*alpha_drag[160])+0.9682458365518543*(f[69]*alpha_drag[111]+f[41]*alpha_drag[91])+0.8660254037844386*f[27]*alpha_drag[83]+0.9682458365518543*(f[4]*alpha_drag[83]+f[17]*alpha_drag[81]+f[10]*alpha_drag[80]); 
  out[28] += (1.472970759092948*(alpha_diffusion[11]*f[31]+f[11]*alpha_diffusion[31])+1.5*(alpha_diffusion[1]*f[11]+f[1]*alpha_diffusion[11])+1.677050983124842*(alpha_diffusion[0]*f[1]+f[0]*alpha_diffusion[1]))*rdvzSq4nu+0.8504200642707612*f[25]*alpha_drag[191]+(0.8504200642707612*f[54]+0.8660254037844386*f[8])*alpha_drag[171]+(0.8660254037844386*f[28]+0.9682458365518543*f[1])*alpha_drag[164]+0.8660254037844386*f[25]*alpha_drag[161]+0.9682458365518543*(f[4]*alpha_drag[161]+f[8]*alpha_drag[160]); 
  out[29] += 1.677050983124842*(alpha_diffusion[31]*f[48]+alpha_diffusion[11]*f[19]+alpha_diffusion[1]*f[5]+alpha_diffusion[0]*f[2])*rdvzSq4nu+0.9682458365518543*(f[67]*alpha_drag[191]+f[39]*alpha_drag[171])+0.8660254037844386*f[29]*alpha_drag[164]+0.9682458365518543*(f[2]*alpha_drag[164]+f[16]*alpha_drag[161]+f[9]*alpha_drag[160])+0.4330127018922193*(alpha_drag[2]*f[29]+alpha_drag[1]*f[28]+alpha_drag[0]*f[14]); 
  out[30] += 1.677050983124842*(alpha_diffusion[31]*f[50]+alpha_diffusion[11]*f[21]+alpha_diffusion[1]*f[6]+alpha_diffusion[0]*f[3])*rdvzSq4nu+0.9682458365518543*(f[69]*alpha_drag[191]+f[41]*alpha_drag[171])+0.8660254037844386*f[30]*alpha_drag[164]+0.9682458365518543*(f[3]*alpha_drag[164]+f[17]*alpha_drag[161]+f[10]*alpha_drag[160])+0.4330127018922193*(f[30]*alpha_drag[83]+f[28]*alpha_drag[81]+f[14]*alpha_drag[80]); 
  out[32] += 5.7282196186948*(alpha_diffusion[31]*f[48]+alpha_diffusion[11]*f[19]+alpha_diffusion[1]*f[5]+alpha_diffusion[0]*f[2])*rdvxSq4nu+1.299038105676658*alpha_drag[2]*f[32]+0.6614378277661477*alpha_drag[31]*f[31]+1.479019945774904*(alpha_drag[1]*f[20]+alpha_drag[0]*f[12])+0.6614378277661477*alpha_drag[11]*f[11]+1.984313483298443*alpha_drag[2]*f[2]+0.6614378277661477*(alpha_drag[1]*f[1]+alpha_drag[0]*f[0]); 
  out[33] += 5.7282196186948*(alpha_diffusion[31]*f[50]+alpha_diffusion[11]*f[21]+alpha_diffusion[1]*f[6]+alpha_diffusion[0]*f[3])*rdvySq4nu+0.6614378277661477*(f[31]*alpha_drag[111]+f[11]*alpha_drag[91])+(1.299038105676658*f[33]+1.984313483298443*f[3])*alpha_drag[83]+(1.479019945774904*f[23]+0.6614378277661477*f[1])*alpha_drag[81]+(1.479019945774904*f[13]+0.6614378277661477*f[0])*alpha_drag[80]; 
  out[34] += 5.7282196186948*(alpha_diffusion[31]*f[54]+alpha_diffusion[11]*f[25]+alpha_diffusion[1]*f[8]+alpha_diffusion[0]*f[4])*rdvzSq4nu+0.6614378277661477*(f[31]*alpha_drag[191]+f[11]*alpha_drag[171])+(1.299038105676658*f[34]+1.984313483298443*f[4])*alpha_drag[164]+(1.479019945774904*f[28]+0.6614378277661477*f[1])*alpha_drag[161]+(1.479019945774904*f[14]+0.6614378277661477*f[0])*alpha_drag[160]; 
  out[35] += 0.3803194146278324*f[36]*alpha_drag[191]+(0.3803194146278324*f[64]+0.3872983346207416*f[15])*alpha_drag[171]+0.4330127018922193*f[35]*alpha_drag[164]+0.3872983346207416*f[36]*alpha_drag[161]+0.4330127018922193*(f[7]*alpha_drag[161]+f[15]*alpha_drag[160])+0.3803194146278324*f[39]*alpha_drag[111]+(0.3803194146278324*f[67]+0.3872983346207416*f[16])*alpha_drag[91]+0.4330127018922193*f[35]*alpha_drag[83]+0.3872983346207416*f[39]*alpha_drag[81]+0.4330127018922193*(f[9]*alpha_drag[81]+f[16]*alpha_drag[80])+0.3803194146278324*alpha_drag[11]*f[69]+(0.3803194146278324*alpha_drag[31]+0.3872983346207416*alpha_drag[1])*f[41]+0.4330127018922193*alpha_drag[2]*f[35]+0.3872983346207416*alpha_drag[11]*f[17]+0.4330127018922193*(alpha_drag[0]*f[17]+alpha_drag[1]*f[10]); 
  out[36] += (0.2581988897471612*f[48]+0.3803194146278324*f[5])*alpha_drag[111]+0.276641667586244*f[19]*alpha_drag[91]+0.4330127018922193*(f[2]*alpha_drag[91]+f[36]*alpha_drag[83])+(0.3803194146278324*f[48]+0.3872983346207416*f[5])*alpha_drag[81]+0.4330127018922193*f[19]*alpha_drag[80]+(0.2581988897471612*alpha_drag[31]+0.3803194146278324*alpha_drag[1])*f[50]+0.4330127018922193*alpha_drag[2]*f[36]+0.3803194146278324*f[6]*alpha_drag[31]+0.276641667586244*alpha_drag[11]*f[21]+0.4330127018922193*(alpha_drag[0]*f[21]+f[3]*alpha_drag[11])+0.3872983346207416*alpha_drag[1]*f[6]; 
  out[37] += (1.472970759092948*(alpha_diffusion[11]*f[50]+f[21]*alpha_diffusion[31])+1.5*(alpha_diffusion[1]*f[21]+f[6]*alpha_diffusion[11])+1.677050983124842*(alpha_diffusion[0]*f[6]+alpha_diffusion[1]*f[3]))*rdvxSq4nu+0.3872983346207416*f[20]*alpha_drag[91]+0.4330127018922193*(f[37]*alpha_drag[83]+f[12]*alpha_drag[81]+f[20]*alpha_drag[80])+0.8504200642707612*alpha_drag[11]*f[64]+0.8660254037844386*alpha_drag[2]*f[37]+0.8504200642707612*alpha_drag[31]*f[36]+0.8660254037844386*(alpha_drag[1]*f[36]+alpha_drag[11]*f[15])+0.9682458365518543*(alpha_drag[0]*f[15]+alpha_drag[1]*f[7]+alpha_drag[2]*f[6]); 
  out[38] += (1.472970759092948*(alpha_diffusion[11]*f[48]+f[19]*alpha_diffusion[31])+1.5*(alpha_diffusion[1]*f[19]+f[5]*alpha_diffusion[11])+1.677050983124842*(alpha_diffusion[0]*f[5]+alpha_diffusion[1]*f[2]))*rdvySq4nu+0.8504200642707612*f[36]*alpha_drag[111]+(0.8504200642707612*f[64]+0.8660254037844386*f[15])*alpha_drag[91]+(0.8660254037844386*f[38]+0.9682458365518543*f[5])*alpha_drag[83]+0.8660254037844386*f[36]*alpha_drag[81]+0.9682458365518543*(f[7]*alpha_drag[81]+f[15]*alpha_drag[80])+0.4330127018922193*alpha_drag[2]*f[38]+0.3872983346207416*alpha_drag[11]*f[23]+0.4330127018922193*(alpha_drag[0]*f[23]+alpha_drag[1]*f[13]); 
  out[39] += (0.2581988897471612*f[48]+0.3803194146278324*f[5])*alpha_drag[191]+0.276641667586244*f[19]*alpha_drag[171]+0.4330127018922193*(f[2]*alpha_drag[171]+f[39]*alpha_drag[164])+(0.3803194146278324*f[48]+0.3872983346207416*f[5])*alpha_drag[161]+0.4330127018922193*f[19]*alpha_drag[160]+(0.2581988897471612*alpha_drag[31]+0.3803194146278324*alpha_drag[1])*f[54]+0.4330127018922193*alpha_drag[2]*f[39]+0.3803194146278324*f[8]*alpha_drag[31]+0.276641667586244*alpha_drag[11]*f[25]+0.4330127018922193*(alpha_drag[0]*f[25]+f[4]*alpha_drag[11])+0.3872983346207416*alpha_drag[1]*f[8]; 
  out[40] += (1.472970759092948*(alpha_diffusion[11]*f[54]+f[25]*alpha_diffusion[31])+1.5*(alpha_diffusion[1]*f[25]+f[8]*alpha_diffusion[11])+1.677050983124842*(alpha_diffusion[0]*f[8]+alpha_diffusion[1]*f[4]))*rdvxSq4nu+0.3872983346207416*f[20]*alpha_drag[171]+0.4330127018922193*(f[40]*alpha_drag[164]+f[12]*alpha_drag[161]+f[20]*alpha_drag[160])+0.8504200642707612*alpha_drag[11]*f[67]+0.8660254037844386*alpha_drag[2]*f[40]+0.8504200642707612*alpha_drag[31]*f[39]+0.8660254037844386*(alpha_drag[1]*f[39]+alpha_drag[11]*f[16])+0.9682458365518543*(alpha_drag[0]*f[16]+alpha_drag[1]*f[9]+alpha_drag[2]*f[8]); 
  out[41] += (0.2581988897471612*f[50]+0.3803194146278324*f[6])*alpha_drag[191]+0.276641667586244*f[21]*alpha_drag[171]+0.4330127018922193*(f[3]*alpha_drag[171]+f[41]*alpha_drag[164])+(0.3803194146278324*f[50]+0.3872983346207416*f[6])*alpha_drag[161]+0.4330127018922193*f[21]*alpha_drag[160]+(0.2581988897471612*f[54]+0.3803194146278324*f[8])*alpha_drag[111]+0.276641667586244*f[25]*alpha_drag[91]+0.4330127018922193*(f[4]*alpha_drag[91]+f[41]*alpha_drag[83])+(0.3803194146278324*f[54]+0.3872983346207416*f[8])*alpha_drag[81]+0.4330127018922193*f[25]*alpha_drag[80]; 
  out[42] += 1.677050983124842*(alpha_diffusion[31]*f[69]+alpha_diffusion[11]*f[41]+alpha_diffusion[1]*f[17]+alpha_diffusion[0]*f[10])*rdvxSq4nu+0.4330127018922193*(f[42]*alpha_drag[164]+f[37]*alpha_drag[161]+f[22]*alpha_drag[160]+f[42]*alpha_drag[83]+f[40]*alpha_drag[81]+f[26]*alpha_drag[80])+0.9682458365518543*(alpha_drag[31]*f[76]+alpha_drag[11]*f[60])+0.8660254037844386*alpha_drag[2]*f[42]+0.9682458365518543*(alpha_drag[1]*f[35]+alpha_drag[0]*f[18]+alpha_drag[2]*f[10]); 
  out[43] += (1.472970759092948*(alpha_diffusion[11]*f[54]+f[25]*alpha_diffusion[31])+1.5*(alpha_diffusion[1]*f[25]+f[8]*alpha_diffusion[11])+1.677050983124842*(alpha_diffusion[0]*f[8]+alpha_diffusion[1]*f[4]))*rdvySq4nu+0.3872983346207416*f[23]*alpha_drag[171]+0.4330127018922193*(f[43]*alpha_drag[164]+f[13]*alpha_drag[161]+f[23]*alpha_drag[160])+0.8504200642707612*f[41]*alpha_drag[111]+(0.8504200642707612*f[69]+0.8660254037844386*f[17])*alpha_drag[91]+(0.8660254037844386*f[43]+0.9682458365518543*f[8])*alpha_drag[83]+0.8660254037844386*f[41]*alpha_drag[81]+0.9682458365518543*(f[10]*alpha_drag[81]+f[17]*alpha_drag[80]); 
  out[44] += 1.677050983124842*(alpha_diffusion[31]*f[67]+alpha_diffusion[11]*f[39]+alpha_diffusion[1]*f[16]+alpha_diffusion[0]*f[9])*rdvySq4nu+0.4330127018922193*(f[44]*alpha_drag[164]+f[38]*alpha_drag[161]+f[24]*alpha_drag[160])+0.9682458365518543*(f[76]*alpha_drag[111]+f[60]*alpha_drag[91])+0.8660254037844386*f[44]*alpha_drag[83]+0.9682458365518543*(f[9]*alpha_drag[83]+f[35]*alpha_drag[81]+f[18]*alpha_drag[80])+0.4330127018922193*(alpha_drag[2]*f[44]+alpha_drag[1]*f[43]+alpha_drag[0]*f[27]); 
  out[45] += (1.472970759092948*(alpha_diffusion[11]*f[48]+f[19]*alpha_diffusion[31])+1.5*(alpha_diffusion[1]*f[19]+f[5]*alpha_diffusion[11])+1.677050983124842*(alpha_diffusion[0]*f[5]+alpha_diffusion[1]*f[2]))*rdvzSq4nu+0.8504200642707612*f[39]*alpha_drag[191]+(0.8504200642707612*f[67]+0.8660254037844386*f[16])*alpha_drag[171]+(0.8660254037844386*f[45]+0.9682458365518543*f[5])*alpha_drag[164]+0.8660254037844386*f[39]*alpha_drag[161]+0.9682458365518543*(f[9]*alpha_drag[161]+f[16]*alpha_drag[160])+0.4330127018922193*alpha_drag[2]*f[45]+0.3872983346207416*alpha_drag[11]*f[28]+0.4330127018922193*(alpha_drag[0]*f[28]+alpha_drag[1]*f[14]); 
  out[46] += (1.472970759092948*(alpha_diffusion[11]*f[50]+f[21]*alpha_diffusion[31])+1.5*(alpha_diffusion[1]*f[21]+f[6]*alpha_diffusion[11])+1.677050983124842*(alpha_diffusion[0]*f[6]+alpha_diffusion[1]*f[3]))*rdvzSq4nu+0.8504200642707612*f[41]*alpha_drag[191]+(0.8504200642707612*f[69]+0.8660254037844386*f[17])*alpha_drag[171]+(0.8660254037844386*f[46]+0.9682458365518543*f[6])*alpha_drag[164]+0.8660254037844386*f[41]*alpha_drag[161]+0.9682458365518543*(f[10]*alpha_drag[161]+f[17]*alpha_drag[160])+0.3872983346207416*f[28]*alpha_drag[91]+0.4330127018922193*(f[46]*alpha_drag[83]+f[14]*alpha_drag[81]+f[28]*alpha_drag[80]); 
  out[47] += 1.677050983124842*(alpha_diffusion[31]*f[64]+alpha_diffusion[11]*f[36]+alpha_diffusion[1]*f[15]+alpha_diffusion[0]*f[7])*rdvzSq4nu+0.9682458365518543*(f[76]*alpha_drag[191]+f[60]*alpha_drag[171])+0.8660254037844386*f[47]*alpha_drag[164]+0.9682458365518543*(f[7]*alpha_drag[164]+f[35]*alpha_drag[161]+f[18]*alpha_drag[160])+0.4330127018922193*(f[47]*alpha_drag[83]+f[45]*alpha_drag[81]+f[29]*alpha_drag[80]+alpha_drag[2]*f[47]+alpha_drag[1]*f[46]+alpha_drag[0]*f[30]); 
  out[48] += 0.4330127018922193*alpha_drag[2]*f[48]+(0.2581988897471612*alpha_drag[11]+0.4330127018922193*alpha_drag[0])*f[31]+(0.2581988897471612*f[11]+0.4330127018922193*f[0])*alpha_drag[31]+0.3803194146278324*(alpha_drag[1]*f[11]+f[1]*alpha_drag[11]); 
  out[49] += (5.031152949374527*(alpha_diffusion[11]*f[48]+f[19]*alpha_diffusion[31])+5.1234753829798*(alpha_diffusion[1]*f[19]+f[5]*alpha_diffusion[11])+5.7282196186948*(alpha_diffusion[0]*f[5]+alpha_diffusion[1]*f[2]))*rdvxSq4nu+1.299038105676658*alpha_drag[2]*f[49]+0.5809475019311124*(alpha_drag[11]*f[31]+f[11]*alpha_drag[31])+1.322875655532295*alpha_drag[11]*f[20]+1.479019945774904*(alpha_drag[0]*f[20]+alpha_drag[1]*f[12])+0.5916079783099616*(alpha_drag[1]*f[11]+f[1]*alpha_drag[11])+1.984313483298443*alpha_drag[2]*f[5]+0.6614378277661477*(alpha_drag[0]*f[1]+f[0]*alpha_drag[1]); 
  out[50] += (0.2581988897471612*f[11]+0.4330127018922193*f[0])*alpha_drag[111]+(0.2581988897471612*f[31]+0.3803194146278324*f[1])*alpha_drag[91]+0.4330127018922193*f[50]*alpha_drag[83]+0.3803194146278324*f[11]*alpha_drag[81]+0.4330127018922193*f[31]*alpha_drag[80]; 
  out[51] += 5.7282196186948*(alpha_diffusion[31]*f[64]+alpha_diffusion[11]*f[36]+alpha_diffusion[1]*f[15]+alpha_diffusion[0]*f[7])*rdvxSq4nu+0.4330127018922193*(f[51]*alpha_drag[83]+f[49]*alpha_drag[81]+f[32]*alpha_drag[80])+1.299038105676658*alpha_drag[2]*f[51]+0.6614378277661477*alpha_drag[31]*f[50]+1.479019945774904*(alpha_drag[1]*f[37]+alpha_drag[0]*f[22])+0.6614378277661477*alpha_drag[11]*f[21]+1.984313483298443*alpha_drag[2]*f[7]+0.6614378277661477*(alpha_drag[1]*f[6]+alpha_drag[0]*f[3]); 
  out[52] += (5.031152949374527*(alpha_diffusion[11]*f[50]+f[21]*alpha_diffusion[31])+5.1234753829798*(alpha_diffusion[1]*f[21]+f[6]*alpha_diffusion[11])+5.7282196186948*(alpha_diffusion[0]*f[6]+alpha_diffusion[1]*f[3]))*rdvySq4nu+0.5809475019311124*f[11]*alpha_drag[111]+(0.5809475019311124*f[31]+1.322875655532295*f[23]+0.5916079783099616*f[1])*alpha_drag[91]+(1.299038105676658*f[52]+1.984313483298443*f[6])*alpha_drag[83]+(1.479019945774904*f[13]+0.5916079783099616*f[11]+0.6614378277661477*f[0])*alpha_drag[81]+(1.479019945774904*f[23]+0.6614378277661477*f[1])*alpha_drag[80]; 
  out[53] += 5.7282196186948*(alpha_diffusion[31]*f[64]+alpha_diffusion[11]*f[36]+alpha_diffusion[1]*f[15]+alpha_diffusion[0]*f[7])*rdvySq4nu+0.6614378277661477*(f[48]*alpha_drag[111]+f[19]*alpha_drag[91])+(1.299038105676658*f[53]+1.984313483298443*f[7])*alpha_drag[83]+(1.479019945774904*f[38]+0.6614378277661477*f[5])*alpha_drag[81]+(1.479019945774904*f[24]+0.6614378277661477*f[2])*alpha_drag[80]+0.4330127018922193*(alpha_drag[2]*f[53]+alpha_drag[1]*f[52]+alpha_drag[0]*f[33]); 
  out[54] += (0.2581988897471612*f[11]+0.4330127018922193*f[0])*alpha_drag[191]+(0.2581988897471612*f[31]+0.3803194146278324*f[1])*alpha_drag[171]+0.4330127018922193*f[54]*alpha_drag[164]+0.3803194146278324*f[11]*alpha_drag[161]+0.4330127018922193*f[31]*alpha_drag[160]; 
  out[55] += 5.7282196186948*(alpha_diffusion[31]*f[67]+alpha_diffusion[11]*f[39]+alpha_diffusion[1]*f[16]+alpha_diffusion[0]*f[9])*rdvxSq4nu+0.4330127018922193*(f[55]*alpha_drag[164]+f[49]*alpha_drag[161]+f[32]*alpha_drag[160])+1.299038105676658*alpha_drag[2]*f[55]+0.6614378277661477*alpha_drag[31]*f[54]+1.479019945774904*(alpha_drag[1]*f[40]+alpha_drag[0]*f[26])+0.6614378277661477*alpha_drag[11]*f[25]+1.984313483298443*alpha_drag[2]*f[9]+0.6614378277661477*(alpha_drag[1]*f[8]+alpha_drag[0]*f[4]); 
  out[56] += 5.7282196186948*(alpha_diffusion[31]*f[69]+alpha_diffusion[11]*f[41]+alpha_diffusion[1]*f[17]+alpha_diffusion[0]*f[10])*rdvySq4nu+0.4330127018922193*(f[56]*alpha_drag[164]+f[52]*alpha_drag[161]+f[33]*alpha_drag[160])+0.6614378277661477*(f[54]*alpha_drag[111]+f[25]*alpha_drag[91])+(1.299038105676658*f[56]+1.984313483298443*f[10])*alpha_drag[83]+(1.479019945774904*f[43]+0.6614378277661477*f[8])*alpha_drag[81]+(1.479019945774904*f[27]+0.6614378277661477*f[4])*alpha_drag[80]; 
  out[57] += (5.031152949374527*(alpha_diffusion[11]*f[54]+f[25]*alpha_diffusion[31])+5.1234753829798*(alpha_diffusion[1]*f[25]+f[8]*alpha_diffusion[11])+5.7282196186948*(alpha_diffusion[0]*f[8]+alpha_diffusion[1]*f[4]))*rdvzSq4nu+0.5809475019311124*f[11]*alpha_drag[191]+(0.5809475019311124*f[31]+1.322875655532295*f[28]+0.5916079783099616*f[1])*alpha_drag[171]+(1.299038105676658*f[57]+1.984313483298443*f[8])*alpha_drag[164]+(1.479019945774904*f[14]+0.5916079783099616*f[11]+0.6614378277661477*f[0])*alpha_drag[161]+(1.479019945774904*f[28]+0.6614378277661477*f[1])*alpha_drag[160]; 
  out[58] += 5.7282196186948*(alpha_diffusion[31]*f[67]+alpha_diffusion[11]*f[39]+alpha_diffusion[1]*f[16]+alpha_diffusion[0]*f[9])*rdvzSq4nu+0.6614378277661477*(f[48]*alpha_drag[191]+f[19]*alpha_drag[171])+(1.299038105676658*f[58]+1.984313483298443*f[9])*alpha_drag[164]+(1.479019945774904*f[45]+0.6614378277661477*f[5])*alpha_drag[161]+(1.479019945774904*f[29]+0.6614378277661477*f[2])*alpha_drag[160]+0.4330127018922193*(alpha_drag[2]*f[58]+alpha_drag[1]*f[57]+alpha_drag[0]*f[34]); 
  out[59] += 5.7282196186948*(alpha_diffusion[31]*f[69]+alpha_diffusion[11]*f[41]+alpha_diffusion[1]*f[17]+alpha_diffusion[0]*f[10])*rdvzSq4nu+0.6614378277661477*(f[50]*alpha_drag[191]+f[21]*alpha_drag[171])+(1.299038105676658*f[59]+1.984313483298443*f[10])*alpha_drag[164]+(1.479019945774904*f[46]+0.6614378277661477*f[6])*alpha_drag[161]+(1.479019945774904*f[30]+0.6614378277661477*f[3])*alpha_drag[160]+0.4330127018922193*(f[59]*alpha_drag[83]+f[57]*alpha_drag[81]+f[34]*alpha_drag[80]); 
  out[60] += (0.2581988897471612*f[64]+0.3803194146278324*f[15])*alpha_drag[191]+0.276641667586244*f[36]*alpha_drag[171]+0.4330127018922193*(f[7]*alpha_drag[171]+f[60]*alpha_drag[164])+(0.3803194146278324*f[64]+0.3872983346207416*f[15])*alpha_drag[161]+0.4330127018922193*f[36]*alpha_drag[160]+(0.2581988897471612*f[67]+0.3803194146278324*f[16])*alpha_drag[111]+0.276641667586244*f[39]*alpha_drag[91]+0.4330127018922193*(f[9]*alpha_drag[91]+f[60]*alpha_drag[83])+(0.3803194146278324*f[67]+0.3872983346207416*f[16])*alpha_drag[81]+0.4330127018922193*f[39]*alpha_drag[80]+(0.2581988897471612*alpha_drag[31]+0.3803194146278324*alpha_drag[1])*f[69]+0.4330127018922193*alpha_drag[2]*f[60]+(0.276641667586244*alpha_drag[11]+0.4330127018922193*alpha_drag[0])*f[41]+f[17]*(0.3803194146278324*alpha_drag[31]+0.3872983346207416*alpha_drag[1])+0.4330127018922193*f[10]*alpha_drag[11]; 
  out[61] += (1.472970759092948*alpha_diffusion[11]*f[69]+(1.472970759092948*alpha_diffusion[31]+1.5*alpha_diffusion[1])*f[41]+1.677050983124842*(alpha_diffusion[0]*f[17]+alpha_diffusion[1]*f[10])+1.5*alpha_diffusion[11]*f[17])*rdvxSq4nu+0.3872983346207416*f[37]*alpha_drag[171]+0.4330127018922193*(f[61]*alpha_drag[164]+f[22]*alpha_drag[161]+f[37]*alpha_drag[160])+0.3872983346207416*f[40]*alpha_drag[91]+0.4330127018922193*(f[61]*alpha_drag[83]+f[26]*alpha_drag[81]+f[40]*alpha_drag[80])+0.8504200642707612*alpha_drag[11]*f[76]+0.8660254037844386*alpha_drag[2]*f[61]+0.8504200642707612*alpha_drag[31]*f[60]+0.8660254037844386*(alpha_drag[1]*f[60]+alpha_drag[11]*f[35])+0.9682458365518543*(alpha_drag[0]*f[35]+alpha_drag[1]*f[18]+alpha_drag[2]*f[17]); 
  out[62] += (1.472970759092948*alpha_diffusion[11]*f[67]+(1.472970759092948*alpha_diffusion[31]+1.5*alpha_diffusion[1])*f[39]+1.677050983124842*(alpha_diffusion[0]*f[16]+alpha_diffusion[1]*f[9])+1.5*alpha_diffusion[11]*f[16])*rdvySq4nu+0.3872983346207416*f[38]*alpha_drag[171]+0.4330127018922193*(f[62]*alpha_drag[164]+f[24]*alpha_drag[161]+f[38]*alpha_drag[160])+0.8504200642707612*f[60]*alpha_drag[111]+(0.8504200642707612*f[76]+0.8660254037844386*f[35])*alpha_drag[91]+(0.8660254037844386*f[62]+0.9682458365518543*f[16])*alpha_drag[83]+0.8660254037844386*f[60]*alpha_drag[81]+0.9682458365518543*(f[18]*alpha_drag[81]+f[35]*alpha_drag[80])+0.4330127018922193*alpha_drag[2]*f[62]+0.3872983346207416*alpha_drag[11]*f[43]+0.4330127018922193*(alpha_drag[0]*f[43]+alpha_drag[1]*f[27]); 
  out[63] += (1.472970759092948*alpha_diffusion[11]*f[64]+(1.472970759092948*alpha_diffusion[31]+1.5*alpha_diffusion[1])*f[36]+1.677050983124842*(alpha_diffusion[0]*f[15]+alpha_diffusion[1]*f[7])+1.5*alpha_diffusion[11]*f[15])*rdvzSq4nu+0.8504200642707612*f[60]*alpha_drag[191]+(0.8504200642707612*f[76]+0.8660254037844386*f[35])*alpha_drag[171]+(0.8660254037844386*f[63]+0.9682458365518543*f[15])*alpha_drag[164]+0.8660254037844386*f[60]*alpha_drag[161]+0.9682458365518543*(f[18]*alpha_drag[161]+f[35]*alpha_drag[160])+0.3872983346207416*f[45]*alpha_drag[91]+0.4330127018922193*(f[63]*alpha_drag[83]+f[29]*alpha_drag[81]+f[45]*alpha_drag[80]+alpha_drag[2]*f[63])+0.3872983346207416*alpha_drag[11]*f[46]+0.4330127018922193*(alpha_drag[0]*f[46]+alpha_drag[1]*f[30]); 
  out[64] += (0.2581988897471612*f[19]+0.4330127018922193*f[2])*alpha_drag[111]+(0.2581988897471612*f[48]+0.3803194146278324*f[5])*alpha_drag[91]+0.4330127018922193*f[64]*alpha_drag[83]+0.3803194146278324*f[19]*alpha_drag[81]+0.4330127018922193*(f[48]*alpha_drag[80]+alpha_drag[2]*f[64])+(0.2581988897471612*alpha_drag[11]+0.4330127018922193*alpha_drag[0])*f[50]+(0.2581988897471612*f[21]+0.4330127018922193*f[3])*alpha_drag[31]+0.3803194146278324*(alpha_drag[1]*f[21]+f[6]*alpha_drag[11]); 
  out[65] += (5.031152949374527*alpha_diffusion[11]*f[64]+(5.031152949374527*alpha_diffusion[31]+5.1234753829798*alpha_diffusion[1])*f[36]+5.7282196186948*(alpha_diffusion[0]*f[15]+alpha_diffusion[1]*f[7])+5.1234753829798*alpha_diffusion[11]*f[15])*rdvxSq4nu+0.3872983346207416*f[49]*alpha_drag[91]+0.4330127018922193*(f[65]*alpha_drag[83]+f[32]*alpha_drag[81]+f[49]*alpha_drag[80])+1.299038105676658*alpha_drag[2]*f[65]+0.5809475019311124*alpha_drag[11]*f[50]+(1.322875655532295*alpha_drag[11]+1.479019945774904*alpha_drag[0])*f[37]+0.5809475019311124*f[21]*alpha_drag[31]+alpha_drag[1]*(1.479019945774904*f[22]+0.5916079783099616*f[21])+1.984313483298443*alpha_drag[2]*f[15]+0.5916079783099616*f[6]*alpha_drag[11]+0.6614378277661477*(alpha_drag[0]*f[6]+alpha_drag[1]*f[3]); 
  out[66] += (5.031152949374527*alpha_diffusion[11]*f[64]+(5.031152949374527*alpha_diffusion[31]+5.1234753829798*alpha_diffusion[1])*f[36]+5.7282196186948*(alpha_diffusion[0]*f[15]+alpha_diffusion[1]*f[7])+5.1234753829798*alpha_diffusion[11]*f[15])*rdvySq4nu+0.5809475019311124*f[19]*alpha_drag[111]+(0.5809475019311124*f[48]+1.322875655532295*f[38]+0.5916079783099616*f[5])*alpha_drag[91]+(1.299038105676658*f[66]+1.984313483298443*f[15])*alpha_drag[83]+(1.479019945774904*f[24]+0.5916079783099616*f[19]+0.6614378277661477*f[2])*alpha_drag[81]+(1.479019945774904*f[38]+0.6614378277661477*f[5])*alpha_drag[80]+0.4330127018922193*alpha_drag[2]*f[66]+0.3872983346207416*alpha_drag[11]*f[52]+0.4330127018922193*(alpha_drag[0]*f[52]+alpha_drag[1]*f[33]); 
  out[67] += (0.2581988897471612*f[19]+0.4330127018922193*f[2])*alpha_drag[191]+(0.2581988897471612*f[48]+0.3803194146278324*f[5])*alpha_drag[171]+0.4330127018922193*f[67]*alpha_drag[164]+0.3803194146278324*f[19]*alpha_drag[161]+0.4330127018922193*(f[48]*alpha_drag[160]+alpha_drag[2]*f[67])+(0.2581988897471612*alpha_drag[11]+0.4330127018922193*alpha_drag[0])*f[54]+(0.2581988897471612*f[25]+0.4330127018922193*f[4])*alpha_drag[31]+0.3803194146278324*(alpha_drag[1]*f[25]+f[8]*alpha_drag[11]); 
  out[68] += (5.031152949374527*alpha_diffusion[11]*f[67]+(5.031152949374527*alpha_diffusion[31]+5.1234753829798*alpha_diffusion[1])*f[39]+5.7282196186948*(alpha_diffusion[0]*f[16]+alpha_diffusion[1]*f[9])+5.1234753829798*alpha_diffusion[11]*f[16])*rdvxSq4nu+0.3872983346207416*f[49]*alpha_drag[171]+0.4330127018922193*(f[68]*alpha_drag[164]+f[32]*alpha_drag[161]+f[49]*alpha_drag[160])+1.299038105676658*alpha_drag[2]*f[68]+0.5809475019311124*alpha_drag[11]*f[54]+(1.322875655532295*alpha_drag[11]+1.479019945774904*alpha_drag[0])*f[40]+0.5809475019311124*f[25]*alpha_drag[31]+alpha_drag[1]*(1.479019945774904*f[26]+0.5916079783099616*f[25])+1.984313483298443*alpha_drag[2]*f[16]+0.5916079783099616*f[8]*alpha_drag[11]+0.6614378277661477*(alpha_drag[0]*f[8]+alpha_drag[1]*f[4]); 
  out[69] += (0.2581988897471612*f[21]+0.4330127018922193*f[3])*alpha_drag[191]+(0.2581988897471612*f[50]+0.3803194146278324*f[6])*alpha_drag[171]+0.4330127018922193*f[69]*alpha_drag[164]+0.3803194146278324*f[21]*alpha_drag[161]+0.4330127018922193*f[50]*alpha_drag[160]+(0.2581988897471612*f[25]+0.4330127018922193*f[4])*alpha_drag[111]+(0.2581988897471612*f[54]+0.3803194146278324*f[8])*alpha_drag[91]+0.4330127018922193*f[69]*alpha_drag[83]+0.3803194146278324*f[25]*alpha_drag[81]+0.4330127018922193*f[54]*alpha_drag[80]; 
  out[70] += 5.7282196186948*(alpha_diffusion[31]*f[76]+alpha_diffusion[11]*f[60]+alpha_diffusion[1]*f[35]+alpha_diffusion[0]*f[18])*rdvxSq4nu+0.4330127018922193*(f[70]*alpha_drag[164]+f[65]*alpha_drag[161]+f[51]*alpha_drag[160]+f[70]*alpha_drag[83]+f[68]*alpha_drag[81]+f[55]*alpha_drag[80])+1.299038105676658*alpha_drag[2]*f[70]+0.6614378277661477*alpha_drag[31]*f[69]+1.479019945774904*(alpha_drag[1]*f[61]+alpha_drag[0]*f[42])+0.6614378277661477*alpha_drag[11]*f[41]+1.984313483298443*alpha_drag[2]*f[18]+0.6614378277661477*(alpha_drag[1]*f[17]+alpha_drag[0]*f[10]); 
  out[71] += (5.031152949374527*alpha_diffusion[11]*f[69]+(5.031152949374527*alpha_diffusion[31]+5.1234753829798*alpha_diffusion[1])*f[41]+5.7282196186948*(alpha_diffusion[0]*f[17]+alpha_diffusion[1]*f[10])+5.1234753829798*alpha_diffusion[11]*f[17])*rdvySq4nu+0.3872983346207416*f[52]*alpha_drag[171]+0.4330127018922193*(f[71]*alpha_drag[164]+f[33]*alpha_drag[161]+f[52]*alpha_drag[160])+0.5809475019311124*f[25]*alpha_drag[111]+(0.5809475019311124*f[54]+1.322875655532295*f[43]+0.5916079783099616*f[8])*alpha_drag[91]+(1.299038105676658*f[71]+1.984313483298443*f[17])*alpha_drag[83]+(1.479019945774904*f[27]+0.5916079783099616*f[25]+0.6614378277661477*f[4])*alpha_drag[81]+(1.479019945774904*f[43]+0.6614378277661477*f[8])*alpha_drag[80]; 
  out[72] += 5.7282196186948*(alpha_diffusion[31]*f[76]+alpha_diffusion[11]*f[60]+alpha_diffusion[1]*f[35]+alpha_diffusion[0]*f[18])*rdvySq4nu+0.4330127018922193*(f[72]*alpha_drag[164]+f[66]*alpha_drag[161]+f[53]*alpha_drag[160])+0.6614378277661477*(f[67]*alpha_drag[111]+f[39]*alpha_drag[91])+(1.299038105676658*f[72]+1.984313483298443*f[18])*alpha_drag[83]+(1.479019945774904*f[62]+0.6614378277661477*f[16])*alpha_drag[81]+(1.479019945774904*f[44]+0.6614378277661477*f[9])*alpha_drag[80]+0.4330127018922193*(alpha_drag[2]*f[72]+alpha_drag[1]*f[71]+alpha_drag[0]*f[56]); 
  out[73] += (5.031152949374527*alpha_diffusion[11]*f[67]+(5.031152949374527*alpha_diffusion[31]+5.1234753829798*alpha_diffusion[1])*f[39]+5.7282196186948*(alpha_diffusion[0]*f[16]+alpha_diffusion[1]*f[9])+5.1234753829798*alpha_diffusion[11]*f[16])*rdvzSq4nu+0.5809475019311124*f[19]*alpha_drag[191]+(0.5809475019311124*f[48]+1.322875655532295*f[45]+0.5916079783099616*f[5])*alpha_drag[171]+(1.299038105676658*f[73]+1.984313483298443*f[16])*alpha_drag[164]+(1.479019945774904*f[29]+0.5916079783099616*f[19]+0.6614378277661477*f[2])*alpha_drag[161]+(1.479019945774904*f[45]+0.6614378277661477*f[5])*alpha_drag[160]+0.4330127018922193*alpha_drag[2]*f[73]+0.3872983346207416*alpha_drag[11]*f[57]+0.4330127018922193*(alpha_drag[0]*f[57]+alpha_drag[1]*f[34]); 
  out[74] += (5.031152949374527*alpha_diffusion[11]*f[69]+(5.031152949374527*alpha_diffusion[31]+5.1234753829798*alpha_diffusion[1])*f[41]+5.7282196186948*(alpha_diffusion[0]*f[17]+alpha_diffusion[1]*f[10])+5.1234753829798*alpha_diffusion[11]*f[17])*rdvzSq4nu+0.5809475019311124*f[21]*alpha_drag[191]+(0.5809475019311124*f[50]+1.322875655532295*f[46]+0.5916079783099616*f[6])*alpha_drag[171]+(1.299038105676658*f[74]+1.984313483298443*f[17])*alpha_drag[164]+(1.479019945774904*f[30]+0.5916079783099616*f[21]+0.6614378277661477*f[3])*alpha_drag[161]+(1.479019945774904*f[46]+0.6614378277661477*f[6])*alpha_drag[160]+0.3872983346207416*f[57]*alpha_drag[91]+0.4330127018922193*(f[74]*alpha_drag[83]+f[34]*alpha_drag[81]+f[57]*alpha_drag[80]); 
  out[75] += 5.7282196186948*(alpha_diffusion[31]*f[76]+alpha_diffusion[11]*f[60]+alpha_diffusion[1]*f[35]+alpha_diffusion[0]*f[18])*rdvzSq4nu+0.6614378277661477*(f[64]*alpha_drag[191]+f[36]*alpha_drag[171])+(1.299038105676658*f[75]+1.984313483298443*f[18])*alpha_drag[164]+(1.479019945774904*f[63]+0.6614378277661477*f[15])*alpha_drag[161]+(1.479019945774904*f[47]+0.6614378277661477*f[7])*alpha_drag[160]+0.4330127018922193*(f[75]*alpha_drag[83]+f[73]*alpha_drag[81]+f[58]*alpha_drag[80]+alpha_drag[2]*f[75]+alpha_drag[1]*f[74]+alpha_drag[0]*f[59]); 
  out[76] += (0.2581988897471612*f[36]+0.4330127018922193*f[7])*alpha_drag[191]+(0.2581988897471612*f[64]+0.3803194146278324*f[15])*alpha_drag[171]+0.4330127018922193*f[76]*alpha_drag[164]+0.3803194146278324*f[36]*alpha_drag[161]+0.4330127018922193*f[64]*alpha_drag[160]+(0.2581988897471612*f[39]+0.4330127018922193*f[9])*alpha_drag[111]+(0.2581988897471612*f[67]+0.3803194146278324*f[16])*alpha_drag[91]+0.4330127018922193*f[76]*alpha_drag[83]+0.3803194146278324*f[39]*alpha_drag[81]+0.4330127018922193*(f[67]*alpha_drag[80]+alpha_drag[2]*f[76])+(0.2581988897471612*alpha_drag[11]+0.4330127018922193*alpha_drag[0])*f[69]+(0.2581988897471612*alpha_drag[31]+0.3803194146278324*alpha_drag[1])*f[41]+0.4330127018922193*f[10]*alpha_drag[31]+0.3803194146278324*alpha_drag[11]*f[17]; 
  out[77] += (5.031152949374527*alpha_diffusion[11]*f[76]+(5.031152949374527*alpha_diffusion[31]+5.1234753829798*alpha_diffusion[1])*f[60]+5.7282196186948*(alpha_diffusion[0]*f[35]+alpha_diffusion[1]*f[18])+5.1234753829798*alpha_diffusion[11]*f[35])*rdvxSq4nu+0.3872983346207416*f[65]*alpha_drag[171]+0.4330127018922193*(f[77]*alpha_drag[164]+f[51]*alpha_drag[161]+f[65]*alpha_drag[160])+0.3872983346207416*f[68]*alpha_drag[91]+0.4330127018922193*(f[77]*alpha_drag[83]+f[55]*alpha_drag[81]+f[68]*alpha_drag[80])+1.299038105676658*alpha_drag[2]*f[77]+alpha_drag[11]*(0.5809475019311124*f[69]+1.322875655532295*f[61])+1.479019945774904*(alpha_drag[0]*f[61]+alpha_drag[1]*f[42])+(0.5809475019311124*alpha_drag[31]+0.5916079783099616*alpha_drag[1])*f[41]+1.984313483298443*alpha_drag[2]*f[35]+0.5916079783099616*alpha_drag[11]*f[17]+0.6614378277661477*(alpha_drag[0]*f[17]+alpha_drag[1]*f[10]); 
  out[78] += (5.031152949374527*alpha_diffusion[11]*f[76]+(5.031152949374527*alpha_diffusion[31]+5.1234753829798*alpha_diffusion[1])*f[60]+5.7282196186948*(alpha_diffusion[0]*f[35]+alpha_diffusion[1]*f[18])+5.1234753829798*alpha_diffusion[11]*f[35])*rdvySq4nu+0.3872983346207416*f[66]*alpha_drag[171]+0.4330127018922193*(f[78]*alpha_drag[164]+f[53]*alpha_drag[161]+f[66]*alpha_drag[160])+0.5809475019311124*f[39]*alpha_drag[111]+(0.5809475019311124*f[67]+1.322875655532295*f[62]+0.5916079783099616*f[16])*alpha_drag[91]+(1.299038105676658*f[78]+1.984313483298443*f[35])*alpha_drag[83]+(1.479019945774904*f[44]+0.5916079783099616*f[39]+0.6614378277661477*f[9])*alpha_drag[81]+(1.479019945774904*f[62]+0.6614378277661477*f[16])*alpha_drag[80]+0.4330127018922193*alpha_drag[2]*f[78]+0.3872983346207416*alpha_drag[11]*f[71]+0.4330127018922193*(alpha_drag[0]*f[71]+alpha_drag[1]*f[56]); 
  out[79] += (5.031152949374527*alpha_diffusion[11]*f[76]+(5.031152949374527*alpha_diffusion[31]+5.1234753829798*alpha_diffusion[1])*f[60]+5.7282196186948*(alpha_diffusion[0]*f[35]+alpha_diffusion[1]*f[18])+5.1234753829798*alpha_diffusion[11]*f[35])*rdvzSq4nu+0.5809475019311124*f[36]*alpha_drag[191]+(0.5809475019311124*f[64]+1.322875655532295*f[63]+0.5916079783099616*f[15])*alpha_drag[171]+(1.299038105676658*f[79]+1.984313483298443*f[35])*alpha_drag[164]+(1.479019945774904*f[47]+0.5916079783099616*f[36]+0.6614378277661477*f[7])*alpha_drag[161]+(1.479019945774904*f[63]+0.6614378277661477*f[15])*alpha_drag[160]+0.3872983346207416*f[73]*alpha_drag[91]+0.4330127018922193*(f[79]*alpha_drag[83]+f[58]*alpha_drag[81]+f[73]*alpha_drag[80]+alpha_drag[2]*f[79])+0.3872983346207416*alpha_drag[11]*f[74]+0.4330127018922193*(alpha_drag[0]*f[74]+alpha_drag[1]*f[59]); 

  return alpha_mid; 

} 
