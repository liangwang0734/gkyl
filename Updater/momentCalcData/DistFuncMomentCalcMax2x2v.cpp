#include <DistFuncMomentCalcModDecl.h> 
void MomentCalc2x2vMax_M0_P1(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[2]*dxv[3]/4; 
  out[0] += 2.0*f[0]*volFact; 
  out[1] += 2.0*f[1]*volFact; 
  out[2] += 2.0*f[2]*volFact; 
} 
void MomentCalc2x2vMax_M0_P2(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[2]*dxv[3]/4; 
  out[0] += 2.0*f[0]*volFact; 
  out[1] += 2.0*f[1]*volFact; 
  out[2] += 2.0*f[2]*volFact; 
  out[3] += 2.0*f[5]*volFact; 
  out[4] += 2.0*f[11]*volFact; 
  out[5] += 2.0*f[12]*volFact; 
} 
void MomentCalc2x2vMax_M1i_P1(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[2]*dxv[3]/4; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  out[0] += 2.0*f[0]*volFact*wx1+0.5773502691896258*f[3]*dv1*volFact; 
  out[1] += 2.0*f[1]*volFact*wx1; 
  out[2] += 2.0*f[2]*volFact*wx1; 
  out[3] += 2.0*f[0]*volFact*wx2+0.5773502691896258*f[4]*dv2*volFact; 
  out[4] += 2.0*f[1]*volFact*wx2; 
  out[5] += 2.0*f[2]*volFact*wx2; 
} 
void MomentCalc2x2vMax_M1i_P2(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[2]*dxv[3]/4; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  out[0] += 2.0*f[0]*volFact*wx1+0.5773502691896258*f[3]*dv1*volFact; 
  out[1] += 2.0*f[1]*volFact*wx1+0.5773502691896258*f[6]*dv1*volFact; 
  out[2] += 2.0*f[2]*volFact*wx1+0.5773502691896258*f[7]*dv1*volFact; 
  out[3] += 2.0*f[5]*volFact*wx1; 
  out[4] += 2.0*f[11]*volFact*wx1; 
  out[5] += 2.0*f[12]*volFact*wx1; 
  out[6] += 2.0*f[0]*volFact*wx2+0.5773502691896258*f[4]*dv2*volFact; 
  out[7] += 2.0*f[1]*volFact*wx2+0.5773502691896258*f[8]*dv2*volFact; 
  out[8] += 2.0*f[2]*volFact*wx2+0.5773502691896258*f[9]*dv2*volFact; 
  out[9] += 2.0*f[5]*volFact*wx2; 
  out[10] += 2.0*f[11]*volFact*wx2; 
  out[11] += 2.0*f[12]*volFact*wx2; 
} 
void MomentCalc2x2vMax_M2ij_P1(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[2]*dxv[3]/4; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  out[0] += 2.0*f[0]*volFact*wx1_sq+1.154700538379252*f[3]*dv1*volFact*wx1+0.1666666666666667*f[0]*dv1_sq*volFact; 
  out[1] += 2.0*f[1]*volFact*wx1_sq+0.1666666666666667*f[1]*dv1_sq*volFact; 
  out[2] += 2.0*f[2]*volFact*wx1_sq+0.1666666666666667*f[2]*dv1_sq*volFact; 
  out[3] += 2.0*f[0]*volFact*wx1*wx2+0.5773502691896258*f[3]*dv1*volFact*wx2+0.5773502691896258*f[4]*dv2*volFact*wx1; 
  out[4] += 2.0*f[1]*volFact*wx1*wx2; 
  out[5] += 2.0*f[2]*volFact*wx1*wx2; 
  out[6] += 2.0*f[0]*volFact*wx2_sq+1.154700538379252*f[4]*dv2*volFact*wx2+0.1666666666666667*f[0]*dv2_sq*volFact; 
  out[7] += 2.0*f[1]*volFact*wx2_sq+0.1666666666666667*f[1]*dv2_sq*volFact; 
  out[8] += 2.0*f[2]*volFact*wx2_sq+0.1666666666666667*f[2]*dv2_sq*volFact; 
} 
void MomentCalc2x2vMax_M2ij_P2(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[2]*dxv[3]/4; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  out[0] += 2.0*f[0]*volFact*wx1_sq+1.154700538379252*f[3]*dv1*volFact*wx1+0.149071198499986*f[13]*dv1_sq*volFact+0.1666666666666667*f[0]*dv1_sq*volFact; 
  out[1] += 2.0*f[1]*volFact*wx1_sq+1.154700538379252*f[6]*dv1*volFact*wx1+0.1666666666666667*f[1]*dv1_sq*volFact; 
  out[2] += 2.0*f[2]*volFact*wx1_sq+1.154700538379252*f[7]*dv1*volFact*wx1+0.1666666666666667*f[2]*dv1_sq*volFact; 
  out[3] += 2.0*f[5]*volFact*wx1_sq+0.1666666666666667*f[5]*dv1_sq*volFact; 
  out[4] += 2.0*f[11]*volFact*wx1_sq+0.1666666666666667*f[11]*dv1_sq*volFact; 
  out[5] += 2.0*f[12]*volFact*wx1_sq+0.1666666666666667*f[12]*dv1_sq*volFact; 
  out[6] += 2.0*f[0]*volFact*wx1*wx2+0.5773502691896258*f[3]*dv1*volFact*wx2+0.5773502691896258*f[4]*dv2*volFact*wx1+0.1666666666666667*f[10]*dv1*dv2*volFact; 
  out[7] += 2.0*f[1]*volFact*wx1*wx2+0.5773502691896258*f[6]*dv1*volFact*wx2+0.5773502691896258*f[8]*dv2*volFact*wx1; 
  out[8] += 2.0*f[2]*volFact*wx1*wx2+0.5773502691896258*f[7]*dv1*volFact*wx2+0.5773502691896258*f[9]*dv2*volFact*wx1; 
  out[9] += 2.0*f[5]*volFact*wx1*wx2; 
  out[10] += 2.0*f[11]*volFact*wx1*wx2; 
  out[11] += 2.0*f[12]*volFact*wx1*wx2; 
  out[12] += 2.0*f[0]*volFact*wx2_sq+1.154700538379252*f[4]*dv2*volFact*wx2+0.149071198499986*f[14]*dv2_sq*volFact+0.1666666666666667*f[0]*dv2_sq*volFact; 
  out[13] += 2.0*f[1]*volFact*wx2_sq+1.154700538379252*f[8]*dv2*volFact*wx2+0.1666666666666667*f[1]*dv2_sq*volFact; 
  out[14] += 2.0*f[2]*volFact*wx2_sq+1.154700538379252*f[9]*dv2*volFact*wx2+0.1666666666666667*f[2]*dv2_sq*volFact; 
  out[15] += 2.0*f[5]*volFact*wx2_sq+0.1666666666666667*f[5]*dv2_sq*volFact; 
  out[16] += 2.0*f[11]*volFact*wx2_sq+0.1666666666666667*f[11]*dv2_sq*volFact; 
  out[17] += 2.0*f[12]*volFact*wx2_sq+0.1666666666666667*f[12]*dv2_sq*volFact; 
} 
void MomentCalc2x2vMax_M2_P1(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[2]*dxv[3]/4; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  out[0] += 2.0*f[0]*volFact*wx2_sq+1.154700538379252*f[4]*dv2*volFact*wx2+2.0*f[0]*volFact*wx1_sq+1.154700538379252*f[3]*dv1*volFact*wx1+0.1666666666666667*f[0]*dv2_sq*volFact+0.1666666666666667*f[0]*dv1_sq*volFact; 
  out[1] += 2.0*f[1]*volFact*wx2_sq+2.0*f[1]*volFact*wx1_sq+0.1666666666666667*f[1]*dv2_sq*volFact+0.1666666666666667*f[1]*dv1_sq*volFact; 
  out[2] += 2.0*f[2]*volFact*wx2_sq+2.0*f[2]*volFact*wx1_sq+0.1666666666666667*f[2]*dv2_sq*volFact+0.1666666666666667*f[2]*dv1_sq*volFact; 
} 
void MomentCalc2x2vMax_M2_P2(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[2]*dxv[3]/4; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  out[0] += 2.0*f[0]*volFact*wx2_sq+1.154700538379252*f[4]*dv2*volFact*wx2+2.0*f[0]*volFact*wx1_sq+1.154700538379252*f[3]*dv1*volFact*wx1+0.149071198499986*f[14]*dv2_sq*volFact+0.1666666666666667*f[0]*dv2_sq*volFact+0.149071198499986*f[13]*dv1_sq*volFact+0.1666666666666667*f[0]*dv1_sq*volFact; 
  out[1] += 2.0*f[1]*volFact*wx2_sq+1.154700538379252*f[8]*dv2*volFact*wx2+2.0*f[1]*volFact*wx1_sq+1.154700538379252*f[6]*dv1*volFact*wx1+0.1666666666666667*f[1]*dv2_sq*volFact+0.1666666666666667*f[1]*dv1_sq*volFact; 
  out[2] += 2.0*f[2]*volFact*wx2_sq+1.154700538379252*f[9]*dv2*volFact*wx2+2.0*f[2]*volFact*wx1_sq+1.154700538379252*f[7]*dv1*volFact*wx1+0.1666666666666667*f[2]*dv2_sq*volFact+0.1666666666666667*f[2]*dv1_sq*volFact; 
  out[3] += 2.0*f[5]*volFact*wx2_sq+2.0*f[5]*volFact*wx1_sq+0.1666666666666667*f[5]*dv2_sq*volFact+0.1666666666666667*f[5]*dv1_sq*volFact; 
  out[4] += 2.0*f[11]*volFact*wx2_sq+2.0*f[11]*volFact*wx1_sq+0.1666666666666667*f[11]*dv2_sq*volFact+0.1666666666666667*f[11]*dv1_sq*volFact; 
  out[5] += 2.0*f[12]*volFact*wx2_sq+2.0*f[12]*volFact*wx1_sq+0.1666666666666667*f[12]*dv2_sq*volFact+0.1666666666666667*f[12]*dv1_sq*volFact; 
} 
