#include <VlasovModDecl.h> 
__host__ __device__ double VlasovSurfElcMag1x2vSer_VX_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *EM, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w:         Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // amax:      amax in global lax flux.
  // EM:        EM field.
  // fl/fr:     Distribution function in left/right cells 
  // outl/outr: output distribution function in left/right cells 
  // returns abs(amid) for use in determining amax in cfl and global lax flux 
  double dv10l = 2/dxvl[1]; 
  double dv10r = 2/dxvr[1]; 
  const double *E0 = &EM[0]; 
  const double dv1 = dxvr[1], wv1 = wr[1]; 
  const double dv2 = dxvr[2], wv2 = wr[2]; 
  const double *B2 = &EM[20]; 

  double Ghat[12]; 
  double favg[12]; 
  double alpha[12]; 

  favg[0] = (-1.870828693386971*fr[18])+1.870828693386971*fl[18]+1.58113883008419*fr[8]+1.58113883008419*fl[8]-1.224744871391589*fr[2]+1.224744871391589*fl[2]+0.7071067811865475*fr[0]+0.7071067811865475*fl[0]; 
  favg[1] = (-1.870828693386971*fr[24])+1.870828693386971*fl[24]+1.58113883008419*fr[12]+1.58113883008419*fl[12]-1.224744871391589*fr[4]+1.224744871391589*fl[4]+0.7071067811865475*fr[1]+0.7071067811865475*fl[1]; 
  favg[2] = (-1.870828693386971*fr[26])+1.870828693386971*fl[26]+1.58113883008419*fr[14]+1.58113883008419*fl[14]-1.224744871391589*fr[6]+1.224744871391589*fl[6]+0.7071067811865475*fr[3]+0.7071067811865475*fl[3]; 
  favg[3] = (-1.870828693386971*fr[30])+1.870828693386971*fl[30]+1.58113883008419*fr[21]+1.58113883008419*fl[21]-1.224744871391589*fr[10]+1.224744871391589*fl[10]+0.7071067811865475*fr[5]+0.7071067811865475*fl[5]; 
  favg[4] = (-1.224744871391589*fr[11])+1.224744871391589*fl[11]+0.7071067811865475*fr[7]+0.7071067811865475*fl[7]; 
  favg[5] = (-1.224744871391589*fr[16])+1.224744871391589*fl[16]+0.7071067811865475*fr[9]+0.7071067811865475*fl[9]; 
  favg[6] = (-1.224744871391589*fr[20])+1.224744871391589*fl[20]+0.7071067811865475*fr[13]+0.7071067811865475*fl[13]; 
  favg[7] = (-1.224744871391589*fr[22])+1.224744871391589*fl[22]+0.7071067811865475*fr[15]+0.7071067811865475*fl[15]; 
  favg[8] = (-1.224744871391589*fr[23])+1.224744871391589*fl[23]+0.7071067811865475*fr[17]+0.7071067811865475*fl[17]; 
  favg[9] = (-1.224744871391589*fr[28])+1.224744871391589*fl[28]+0.7071067811865475*fr[19]+0.7071067811865475*fl[19]; 
  favg[10] = (-1.224744871391589*fr[29])+1.224744871391589*fl[29]+0.7071067811865475*fr[25]+0.7071067811865475*fl[25]; 
  favg[11] = (-1.224744871391589*fr[31])+1.224744871391589*fl[31]+0.7071067811865475*fr[27]+0.7071067811865475*fl[27]; 

  alpha[0] = 1.414213562373095*(B2[0]*wv2+E0[0]); 
  alpha[1] = 1.414213562373095*(B2[1]*wv2+E0[1]); 
  alpha[2] = 0.408248290463863*B2[0]*dv2; 
  alpha[3] = 0.408248290463863*B2[1]*dv2; 
  alpha[4] = 1.414213562373095*(B2[2]*wv2+E0[2]); 
  alpha[6] = 0.408248290463863*B2[2]*dv2; 
  alpha[8] = 1.414213562373095*(B2[3]*wv2+E0[3]); 
  alpha[10] = 0.408248290463863*B2[3]*dv2; 

  const double amid = 0.5*alpha[0]-0.5590169943749475*alpha[4]; 

  Ghat[0] = 0.3535533905932737*(2.645751311064591*(fr[18]+fl[18])+2.23606797749979*(fl[8]-1.0*fr[8])+1.732050807568877*(fr[2]+fl[2])-1.0*fr[0]+fl[0])*amax+0.25*(alpha[10]*favg[10]+alpha[8]*favg[8]+alpha[6]*favg[6]+alpha[4]*favg[4]+alpha[3]*favg[3]+alpha[2]*favg[2]+alpha[1]*favg[1]+alpha[0]*favg[0]); 
  Ghat[1] = 0.1178511301977579*(1.732050807568877*(4.58257569495584*(fr[24]+fl[24])+3.872983346207417*(fl[12]-1.0*fr[12])+3.0*(fr[4]+fl[4]))+3.0*(fl[1]-1.0*fr[1]))*amax+0.002380952380952381*(92.22255689363637*(alpha[6]*favg[10]+favg[6]*alpha[10]+alpha[4]*favg[8]+favg[4]*alpha[8])+93.91485505499116*(alpha[3]*favg[6]+favg[3]*alpha[6]+alpha[1]*favg[4]+favg[1]*alpha[4])+105.0*(alpha[2]*favg[3]+favg[2]*alpha[3]+alpha[0]*favg[1]+favg[0]*alpha[1])); 
  Ghat[2] = 0.1178511301977579*(1.732050807568877*(4.58257569495584*(fr[26]+fl[26])+3.872983346207417*(fl[14]-1.0*fr[14])+3.0*(fr[6]+fl[6]))+3.0*(fl[3]-1.0*fr[3]))*amax+0.002380952380952381*(105.0*(alpha[8]*favg[10]+favg[8]*alpha[10])+93.91485505499116*alpha[3]*favg[7]+105.0*(alpha[4]*favg[6]+favg[4]*alpha[6])+93.91485505499116*alpha[2]*favg[5]+105.0*(alpha[1]*favg[3]+favg[1]*alpha[3]+alpha[0]*favg[2]+favg[0]*alpha[2])); 
  Ghat[3] = 0.3535533905932737*(2.645751311064591*(fr[30]+fl[30])+2.23606797749979*(fl[21]-1.0*fr[21])+1.732050807568877*(fr[10]+fl[10])-1.0*fr[5]+fl[5])*amax+0.002380952380952381*(92.22255689363638*(alpha[4]*favg[10]+favg[4]*alpha[10])+92.2225568936364*(alpha[6]*favg[8]+favg[6]*alpha[8])+(84.0*alpha[6]+93.91485505499116*alpha[2])*favg[7]+93.91485505499116*(alpha[1]*favg[6]+favg[1]*alpha[6]+alpha[3]*(favg[5]+favg[4])+favg[3]*alpha[4])+105.0*(alpha[0]*favg[3]+favg[0]*alpha[3]+alpha[1]*favg[2]+favg[1]*alpha[2])); 
  Ghat[4] = 0.07071067811865475*(8.660254037844387*(fr[11]+fl[11])+5.0*(fl[7]-1.0*fr[7]))*amax+0.002380952380952381*((62.60990336999411*alpha[10]+92.22255689363638*alpha[3])*favg[10]+92.22255689363638*favg[3]*alpha[10]+(62.60990336999411*alpha[8]+92.22255689363637*alpha[1])*favg[8]+92.22255689363637*favg[1]*alpha[8]+(67.0820393249937*alpha[6]+105.0*alpha[2])*favg[6]+105.0*favg[2]*alpha[6]+(67.0820393249937*alpha[4]+105.0*alpha[0])*favg[4]+105.0*favg[0]*alpha[4]+93.91485505499116*(alpha[3]*favg[3]+alpha[1]*favg[1])); 
  Ghat[5] = 0.07071067811865475*(8.660254037844387*(fr[16]+fl[16])+5.0*(fl[9]-1.0*fr[9]))*amax+0.002380952380952381*(92.22255689363638*alpha[3]*favg[11]+93.91485505499116*alpha[10]*favg[10]+92.22255689363637*alpha[2]*favg[9]+105.0*alpha[1]*favg[7]+93.91485505499116*alpha[6]*favg[6]+105.0*alpha[0]*favg[5]+93.91485505499116*(alpha[3]*favg[3]+alpha[2]*favg[2])); 
  Ghat[6] = 0.07071067811865475*(8.660254037844387*(fr[20]+fl[20])+5.0*(fl[13]-1.0*fr[13]))*amax+7.936507936507936e-4*((187.8297101099824*alpha[8]+276.6676706809091*alpha[1])*favg[10]+(187.8297101099824*favg[8]+247.4590875276153*favg[7]+276.6676706809091*favg[1])*alpha[10]+276.6676706809092*(alpha[3]*favg[8]+favg[3]*alpha[8])+252.0*alpha[3]*favg[7]+(201.2461179749811*alpha[4]+315.0*alpha[0])*favg[6]+(281.7445651649735*favg[5]+201.2461179749811*favg[4]+315.0*favg[0])*alpha[6]+314.9999999999999*(alpha[2]*favg[4]+favg[2]*alpha[4])+281.7445651649734*(alpha[1]*favg[3]+favg[1]*alpha[3])); 
  Ghat[7] = 0.07071067811865475*(8.660254037844387*(fr[22]+fl[22])+5.0*(fl[15]-1.0*fr[15]))*amax+0.002380952380952381*((82.48636250920512*alpha[6]+92.22255689363637*alpha[2])*favg[11]+82.48636250920512*(alpha[6]*favg[10]+favg[6]*alpha[10])+92.2225568936364*alpha[3]*favg[9]+(93.91485505499116*alpha[4]+105.0*alpha[0])*favg[7]+84.0*(alpha[3]*favg[6]+favg[3]*alpha[6])+105.0*alpha[1]*favg[5]+93.91485505499116*(alpha[2]*favg[3]+favg[2]*alpha[3])); 
  Ghat[8] = 0.05050762722761053*(12.12435565298214*(fr[23]+fl[23])+7.0*(fl[17]-1.0*fr[17]))*amax+7.936507936507936e-4*((187.8297101099824*alpha[6]+315.0*alpha[2])*favg[10]+(187.8297101099824*favg[6]+315.0*favg[2])*alpha[10]+(187.8297101099823*alpha[4]+315.0*alpha[0])*favg[8]+(187.8297101099823*favg[4]+315.0*favg[0])*alpha[8]+276.6676706809092*(alpha[3]*favg[6]+favg[3]*alpha[6])+276.6676706809091*(alpha[1]*favg[4]+favg[1]*alpha[4])); 
  Ghat[9] = 0.05050762722761053*(12.12435565298214*(fr[28]+fl[28])+7.0*(fl[19]-1.0*fr[19]))*amax+0.002380952380952381*(105.0*alpha[1]*favg[11]+105.0*alpha[0]*favg[9]+92.2225568936364*alpha[3]*favg[7]+92.22255689363637*alpha[2]*favg[5]); 
  Ghat[10] = 0.05050762722761053*(12.12435565298214*(fr[29]+fl[29])+7.0*(fl[25]-1.0*fr[25]))*amax+7.936507936507936e-4*((187.8297101099823*alpha[4]+315.0*alpha[0])*favg[10]+(281.7445651649735*favg[5]+187.8297101099823*favg[4]+315.0*favg[0])*alpha[10]+(187.8297101099824*alpha[6]+315.0*alpha[2])*favg[8]+(187.8297101099824*favg[6]+315.0*favg[2])*alpha[8]+247.4590875276153*alpha[6]*favg[7]+276.6676706809091*(alpha[1]*favg[6]+favg[1]*alpha[6])+276.6676706809092*(alpha[3]*favg[4]+favg[3]*alpha[4])); 
  Ghat[11] = 0.05050762722761053*(12.12435565298214*(fr[31]+fl[31])+7.0*(fl[27]-1.0*fr[27]))*amax+0.002380952380952381*((93.91485505499116*alpha[4]+105.0*alpha[0])*favg[11]+105.0*alpha[1]*favg[9]+(82.48636250920512*alpha[6]+92.22255689363637*alpha[2])*favg[7]+92.22255689363638*alpha[3]*favg[5]); 

  outr[0] += 0.7071067811865475*Ghat[0]*dv10r; 
  outr[1] += 0.7071067811865475*Ghat[1]*dv10r; 
  outr[2] += -1.224744871391589*Ghat[0]*dv10r; 
  outr[3] += 0.7071067811865475*Ghat[2]*dv10r; 
  outr[4] += -1.224744871391589*Ghat[1]*dv10r; 
  outr[5] += 0.7071067811865475*Ghat[3]*dv10r; 
  outr[6] += -1.224744871391589*Ghat[2]*dv10r; 
  outr[7] += 0.7071067811865475*Ghat[4]*dv10r; 
  outr[8] += 1.58113883008419*Ghat[0]*dv10r; 
  outr[9] += 0.7071067811865475*Ghat[5]*dv10r; 
  outr[10] += -1.224744871391589*Ghat[3]*dv10r; 
  outr[11] += -1.224744871391589*Ghat[4]*dv10r; 
  outr[12] += 1.58113883008419*Ghat[1]*dv10r; 
  outr[13] += 0.7071067811865475*Ghat[6]*dv10r; 
  outr[14] += 1.58113883008419*Ghat[2]*dv10r; 
  outr[15] += 0.7071067811865475*Ghat[7]*dv10r; 
  outr[16] += -1.224744871391589*Ghat[5]*dv10r; 
  outr[17] += 0.7071067811865475*Ghat[8]*dv10r; 
  outr[18] += -1.870828693386971*Ghat[0]*dv10r; 
  outr[19] += 0.7071067811865475*Ghat[9]*dv10r; 
  outr[20] += -1.224744871391589*Ghat[6]*dv10r; 
  outr[21] += 1.58113883008419*Ghat[3]*dv10r; 
  outr[22] += -1.224744871391589*Ghat[7]*dv10r; 
  outr[23] += -1.224744871391589*Ghat[8]*dv10r; 
  outr[24] += -1.870828693386971*Ghat[1]*dv10r; 
  outr[25] += 0.7071067811865475*Ghat[10]*dv10r; 
  outr[26] += -1.870828693386971*Ghat[2]*dv10r; 
  outr[27] += 0.7071067811865475*Ghat[11]*dv10r; 
  outr[28] += -1.224744871391589*Ghat[9]*dv10r; 
  outr[29] += -1.224744871391589*Ghat[10]*dv10r; 
  outr[30] += -1.870828693386971*Ghat[3]*dv10r; 
  outr[31] += -1.224744871391589*Ghat[11]*dv10r; 

  outl[0] += -0.7071067811865475*Ghat[0]*dv10l; 
  outl[1] += -0.7071067811865475*Ghat[1]*dv10l; 
  outl[2] += -1.224744871391589*Ghat[0]*dv10l; 
  outl[3] += -0.7071067811865475*Ghat[2]*dv10l; 
  outl[4] += -1.224744871391589*Ghat[1]*dv10l; 
  outl[5] += -0.7071067811865475*Ghat[3]*dv10l; 
  outl[6] += -1.224744871391589*Ghat[2]*dv10l; 
  outl[7] += -0.7071067811865475*Ghat[4]*dv10l; 
  outl[8] += -1.58113883008419*Ghat[0]*dv10l; 
  outl[9] += -0.7071067811865475*Ghat[5]*dv10l; 
  outl[10] += -1.224744871391589*Ghat[3]*dv10l; 
  outl[11] += -1.224744871391589*Ghat[4]*dv10l; 
  outl[12] += -1.58113883008419*Ghat[1]*dv10l; 
  outl[13] += -0.7071067811865475*Ghat[6]*dv10l; 
  outl[14] += -1.58113883008419*Ghat[2]*dv10l; 
  outl[15] += -0.7071067811865475*Ghat[7]*dv10l; 
  outl[16] += -1.224744871391589*Ghat[5]*dv10l; 
  outl[17] += -0.7071067811865475*Ghat[8]*dv10l; 
  outl[18] += -1.870828693386971*Ghat[0]*dv10l; 
  outl[19] += -0.7071067811865475*Ghat[9]*dv10l; 
  outl[20] += -1.224744871391589*Ghat[6]*dv10l; 
  outl[21] += -1.58113883008419*Ghat[3]*dv10l; 
  outl[22] += -1.224744871391589*Ghat[7]*dv10l; 
  outl[23] += -1.224744871391589*Ghat[8]*dv10l; 
  outl[24] += -1.870828693386971*Ghat[1]*dv10l; 
  outl[25] += -0.7071067811865475*Ghat[10]*dv10l; 
  outl[26] += -1.870828693386971*Ghat[2]*dv10l; 
  outl[27] += -0.7071067811865475*Ghat[11]*dv10l; 
  outl[28] += -1.224744871391589*Ghat[9]*dv10l; 
  outl[29] += -1.224744871391589*Ghat[10]*dv10l; 
  outl[30] += -1.870828693386971*Ghat[3]*dv10l; 
  outl[31] += -1.224744871391589*Ghat[11]*dv10l; 

  return std::abs(amid); 
} 
__host__ __device__ double VlasovSurfElcMag1x2vSer_VY_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *EM, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w:         Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // amax:      amax in global lax flux.
  // EM:        EM field.
  // fl/fr:     Distribution function in left/right cells 
  // outl/outr: output distribution function in left/right cells 
  // returns abs(amid) for use in determining amax in cfl and global lax flux 
  double dv11l = 2/dxvl[2]; 
  double dv11r = 2/dxvr[2]; 
  const double *E1 = &EM[4]; 
  const double dv1 = dxvr[1], wv1 = wr[1]; 
  const double dv2 = dxvr[2], wv2 = wr[2]; 
  const double *B2 = &EM[20]; 

  double Ghat[12]; 
  double favg[12]; 
  double alpha[12]; 

  favg[0] = (-1.870828693386971*fr[19])+1.870828693386971*fl[19]+1.58113883008419*fr[9]+1.58113883008419*fl[9]-1.224744871391589*fr[3]+1.224744871391589*fl[3]+0.7071067811865475*fr[0]+0.7071067811865475*fl[0]; 
  favg[1] = (-1.870828693386971*fr[27])+1.870828693386971*fl[27]+1.58113883008419*fr[15]+1.58113883008419*fl[15]-1.224744871391589*fr[5]+1.224744871391589*fl[5]+0.7071067811865475*fr[1]+0.7071067811865475*fl[1]; 
  favg[2] = (-1.870828693386971*fr[28])+1.870828693386971*fl[28]+1.58113883008419*fr[16]+1.58113883008419*fl[16]-1.224744871391589*fr[6]+1.224744871391589*fl[6]+0.7071067811865475*fr[2]+0.7071067811865475*fl[2]; 
  favg[3] = (-1.870828693386971*fr[31])+1.870828693386971*fl[31]+1.58113883008419*fr[22]+1.58113883008419*fl[22]-1.224744871391589*fr[10]+1.224744871391589*fl[10]+0.7071067811865475*fr[4]+0.7071067811865475*fl[4]; 
  favg[4] = (-1.224744871391589*fr[13])+1.224744871391589*fl[13]+0.7071067811865475*fr[7]+0.7071067811865475*fl[7]; 
  favg[5] = (-1.224744871391589*fr[14])+1.224744871391589*fl[14]+0.7071067811865475*fr[8]+0.7071067811865475*fl[8]; 
  favg[6] = (-1.224744871391589*fr[20])+1.224744871391589*fl[20]+0.7071067811865475*fr[11]+0.7071067811865475*fl[11]; 
  favg[7] = (-1.224744871391589*fr[21])+1.224744871391589*fl[21]+0.7071067811865475*fr[12]+0.7071067811865475*fl[12]; 
  favg[8] = (-1.224744871391589*fr[25])+1.224744871391589*fl[25]+0.7071067811865475*fr[17]+0.7071067811865475*fl[17]; 
  favg[9] = (-1.224744871391589*fr[26])+1.224744871391589*fl[26]+0.7071067811865475*fr[18]+0.7071067811865475*fl[18]; 
  favg[10] = (-1.224744871391589*fr[29])+1.224744871391589*fl[29]+0.7071067811865475*fr[23]+0.7071067811865475*fl[23]; 
  favg[11] = (-1.224744871391589*fr[30])+1.224744871391589*fl[30]+0.7071067811865475*fr[24]+0.7071067811865475*fl[24]; 

  alpha[0] = 1.414213562373095*E1[0]-1.414213562373095*B2[0]*wv1; 
  alpha[1] = 1.414213562373095*E1[1]-1.414213562373095*B2[1]*wv1; 
  alpha[2] = -0.408248290463863*B2[0]*dv1; 
  alpha[3] = -0.408248290463863*B2[1]*dv1; 
  alpha[4] = 1.414213562373095*E1[2]-1.414213562373095*B2[2]*wv1; 
  alpha[6] = -0.408248290463863*B2[2]*dv1; 
  alpha[8] = 1.414213562373095*E1[3]-1.414213562373095*B2[3]*wv1; 
  alpha[10] = -0.408248290463863*B2[3]*dv1; 

  const double amid = 0.5*alpha[0]-0.5590169943749475*alpha[4]; 

  Ghat[0] = 0.3535533905932737*(2.645751311064591*(fr[19]+fl[19])+2.23606797749979*(fl[9]-1.0*fr[9])+1.732050807568877*(fr[3]+fl[3])-1.0*fr[0]+fl[0])*amax+0.25*(alpha[10]*favg[10]+alpha[8]*favg[8]+alpha[6]*favg[6]+alpha[4]*favg[4]+alpha[3]*favg[3]+alpha[2]*favg[2]+alpha[1]*favg[1]+alpha[0]*favg[0]); 
  Ghat[1] = 0.1178511301977579*(1.732050807568877*(4.58257569495584*(fr[27]+fl[27])+3.872983346207417*(fl[15]-1.0*fr[15])+3.0*(fr[5]+fl[5]))+3.0*(fl[1]-1.0*fr[1]))*amax+0.002380952380952381*(92.22255689363637*(alpha[6]*favg[10]+favg[6]*alpha[10]+alpha[4]*favg[8]+favg[4]*alpha[8])+93.91485505499116*(alpha[3]*favg[6]+favg[3]*alpha[6]+alpha[1]*favg[4]+favg[1]*alpha[4])+105.0*(alpha[2]*favg[3]+favg[2]*alpha[3]+alpha[0]*favg[1]+favg[0]*alpha[1])); 
  Ghat[2] = 0.1178511301977579*(1.732050807568877*(4.58257569495584*(fr[28]+fl[28])+3.872983346207417*(fl[16]-1.0*fr[16])+3.0*(fr[6]+fl[6]))+3.0*(fl[2]-1.0*fr[2]))*amax+0.002380952380952381*(105.0*(alpha[8]*favg[10]+favg[8]*alpha[10])+93.91485505499116*alpha[3]*favg[7]+105.0*(alpha[4]*favg[6]+favg[4]*alpha[6])+93.91485505499116*alpha[2]*favg[5]+105.0*(alpha[1]*favg[3]+favg[1]*alpha[3]+alpha[0]*favg[2]+favg[0]*alpha[2])); 
  Ghat[3] = 0.3535533905932737*(2.645751311064591*(fr[31]+fl[31])+2.23606797749979*(fl[22]-1.0*fr[22])+1.732050807568877*(fr[10]+fl[10])-1.0*fr[4]+fl[4])*amax+0.002380952380952381*(92.22255689363638*(alpha[4]*favg[10]+favg[4]*alpha[10])+92.2225568936364*(alpha[6]*favg[8]+favg[6]*alpha[8])+(84.0*alpha[6]+93.91485505499116*alpha[2])*favg[7]+93.91485505499116*(alpha[1]*favg[6]+favg[1]*alpha[6]+alpha[3]*(favg[5]+favg[4])+favg[3]*alpha[4])+105.0*(alpha[0]*favg[3]+favg[0]*alpha[3]+alpha[1]*favg[2]+favg[1]*alpha[2])); 
  Ghat[4] = 0.07071067811865475*(8.660254037844387*(fr[13]+fl[13])+5.0*(fl[7]-1.0*fr[7]))*amax+0.002380952380952381*((62.60990336999411*alpha[10]+92.22255689363638*alpha[3])*favg[10]+92.22255689363638*favg[3]*alpha[10]+(62.60990336999411*alpha[8]+92.22255689363637*alpha[1])*favg[8]+92.22255689363637*favg[1]*alpha[8]+(67.0820393249937*alpha[6]+105.0*alpha[2])*favg[6]+105.0*favg[2]*alpha[6]+(67.0820393249937*alpha[4]+105.0*alpha[0])*favg[4]+105.0*favg[0]*alpha[4]+93.91485505499116*(alpha[3]*favg[3]+alpha[1]*favg[1])); 
  Ghat[5] = 0.07071067811865475*(8.660254037844387*(fr[14]+fl[14])+5.0*(fl[8]-1.0*fr[8]))*amax+0.002380952380952381*(92.22255689363638*alpha[3]*favg[11]+93.91485505499116*alpha[10]*favg[10]+92.22255689363637*alpha[2]*favg[9]+105.0*alpha[1]*favg[7]+93.91485505499116*alpha[6]*favg[6]+105.0*alpha[0]*favg[5]+93.91485505499116*(alpha[3]*favg[3]+alpha[2]*favg[2])); 
  Ghat[6] = 0.07071067811865475*(8.660254037844387*(fr[20]+fl[20])+5.0*(fl[11]-1.0*fr[11]))*amax+7.936507936507936e-4*((187.8297101099824*alpha[8]+276.6676706809091*alpha[1])*favg[10]+(187.8297101099824*favg[8]+247.4590875276153*favg[7]+276.6676706809091*favg[1])*alpha[10]+276.6676706809092*(alpha[3]*favg[8]+favg[3]*alpha[8])+252.0*alpha[3]*favg[7]+(201.2461179749811*alpha[4]+315.0*alpha[0])*favg[6]+(281.7445651649735*favg[5]+201.2461179749811*favg[4]+315.0*favg[0])*alpha[6]+314.9999999999999*(alpha[2]*favg[4]+favg[2]*alpha[4])+281.7445651649734*(alpha[1]*favg[3]+favg[1]*alpha[3])); 
  Ghat[7] = 0.07071067811865475*(8.660254037844387*(fr[21]+fl[21])+5.0*(fl[12]-1.0*fr[12]))*amax+0.002380952380952381*((82.48636250920512*alpha[6]+92.22255689363637*alpha[2])*favg[11]+82.48636250920512*(alpha[6]*favg[10]+favg[6]*alpha[10])+92.2225568936364*alpha[3]*favg[9]+(93.91485505499116*alpha[4]+105.0*alpha[0])*favg[7]+84.0*(alpha[3]*favg[6]+favg[3]*alpha[6])+105.0*alpha[1]*favg[5]+93.91485505499116*(alpha[2]*favg[3]+favg[2]*alpha[3])); 
  Ghat[8] = 0.05050762722761053*(12.12435565298214*(fr[25]+fl[25])+7.0*(fl[17]-1.0*fr[17]))*amax+7.936507936507936e-4*((187.8297101099824*alpha[6]+315.0*alpha[2])*favg[10]+(187.8297101099824*favg[6]+315.0*favg[2])*alpha[10]+(187.8297101099823*alpha[4]+315.0*alpha[0])*favg[8]+(187.8297101099823*favg[4]+315.0*favg[0])*alpha[8]+276.6676706809092*(alpha[3]*favg[6]+favg[3]*alpha[6])+276.6676706809091*(alpha[1]*favg[4]+favg[1]*alpha[4])); 
  Ghat[9] = 0.05050762722761053*(12.12435565298214*(fr[26]+fl[26])+7.0*(fl[18]-1.0*fr[18]))*amax+0.002380952380952381*(105.0*alpha[1]*favg[11]+105.0*alpha[0]*favg[9]+92.2225568936364*alpha[3]*favg[7]+92.22255689363637*alpha[2]*favg[5]); 
  Ghat[10] = 0.05050762722761053*(12.12435565298214*(fr[29]+fl[29])+7.0*(fl[23]-1.0*fr[23]))*amax+7.936507936507936e-4*((187.8297101099823*alpha[4]+315.0*alpha[0])*favg[10]+(281.7445651649735*favg[5]+187.8297101099823*favg[4]+315.0*favg[0])*alpha[10]+(187.8297101099824*alpha[6]+315.0*alpha[2])*favg[8]+(187.8297101099824*favg[6]+315.0*favg[2])*alpha[8]+247.4590875276153*alpha[6]*favg[7]+276.6676706809091*(alpha[1]*favg[6]+favg[1]*alpha[6])+276.6676706809092*(alpha[3]*favg[4]+favg[3]*alpha[4])); 
  Ghat[11] = 0.05050762722761053*(12.12435565298214*(fr[30]+fl[30])+7.0*(fl[24]-1.0*fr[24]))*amax+0.002380952380952381*((93.91485505499116*alpha[4]+105.0*alpha[0])*favg[11]+105.0*alpha[1]*favg[9]+(82.48636250920512*alpha[6]+92.22255689363637*alpha[2])*favg[7]+92.22255689363638*alpha[3]*favg[5]); 

  outr[0] += 0.7071067811865475*Ghat[0]*dv11r; 
  outr[1] += 0.7071067811865475*Ghat[1]*dv11r; 
  outr[2] += 0.7071067811865475*Ghat[2]*dv11r; 
  outr[3] += -1.224744871391589*Ghat[0]*dv11r; 
  outr[4] += 0.7071067811865475*Ghat[3]*dv11r; 
  outr[5] += -1.224744871391589*Ghat[1]*dv11r; 
  outr[6] += -1.224744871391589*Ghat[2]*dv11r; 
  outr[7] += 0.7071067811865475*Ghat[4]*dv11r; 
  outr[8] += 0.7071067811865475*Ghat[5]*dv11r; 
  outr[9] += 1.58113883008419*Ghat[0]*dv11r; 
  outr[10] += -1.224744871391589*Ghat[3]*dv11r; 
  outr[11] += 0.7071067811865475*Ghat[6]*dv11r; 
  outr[12] += 0.7071067811865475*Ghat[7]*dv11r; 
  outr[13] += -1.224744871391589*Ghat[4]*dv11r; 
  outr[14] += -1.224744871391589*Ghat[5]*dv11r; 
  outr[15] += 1.58113883008419*Ghat[1]*dv11r; 
  outr[16] += 1.58113883008419*Ghat[2]*dv11r; 
  outr[17] += 0.7071067811865475*Ghat[8]*dv11r; 
  outr[18] += 0.7071067811865475*Ghat[9]*dv11r; 
  outr[19] += -1.870828693386971*Ghat[0]*dv11r; 
  outr[20] += -1.224744871391589*Ghat[6]*dv11r; 
  outr[21] += -1.224744871391589*Ghat[7]*dv11r; 
  outr[22] += 1.58113883008419*Ghat[3]*dv11r; 
  outr[23] += 0.7071067811865475*Ghat[10]*dv11r; 
  outr[24] += 0.7071067811865475*Ghat[11]*dv11r; 
  outr[25] += -1.224744871391589*Ghat[8]*dv11r; 
  outr[26] += -1.224744871391589*Ghat[9]*dv11r; 
  outr[27] += -1.870828693386971*Ghat[1]*dv11r; 
  outr[28] += -1.870828693386971*Ghat[2]*dv11r; 
  outr[29] += -1.224744871391589*Ghat[10]*dv11r; 
  outr[30] += -1.224744871391589*Ghat[11]*dv11r; 
  outr[31] += -1.870828693386971*Ghat[3]*dv11r; 

  outl[0] += -0.7071067811865475*Ghat[0]*dv11l; 
  outl[1] += -0.7071067811865475*Ghat[1]*dv11l; 
  outl[2] += -0.7071067811865475*Ghat[2]*dv11l; 
  outl[3] += -1.224744871391589*Ghat[0]*dv11l; 
  outl[4] += -0.7071067811865475*Ghat[3]*dv11l; 
  outl[5] += -1.224744871391589*Ghat[1]*dv11l; 
  outl[6] += -1.224744871391589*Ghat[2]*dv11l; 
  outl[7] += -0.7071067811865475*Ghat[4]*dv11l; 
  outl[8] += -0.7071067811865475*Ghat[5]*dv11l; 
  outl[9] += -1.58113883008419*Ghat[0]*dv11l; 
  outl[10] += -1.224744871391589*Ghat[3]*dv11l; 
  outl[11] += -0.7071067811865475*Ghat[6]*dv11l; 
  outl[12] += -0.7071067811865475*Ghat[7]*dv11l; 
  outl[13] += -1.224744871391589*Ghat[4]*dv11l; 
  outl[14] += -1.224744871391589*Ghat[5]*dv11l; 
  outl[15] += -1.58113883008419*Ghat[1]*dv11l; 
  outl[16] += -1.58113883008419*Ghat[2]*dv11l; 
  outl[17] += -0.7071067811865475*Ghat[8]*dv11l; 
  outl[18] += -0.7071067811865475*Ghat[9]*dv11l; 
  outl[19] += -1.870828693386971*Ghat[0]*dv11l; 
  outl[20] += -1.224744871391589*Ghat[6]*dv11l; 
  outl[21] += -1.224744871391589*Ghat[7]*dv11l; 
  outl[22] += -1.58113883008419*Ghat[3]*dv11l; 
  outl[23] += -0.7071067811865475*Ghat[10]*dv11l; 
  outl[24] += -0.7071067811865475*Ghat[11]*dv11l; 
  outl[25] += -1.224744871391589*Ghat[8]*dv11l; 
  outl[26] += -1.224744871391589*Ghat[9]*dv11l; 
  outl[27] += -1.870828693386971*Ghat[1]*dv11l; 
  outl[28] += -1.870828693386971*Ghat[2]*dv11l; 
  outl[29] += -1.224744871391589*Ghat[10]*dv11l; 
  outl[30] += -1.224744871391589*Ghat[11]*dv11l; 
  outl[31] += -1.870828693386971*Ghat[3]*dv11l; 

  return std::abs(amid); 
} 
