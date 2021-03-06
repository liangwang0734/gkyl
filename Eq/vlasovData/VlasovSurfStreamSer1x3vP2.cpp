#include <VlasovModDecl.h> 
__host__ __device__ void VlasovSurfStream1x3vSer_X_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. fl/fr: Distribution function in left/right cells 
// outl/outr: Incremented distribution function in left/right cells 
  double dvxl = dxvl[1]; 
  double wxl = wl[1]; 

  double dvxr = dxvr[1]; 
  double wxr = wr[1]; 

  double dxl = 1.0/dxvl[0]; 
  double dxr = 1.0/dxvr[0]; 

  double incr[48]; 

  if (wxr>0) { 
  incr[0] = (2.23606797749979*fl[11]+1.732050807568877*fl[1]+fl[0])*wxl+(0.6454972243679028*fl[19]+0.5*fl[5]+0.2886751345948129*fl[2])*dvxl; 
  incr[1] = ((-3.872983346207417*fl[11])-3.0*fl[1]-1.732050807568877*fl[0])*wxl+((-1.118033988749895*fl[19])-0.8660254037844386*fl[5]-0.5*fl[2])*dvxl; 
  incr[2] = (2.23606797749979*fl[19]+1.732050807568877*fl[5]+fl[2])*wxl+(0.447213595499958*fl[20]+0.2581988897471612*fl[12]+0.6454972243679029*fl[11]+0.5*fl[1]+0.2886751345948129*fl[0])*dvxl; 
  incr[3] = (2.23606797749979*fl[21]+1.732050807568877*fl[6]+fl[3])*wxl+(0.6454972243679029*fl[32]+0.5*fl[15]+0.2886751345948129*fl[7])*dvxl; 
  incr[4] = (2.23606797749979*fl[25]+1.732050807568877*fl[8]+fl[4])*wxl+(0.6454972243679029*fl[35]+0.5*fl[16]+0.2886751345948129*fl[9])*dvxl; 
  incr[5] = ((-3.872983346207417*fl[19])-3.0*fl[5]-1.732050807568877*fl[2])*wxl+((-0.7745966692414834*fl[20])-0.4472135954999579*fl[12]-1.118033988749895*fl[11]-0.8660254037844386*fl[1]-0.5*fl[0])*dvxl; 
  incr[6] = ((-3.872983346207417*fl[21])-3.0*fl[6]-1.732050807568877*fl[3])*wxl+((-1.118033988749895*fl[32])-0.8660254037844386*fl[15]-0.5*fl[7])*dvxl; 
  incr[7] = (2.23606797749979*fl[32]+1.732050807568877*fl[15]+fl[7])*wxl+(0.4472135954999579*fl[33]+0.2581988897471611*fl[22]+0.6454972243679028*fl[21]+0.5*fl[6]+0.2886751345948129*fl[3])*dvxl; 
  incr[8] = ((-3.872983346207417*fl[25])-3.0*fl[8]-1.732050807568877*fl[4])*wxl+((-1.118033988749895*fl[35])-0.8660254037844386*fl[16]-0.5*fl[9])*dvxl; 
  incr[9] = (2.23606797749979*fl[35]+1.732050807568877*fl[16]+fl[9])*wxl+(0.4472135954999579*fl[36]+0.2581988897471611*fl[26]+0.6454972243679028*fl[25]+0.5*fl[8]+0.2886751345948129*fl[4])*dvxl; 
  incr[10] = (2.23606797749979*fl[37]+1.732050807568877*fl[17]+fl[10])*wxl+(0.6454972243679028*fl[44]+0.5*fl[31]+0.2886751345948129*fl[18])*dvxl; 
  incr[11] = (5.0*fl[11]+3.872983346207417*fl[1]+2.23606797749979*fl[0])*wxl+(1.443375672974065*fl[19]+1.118033988749895*fl[5]+0.6454972243679029*fl[2])*dvxl; 
  incr[12] = (1.732050807568877*fl[20]+fl[12])*wxl+(0.5773502691896257*fl[19]+0.4472135954999579*fl[5]+0.2581988897471612*fl[2])*dvxl; 
  incr[13] = (1.732050807568877*fl[23]+fl[13])*wxl+(0.5*fl[34]+0.2886751345948129*fl[24])*dvxl; 
  incr[14] = (1.732050807568877*fl[28]+fl[14])*wxl+(0.5*fl[41]+0.2886751345948129*fl[29])*dvxl; 
  incr[15] = ((-3.872983346207417*fl[32])-3.0*fl[15]-1.732050807568877*fl[7])*wxl+((-0.7745966692414833*fl[33])-0.447213595499958*fl[22]-1.118033988749895*fl[21]-0.8660254037844386*fl[6]-0.5*fl[3])*dvxl; 
  incr[16] = ((-3.872983346207417*fl[35])-3.0*fl[16]-1.732050807568877*fl[9])*wxl+((-0.7745966692414833*fl[36])-0.447213595499958*fl[26]-1.118033988749895*fl[25]-0.8660254037844386*fl[8]-0.5*fl[4])*dvxl; 
  incr[17] = ((-3.872983346207417*fl[37])-3.0*fl[17]-1.732050807568877*fl[10])*wxl+((-1.118033988749895*fl[44])-0.8660254037844386*fl[31]-0.5*fl[18])*dvxl; 
  incr[18] = (2.23606797749979*fl[44]+1.732050807568877*fl[31]+fl[18])*wxl+(0.447213595499958*fl[45]+0.2581988897471612*fl[38]+0.6454972243679029*fl[37]+0.5*fl[17]+0.2886751345948129*fl[10])*dvxl; 
  incr[19] = (5.0*fl[19]+3.872983346207417*fl[5]+2.23606797749979*fl[2])*wxl+(fl[20]+0.5773502691896257*fl[12]+1.443375672974065*fl[11]+1.118033988749895*fl[1]+0.6454972243679028*fl[0])*dvxl; 
  incr[20] = ((-3.0*fl[20])-1.732050807568877*fl[12])*wxl+((-1.0*fl[19])-0.7745966692414834*fl[5]-0.447213595499958*fl[2])*dvxl; 
  incr[21] = (5.0*fl[21]+3.872983346207417*fl[6]+2.23606797749979*fl[3])*wxl+(1.443375672974065*fl[32]+1.118033988749895*fl[15]+0.6454972243679028*fl[7])*dvxl; 
  incr[22] = (1.732050807568877*fl[33]+fl[22])*wxl+(0.5773502691896257*fl[32]+0.447213595499958*fl[15]+0.2581988897471611*fl[7])*dvxl; 
  incr[23] = ((-3.0*fl[23])-1.732050807568877*fl[13])*wxl+((-0.8660254037844387*fl[34])-0.5*fl[24])*dvxl; 
  incr[24] = (1.732050807568877*fl[34]+fl[24])*wxl+(0.5*fl[23]+0.2886751345948129*fl[13])*dvxl; 
  incr[25] = (5.0*fl[25]+3.872983346207417*fl[8]+2.23606797749979*fl[4])*wxl+(1.443375672974065*fl[35]+1.118033988749895*fl[16]+0.6454972243679028*fl[9])*dvxl; 
  incr[26] = (1.732050807568877*fl[36]+fl[26])*wxl+(0.5773502691896257*fl[35]+0.447213595499958*fl[16]+0.2581988897471611*fl[9])*dvxl; 
  incr[27] = (1.732050807568877*fl[39]+fl[27])*wxl+(0.5*fl[46]+0.2886751345948129*fl[40])*dvxl; 
  incr[28] = ((-3.0*fl[28])-1.732050807568877*fl[14])*wxl+((-0.8660254037844387*fl[41])-0.5*fl[29])*dvxl; 
  incr[29] = (1.732050807568877*fl[41]+fl[29])*wxl+(0.5*fl[28]+0.2886751345948129*fl[14])*dvxl; 
  incr[30] = (1.732050807568877*fl[42]+fl[30])*wxl+(0.5*fl[47]+0.2886751345948129*fl[43])*dvxl; 
  incr[31] = ((-3.872983346207417*fl[44])-3.0*fl[31]-1.732050807568877*fl[18])*wxl+((-0.7745966692414834*fl[45])-0.4472135954999579*fl[38]-1.118033988749895*fl[37]-0.8660254037844386*fl[17]-0.5*fl[10])*dvxl; 
  incr[32] = (5.0*fl[32]+3.872983346207417*fl[15]+2.23606797749979*fl[7])*wxl+(fl[33]+0.5773502691896257*fl[22]+1.443375672974065*fl[21]+1.118033988749895*fl[6]+0.6454972243679029*fl[3])*dvxl; 
  incr[33] = ((-3.0*fl[33])-1.732050807568877*fl[22])*wxl+((-1.0*fl[32])-0.7745966692414833*fl[15]-0.4472135954999579*fl[7])*dvxl; 
  incr[34] = ((-3.0*fl[34])-1.732050807568877*fl[24])*wxl+((-0.8660254037844387*fl[23])-0.5*fl[13])*dvxl; 
  incr[35] = (5.0*fl[35]+3.872983346207417*fl[16]+2.23606797749979*fl[9])*wxl+(fl[36]+0.5773502691896257*fl[26]+1.443375672974065*fl[25]+1.118033988749895*fl[8]+0.6454972243679029*fl[4])*dvxl; 
  incr[36] = ((-3.0*fl[36])-1.732050807568877*fl[26])*wxl+((-1.0*fl[35])-0.7745966692414833*fl[16]-0.4472135954999579*fl[9])*dvxl; 
  incr[37] = (5.0*fl[37]+3.872983346207417*fl[17]+2.23606797749979*fl[10])*wxl+(1.443375672974065*fl[44]+1.118033988749895*fl[31]+0.6454972243679029*fl[18])*dvxl; 
  incr[38] = (1.732050807568877*fl[45]+fl[38])*wxl+(0.5773502691896257*fl[44]+0.4472135954999579*fl[31]+0.2581988897471612*fl[18])*dvxl; 
  incr[39] = ((-3.0*fl[39])-1.732050807568877*fl[27])*wxl+((-0.8660254037844387*fl[46])-0.5*fl[40])*dvxl; 
  incr[40] = (1.732050807568877*fl[46]+fl[40])*wxl+(0.5*fl[39]+0.2886751345948129*fl[27])*dvxl; 
  incr[41] = ((-3.0*fl[41])-1.732050807568877*fl[29])*wxl+((-0.8660254037844387*fl[28])-0.5*fl[14])*dvxl; 
  incr[42] = ((-3.0*fl[42])-1.732050807568877*fl[30])*wxl+((-0.8660254037844387*fl[47])-0.5*fl[43])*dvxl; 
  incr[43] = (1.732050807568877*fl[47]+fl[43])*wxl+(0.5*fl[42]+0.2886751345948129*fl[30])*dvxl; 
  incr[44] = (5.0*fl[44]+3.872983346207417*fl[31]+2.23606797749979*fl[18])*wxl+(fl[45]+0.5773502691896257*fl[38]+1.443375672974065*fl[37]+1.118033988749895*fl[17]+0.6454972243679028*fl[10])*dvxl; 
  incr[45] = ((-3.0*fl[45])-1.732050807568877*fl[38])*wxl+((-1.0*fl[44])-0.7745966692414834*fl[31]-0.447213595499958*fl[18])*dvxl; 
  incr[46] = ((-3.0*fl[46])-1.732050807568877*fl[40])*wxl+((-0.8660254037844387*fl[39])-0.5*fl[27])*dvxl; 
  incr[47] = ((-3.0*fl[47])-1.732050807568877*fl[43])*wxl+((-0.8660254037844387*fl[42])-0.5*fl[30])*dvxl; 

  outr[0] += incr[0]*dxr; 
  outr[1] += incr[1]*dxr; 
  outr[2] += incr[2]*dxr; 
  outr[3] += incr[3]*dxr; 
  outr[4] += incr[4]*dxr; 
  outr[5] += incr[5]*dxr; 
  outr[6] += incr[6]*dxr; 
  outr[7] += incr[7]*dxr; 
  outr[8] += incr[8]*dxr; 
  outr[9] += incr[9]*dxr; 
  outr[10] += incr[10]*dxr; 
  outr[11] += incr[11]*dxr; 
  outr[12] += incr[12]*dxr; 
  outr[13] += incr[13]*dxr; 
  outr[14] += incr[14]*dxr; 
  outr[15] += incr[15]*dxr; 
  outr[16] += incr[16]*dxr; 
  outr[17] += incr[17]*dxr; 
  outr[18] += incr[18]*dxr; 
  outr[19] += incr[19]*dxr; 
  outr[20] += incr[20]*dxr; 
  outr[21] += incr[21]*dxr; 
  outr[22] += incr[22]*dxr; 
  outr[23] += incr[23]*dxr; 
  outr[24] += incr[24]*dxr; 
  outr[25] += incr[25]*dxr; 
  outr[26] += incr[26]*dxr; 
  outr[27] += incr[27]*dxr; 
  outr[28] += incr[28]*dxr; 
  outr[29] += incr[29]*dxr; 
  outr[30] += incr[30]*dxr; 
  outr[31] += incr[31]*dxr; 
  outr[32] += incr[32]*dxr; 
  outr[33] += incr[33]*dxr; 
  outr[34] += incr[34]*dxr; 
  outr[35] += incr[35]*dxr; 
  outr[36] += incr[36]*dxr; 
  outr[37] += incr[37]*dxr; 
  outr[38] += incr[38]*dxr; 
  outr[39] += incr[39]*dxr; 
  outr[40] += incr[40]*dxr; 
  outr[41] += incr[41]*dxr; 
  outr[42] += incr[42]*dxr; 
  outr[43] += incr[43]*dxr; 
  outr[44] += incr[44]*dxr; 
  outr[45] += incr[45]*dxr; 
  outr[46] += incr[46]*dxr; 
  outr[47] += incr[47]*dxr; 

  outl[0] += -1.0*incr[0]*dxl; 
  outl[1] += incr[1]*dxl; 
  outl[2] += -1.0*incr[2]*dxl; 
  outl[3] += -1.0*incr[3]*dxl; 
  outl[4] += -1.0*incr[4]*dxl; 
  outl[5] += incr[5]*dxl; 
  outl[6] += incr[6]*dxl; 
  outl[7] += -1.0*incr[7]*dxl; 
  outl[8] += incr[8]*dxl; 
  outl[9] += -1.0*incr[9]*dxl; 
  outl[10] += -1.0*incr[10]*dxl; 
  outl[11] += -1.0*incr[11]*dxl; 
  outl[12] += -1.0*incr[12]*dxl; 
  outl[13] += -1.0*incr[13]*dxl; 
  outl[14] += -1.0*incr[14]*dxl; 
  outl[15] += incr[15]*dxl; 
  outl[16] += incr[16]*dxl; 
  outl[17] += incr[17]*dxl; 
  outl[18] += -1.0*incr[18]*dxl; 
  outl[19] += -1.0*incr[19]*dxl; 
  outl[20] += incr[20]*dxl; 
  outl[21] += -1.0*incr[21]*dxl; 
  outl[22] += -1.0*incr[22]*dxl; 
  outl[23] += incr[23]*dxl; 
  outl[24] += -1.0*incr[24]*dxl; 
  outl[25] += -1.0*incr[25]*dxl; 
  outl[26] += -1.0*incr[26]*dxl; 
  outl[27] += -1.0*incr[27]*dxl; 
  outl[28] += incr[28]*dxl; 
  outl[29] += -1.0*incr[29]*dxl; 
  outl[30] += -1.0*incr[30]*dxl; 
  outl[31] += incr[31]*dxl; 
  outl[32] += -1.0*incr[32]*dxl; 
  outl[33] += incr[33]*dxl; 
  outl[34] += incr[34]*dxl; 
  outl[35] += -1.0*incr[35]*dxl; 
  outl[36] += incr[36]*dxl; 
  outl[37] += -1.0*incr[37]*dxl; 
  outl[38] += -1.0*incr[38]*dxl; 
  outl[39] += incr[39]*dxl; 
  outl[40] += -1.0*incr[40]*dxl; 
  outl[41] += incr[41]*dxl; 
  outl[42] += incr[42]*dxl; 
  outl[43] += -1.0*incr[43]*dxl; 
  outl[44] += -1.0*incr[44]*dxl; 
  outl[45] += incr[45]*dxl; 
  outl[46] += incr[46]*dxl; 
  outl[47] += incr[47]*dxl; 
  } else { 
  incr[0] = (2.23606797749979*fr[11]-1.732050807568877*fr[1]+fr[0])*wxr+(0.6454972243679028*fr[19]-0.5*fr[5]+0.2886751345948129*fr[2])*dvxr; 
  incr[1] = ((-3.872983346207417*fr[11])+3.0*fr[1]-1.732050807568877*fr[0])*wxr+((-1.118033988749895*fr[19])+0.8660254037844386*fr[5]-0.5*fr[2])*dvxr; 
  incr[2] = (2.23606797749979*fr[19]-1.732050807568877*fr[5]+fr[2])*wxr+((-0.447213595499958*fr[20])+0.2581988897471612*fr[12]+0.6454972243679029*fr[11]-0.5*fr[1]+0.2886751345948129*fr[0])*dvxr; 
  incr[3] = (2.23606797749979*fr[21]-1.732050807568877*fr[6]+fr[3])*wxr+(0.6454972243679029*fr[32]-0.5*fr[15]+0.2886751345948129*fr[7])*dvxr; 
  incr[4] = (2.23606797749979*fr[25]-1.732050807568877*fr[8]+fr[4])*wxr+(0.6454972243679029*fr[35]-0.5*fr[16]+0.2886751345948129*fr[9])*dvxr; 
  incr[5] = ((-3.872983346207417*fr[19])+3.0*fr[5]-1.732050807568877*fr[2])*wxr+(0.7745966692414834*fr[20]-0.4472135954999579*fr[12]-1.118033988749895*fr[11]+0.8660254037844386*fr[1]-0.5*fr[0])*dvxr; 
  incr[6] = ((-3.872983346207417*fr[21])+3.0*fr[6]-1.732050807568877*fr[3])*wxr+((-1.118033988749895*fr[32])+0.8660254037844386*fr[15]-0.5*fr[7])*dvxr; 
  incr[7] = (2.23606797749979*fr[32]-1.732050807568877*fr[15]+fr[7])*wxr+((-0.4472135954999579*fr[33])+0.2581988897471611*fr[22]+0.6454972243679028*fr[21]-0.5*fr[6]+0.2886751345948129*fr[3])*dvxr; 
  incr[8] = ((-3.872983346207417*fr[25])+3.0*fr[8]-1.732050807568877*fr[4])*wxr+((-1.118033988749895*fr[35])+0.8660254037844386*fr[16]-0.5*fr[9])*dvxr; 
  incr[9] = (2.23606797749979*fr[35]-1.732050807568877*fr[16]+fr[9])*wxr+((-0.4472135954999579*fr[36])+0.2581988897471611*fr[26]+0.6454972243679028*fr[25]-0.5*fr[8]+0.2886751345948129*fr[4])*dvxr; 
  incr[10] = (2.23606797749979*fr[37]-1.732050807568877*fr[17]+fr[10])*wxr+(0.6454972243679028*fr[44]-0.5*fr[31]+0.2886751345948129*fr[18])*dvxr; 
  incr[11] = (5.0*fr[11]-3.872983346207417*fr[1]+2.23606797749979*fr[0])*wxr+(1.443375672974065*fr[19]-1.118033988749895*fr[5]+0.6454972243679029*fr[2])*dvxr; 
  incr[12] = (fr[12]-1.732050807568877*fr[20])*wxr+(0.5773502691896257*fr[19]-0.4472135954999579*fr[5]+0.2581988897471612*fr[2])*dvxr; 
  incr[13] = (fr[13]-1.732050807568877*fr[23])*wxr+(0.2886751345948129*fr[24]-0.5*fr[34])*dvxr; 
  incr[14] = (fr[14]-1.732050807568877*fr[28])*wxr+(0.2886751345948129*fr[29]-0.5*fr[41])*dvxr; 
  incr[15] = ((-3.872983346207417*fr[32])+3.0*fr[15]-1.732050807568877*fr[7])*wxr+(0.7745966692414833*fr[33]-0.447213595499958*fr[22]-1.118033988749895*fr[21]+0.8660254037844386*fr[6]-0.5*fr[3])*dvxr; 
  incr[16] = ((-3.872983346207417*fr[35])+3.0*fr[16]-1.732050807568877*fr[9])*wxr+(0.7745966692414833*fr[36]-0.447213595499958*fr[26]-1.118033988749895*fr[25]+0.8660254037844386*fr[8]-0.5*fr[4])*dvxr; 
  incr[17] = ((-3.872983346207417*fr[37])+3.0*fr[17]-1.732050807568877*fr[10])*wxr+((-1.118033988749895*fr[44])+0.8660254037844386*fr[31]-0.5*fr[18])*dvxr; 
  incr[18] = (2.23606797749979*fr[44]-1.732050807568877*fr[31]+fr[18])*wxr+((-0.447213595499958*fr[45])+0.2581988897471612*fr[38]+0.6454972243679029*fr[37]-0.5*fr[17]+0.2886751345948129*fr[10])*dvxr; 
  incr[19] = (5.0*fr[19]-3.872983346207417*fr[5]+2.23606797749979*fr[2])*wxr+((-1.0*fr[20])+0.5773502691896257*fr[12]+1.443375672974065*fr[11]-1.118033988749895*fr[1]+0.6454972243679028*fr[0])*dvxr; 
  incr[20] = (3.0*fr[20]-1.732050807568877*fr[12])*wxr+((-1.0*fr[19])+0.7745966692414834*fr[5]-0.447213595499958*fr[2])*dvxr; 
  incr[21] = (5.0*fr[21]-3.872983346207417*fr[6]+2.23606797749979*fr[3])*wxr+(1.443375672974065*fr[32]-1.118033988749895*fr[15]+0.6454972243679028*fr[7])*dvxr; 
  incr[22] = (fr[22]-1.732050807568877*fr[33])*wxr+(0.5773502691896257*fr[32]-0.447213595499958*fr[15]+0.2581988897471611*fr[7])*dvxr; 
  incr[23] = (3.0*fr[23]-1.732050807568877*fr[13])*wxr+(0.8660254037844387*fr[34]-0.5*fr[24])*dvxr; 
  incr[24] = (fr[24]-1.732050807568877*fr[34])*wxr+(0.2886751345948129*fr[13]-0.5*fr[23])*dvxr; 
  incr[25] = (5.0*fr[25]-3.872983346207417*fr[8]+2.23606797749979*fr[4])*wxr+(1.443375672974065*fr[35]-1.118033988749895*fr[16]+0.6454972243679028*fr[9])*dvxr; 
  incr[26] = (fr[26]-1.732050807568877*fr[36])*wxr+(0.5773502691896257*fr[35]-0.447213595499958*fr[16]+0.2581988897471611*fr[9])*dvxr; 
  incr[27] = (fr[27]-1.732050807568877*fr[39])*wxr+(0.2886751345948129*fr[40]-0.5*fr[46])*dvxr; 
  incr[28] = (3.0*fr[28]-1.732050807568877*fr[14])*wxr+(0.8660254037844387*fr[41]-0.5*fr[29])*dvxr; 
  incr[29] = (fr[29]-1.732050807568877*fr[41])*wxr+(0.2886751345948129*fr[14]-0.5*fr[28])*dvxr; 
  incr[30] = (fr[30]-1.732050807568877*fr[42])*wxr+(0.2886751345948129*fr[43]-0.5*fr[47])*dvxr; 
  incr[31] = ((-3.872983346207417*fr[44])+3.0*fr[31]-1.732050807568877*fr[18])*wxr+(0.7745966692414834*fr[45]-0.4472135954999579*fr[38]-1.118033988749895*fr[37]+0.8660254037844386*fr[17]-0.5*fr[10])*dvxr; 
  incr[32] = (5.0*fr[32]-3.872983346207417*fr[15]+2.23606797749979*fr[7])*wxr+((-1.0*fr[33])+0.5773502691896257*fr[22]+1.443375672974065*fr[21]-1.118033988749895*fr[6]+0.6454972243679029*fr[3])*dvxr; 
  incr[33] = (3.0*fr[33]-1.732050807568877*fr[22])*wxr+((-1.0*fr[32])+0.7745966692414833*fr[15]-0.4472135954999579*fr[7])*dvxr; 
  incr[34] = (3.0*fr[34]-1.732050807568877*fr[24])*wxr+(0.8660254037844387*fr[23]-0.5*fr[13])*dvxr; 
  incr[35] = (5.0*fr[35]-3.872983346207417*fr[16]+2.23606797749979*fr[9])*wxr+((-1.0*fr[36])+0.5773502691896257*fr[26]+1.443375672974065*fr[25]-1.118033988749895*fr[8]+0.6454972243679029*fr[4])*dvxr; 
  incr[36] = (3.0*fr[36]-1.732050807568877*fr[26])*wxr+((-1.0*fr[35])+0.7745966692414833*fr[16]-0.4472135954999579*fr[9])*dvxr; 
  incr[37] = (5.0*fr[37]-3.872983346207417*fr[17]+2.23606797749979*fr[10])*wxr+(1.443375672974065*fr[44]-1.118033988749895*fr[31]+0.6454972243679029*fr[18])*dvxr; 
  incr[38] = (fr[38]-1.732050807568877*fr[45])*wxr+(0.5773502691896257*fr[44]-0.4472135954999579*fr[31]+0.2581988897471612*fr[18])*dvxr; 
  incr[39] = (3.0*fr[39]-1.732050807568877*fr[27])*wxr+(0.8660254037844387*fr[46]-0.5*fr[40])*dvxr; 
  incr[40] = (fr[40]-1.732050807568877*fr[46])*wxr+(0.2886751345948129*fr[27]-0.5*fr[39])*dvxr; 
  incr[41] = (3.0*fr[41]-1.732050807568877*fr[29])*wxr+(0.8660254037844387*fr[28]-0.5*fr[14])*dvxr; 
  incr[42] = (3.0*fr[42]-1.732050807568877*fr[30])*wxr+(0.8660254037844387*fr[47]-0.5*fr[43])*dvxr; 
  incr[43] = (fr[43]-1.732050807568877*fr[47])*wxr+(0.2886751345948129*fr[30]-0.5*fr[42])*dvxr; 
  incr[44] = (5.0*fr[44]-3.872983346207417*fr[31]+2.23606797749979*fr[18])*wxr+((-1.0*fr[45])+0.5773502691896257*fr[38]+1.443375672974065*fr[37]-1.118033988749895*fr[17]+0.6454972243679028*fr[10])*dvxr; 
  incr[45] = (3.0*fr[45]-1.732050807568877*fr[38])*wxr+((-1.0*fr[44])+0.7745966692414834*fr[31]-0.447213595499958*fr[18])*dvxr; 
  incr[46] = (3.0*fr[46]-1.732050807568877*fr[40])*wxr+(0.8660254037844387*fr[39]-0.5*fr[27])*dvxr; 
  incr[47] = (3.0*fr[47]-1.732050807568877*fr[43])*wxr+(0.8660254037844387*fr[42]-0.5*fr[30])*dvxr; 

  outr[0] += incr[0]*dxr; 
  outr[1] += incr[1]*dxr; 
  outr[2] += incr[2]*dxr; 
  outr[3] += incr[3]*dxr; 
  outr[4] += incr[4]*dxr; 
  outr[5] += incr[5]*dxr; 
  outr[6] += incr[6]*dxr; 
  outr[7] += incr[7]*dxr; 
  outr[8] += incr[8]*dxr; 
  outr[9] += incr[9]*dxr; 
  outr[10] += incr[10]*dxr; 
  outr[11] += incr[11]*dxr; 
  outr[12] += incr[12]*dxr; 
  outr[13] += incr[13]*dxr; 
  outr[14] += incr[14]*dxr; 
  outr[15] += incr[15]*dxr; 
  outr[16] += incr[16]*dxr; 
  outr[17] += incr[17]*dxr; 
  outr[18] += incr[18]*dxr; 
  outr[19] += incr[19]*dxr; 
  outr[20] += incr[20]*dxr; 
  outr[21] += incr[21]*dxr; 
  outr[22] += incr[22]*dxr; 
  outr[23] += incr[23]*dxr; 
  outr[24] += incr[24]*dxr; 
  outr[25] += incr[25]*dxr; 
  outr[26] += incr[26]*dxr; 
  outr[27] += incr[27]*dxr; 
  outr[28] += incr[28]*dxr; 
  outr[29] += incr[29]*dxr; 
  outr[30] += incr[30]*dxr; 
  outr[31] += incr[31]*dxr; 
  outr[32] += incr[32]*dxr; 
  outr[33] += incr[33]*dxr; 
  outr[34] += incr[34]*dxr; 
  outr[35] += incr[35]*dxr; 
  outr[36] += incr[36]*dxr; 
  outr[37] += incr[37]*dxr; 
  outr[38] += incr[38]*dxr; 
  outr[39] += incr[39]*dxr; 
  outr[40] += incr[40]*dxr; 
  outr[41] += incr[41]*dxr; 
  outr[42] += incr[42]*dxr; 
  outr[43] += incr[43]*dxr; 
  outr[44] += incr[44]*dxr; 
  outr[45] += incr[45]*dxr; 
  outr[46] += incr[46]*dxr; 
  outr[47] += incr[47]*dxr; 

  outl[0] += -1.0*incr[0]*dxl; 
  outl[1] += incr[1]*dxl; 
  outl[2] += -1.0*incr[2]*dxl; 
  outl[3] += -1.0*incr[3]*dxl; 
  outl[4] += -1.0*incr[4]*dxl; 
  outl[5] += incr[5]*dxl; 
  outl[6] += incr[6]*dxl; 
  outl[7] += -1.0*incr[7]*dxl; 
  outl[8] += incr[8]*dxl; 
  outl[9] += -1.0*incr[9]*dxl; 
  outl[10] += -1.0*incr[10]*dxl; 
  outl[11] += -1.0*incr[11]*dxl; 
  outl[12] += -1.0*incr[12]*dxl; 
  outl[13] += -1.0*incr[13]*dxl; 
  outl[14] += -1.0*incr[14]*dxl; 
  outl[15] += incr[15]*dxl; 
  outl[16] += incr[16]*dxl; 
  outl[17] += incr[17]*dxl; 
  outl[18] += -1.0*incr[18]*dxl; 
  outl[19] += -1.0*incr[19]*dxl; 
  outl[20] += incr[20]*dxl; 
  outl[21] += -1.0*incr[21]*dxl; 
  outl[22] += -1.0*incr[22]*dxl; 
  outl[23] += incr[23]*dxl; 
  outl[24] += -1.0*incr[24]*dxl; 
  outl[25] += -1.0*incr[25]*dxl; 
  outl[26] += -1.0*incr[26]*dxl; 
  outl[27] += -1.0*incr[27]*dxl; 
  outl[28] += incr[28]*dxl; 
  outl[29] += -1.0*incr[29]*dxl; 
  outl[30] += -1.0*incr[30]*dxl; 
  outl[31] += incr[31]*dxl; 
  outl[32] += -1.0*incr[32]*dxl; 
  outl[33] += incr[33]*dxl; 
  outl[34] += incr[34]*dxl; 
  outl[35] += -1.0*incr[35]*dxl; 
  outl[36] += incr[36]*dxl; 
  outl[37] += -1.0*incr[37]*dxl; 
  outl[38] += -1.0*incr[38]*dxl; 
  outl[39] += incr[39]*dxl; 
  outl[40] += -1.0*incr[40]*dxl; 
  outl[41] += incr[41]*dxl; 
  outl[42] += incr[42]*dxl; 
  outl[43] += -1.0*incr[43]*dxl; 
  outl[44] += -1.0*incr[44]*dxl; 
  outl[45] += incr[45]*dxl; 
  outl[46] += incr[46]*dxl; 
  outl[47] += incr[47]*dxl; 
  } 
} 
