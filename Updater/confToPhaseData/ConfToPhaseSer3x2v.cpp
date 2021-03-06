#include <ConfToPhaseModDecl.h>
void accumulateConfToPhase3x2vSer_P1(const double fact, const double *fconf, double *fphase) 
{ 
  fphase[0] += 2.0*fconf[0]*fact; 
  fphase[1] += 2.0*fconf[1]*fact; 
  fphase[2] += 2.0*fconf[2]*fact; 
  fphase[3] += 2.0*fconf[3]*fact; 
  fphase[6] += 2.0*fconf[4]*fact; 
  fphase[7] += 2.0*fconf[5]*fact; 
  fphase[8] += 2.0*fconf[6]*fact; 
  fphase[16] += 2.0*fconf[7]*fact; 
} 
void accumulateConfToPhase3x2vSer_P2(const double fact, const double *fconf, double *fphase) 
{ 
  fphase[0] += 2.0*fconf[0]*fact; 
  fphase[1] += 2.0*fconf[1]*fact; 
  fphase[2] += 2.0*fconf[2]*fact; 
  fphase[3] += 2.0*fconf[3]*fact; 
  fphase[6] += 2.0*fconf[4]*fact; 
  fphase[7] += 2.0*fconf[5]*fact; 
  fphase[8] += 2.0*fconf[6]*fact; 
  fphase[16] += 2.0*fconf[7]*fact; 
  fphase[17] += 2.0*fconf[8]*fact; 
  fphase[18] += 2.0*fconf[9]*fact; 
  fphase[21] += 2.0*fconf[10]*fact; 
  fphase[31] += 2.0*fconf[11]*fact; 
  fphase[32] += 2.0*fconf[12]*fact; 
  fphase[33] += 2.0*fconf[13]*fact; 
  fphase[34] += 2.0*fconf[14]*fact; 
  fphase[35] += 2.0*fconf[15]*fact; 
  fphase[36] += 2.0*fconf[16]*fact; 
  fphase[56] += 2.0*fconf[17]*fact; 
  fphase[57] += 2.0*fconf[18]*fact; 
  fphase[58] += 2.0*fconf[19]*fact; 
} 
void assignConfToPhase3x2vSer_P1(const double fact, const double *fconf, double *fphase) 
{ 
  fphase[0] = 2.0*fconf[0]*fact; 
  fphase[1] = 2.0*fconf[1]*fact; 
  fphase[2] = 2.0*fconf[2]*fact; 
  fphase[3] = 2.0*fconf[3]*fact; 
  fphase[4] = 0.0; 
  fphase[5] = 0.0; 
  fphase[6] = 2.0*fconf[4]*fact; 
  fphase[7] = 2.0*fconf[5]*fact; 
  fphase[8] = 2.0*fconf[6]*fact; 
  fphase[9] = 0.0; 
  fphase[10] = 0.0; 
  fphase[11] = 0.0; 
  fphase[12] = 0.0; 
  fphase[13] = 0.0; 
  fphase[14] = 0.0; 
  fphase[15] = 0.0; 
  fphase[16] = 2.0*fconf[7]*fact; 
  fphase[17] = 0.0; 
  fphase[18] = 0.0; 
  fphase[19] = 0.0; 
  fphase[20] = 0.0; 
  fphase[21] = 0.0; 
  fphase[22] = 0.0; 
  fphase[23] = 0.0; 
  fphase[24] = 0.0; 
  fphase[25] = 0.0; 
  fphase[26] = 0.0; 
  fphase[27] = 0.0; 
  fphase[28] = 0.0; 
  fphase[29] = 0.0; 
  fphase[30] = 0.0; 
  fphase[31] = 0.0; 
} 
void assignConfToPhase3x2vSer_P2(const double fact, const double *fconf, double *fphase) 
{ 
  fphase[0] = 2.0*fconf[0]*fact; 
  fphase[1] = 2.0*fconf[1]*fact; 
  fphase[2] = 2.0*fconf[2]*fact; 
  fphase[3] = 2.0*fconf[3]*fact; 
  fphase[4] = 0.0; 
  fphase[5] = 0.0; 
  fphase[6] = 2.0*fconf[4]*fact; 
  fphase[7] = 2.0*fconf[5]*fact; 
  fphase[8] = 2.0*fconf[6]*fact; 
  fphase[9] = 0.0; 
  fphase[10] = 0.0; 
  fphase[11] = 0.0; 
  fphase[12] = 0.0; 
  fphase[13] = 0.0; 
  fphase[14] = 0.0; 
  fphase[15] = 0.0; 
  fphase[16] = 2.0*fconf[7]*fact; 
  fphase[17] = 2.0*fconf[8]*fact; 
  fphase[18] = 2.0*fconf[9]*fact; 
  fphase[19] = 0.0; 
  fphase[20] = 0.0; 
  fphase[21] = 2.0*fconf[10]*fact; 
  fphase[22] = 0.0; 
  fphase[23] = 0.0; 
  fphase[24] = 0.0; 
  fphase[25] = 0.0; 
  fphase[26] = 0.0; 
  fphase[27] = 0.0; 
  fphase[28] = 0.0; 
  fphase[29] = 0.0; 
  fphase[30] = 0.0; 
  fphase[31] = 2.0*fconf[11]*fact; 
  fphase[32] = 2.0*fconf[12]*fact; 
  fphase[33] = 2.0*fconf[13]*fact; 
  fphase[34] = 2.0*fconf[14]*fact; 
  fphase[35] = 2.0*fconf[15]*fact; 
  fphase[36] = 2.0*fconf[16]*fact; 
  fphase[37] = 0.0; 
  fphase[38] = 0.0; 
  fphase[39] = 0.0; 
  fphase[40] = 0.0; 
  fphase[41] = 0.0; 
  fphase[42] = 0.0; 
  fphase[43] = 0.0; 
  fphase[44] = 0.0; 
  fphase[45] = 0.0; 
  fphase[46] = 0.0; 
  fphase[47] = 0.0; 
  fphase[48] = 0.0; 
  fphase[49] = 0.0; 
  fphase[50] = 0.0; 
  fphase[51] = 0.0; 
  fphase[52] = 0.0; 
  fphase[53] = 0.0; 
  fphase[54] = 0.0; 
  fphase[55] = 0.0; 
  fphase[56] = 2.0*fconf[17]*fact; 
  fphase[57] = 2.0*fconf[18]*fact; 
  fphase[58] = 2.0*fconf[19]*fact; 
  fphase[59] = 0.0; 
  fphase[60] = 0.0; 
  fphase[61] = 0.0; 
  fphase[62] = 0.0; 
  fphase[63] = 0.0; 
  fphase[64] = 0.0; 
  fphase[65] = 0.0; 
  fphase[66] = 0.0; 
  fphase[67] = 0.0; 
  fphase[68] = 0.0; 
  fphase[69] = 0.0; 
  fphase[70] = 0.0; 
  fphase[71] = 0.0; 
  fphase[72] = 0.0; 
  fphase[73] = 0.0; 
  fphase[74] = 0.0; 
  fphase[75] = 0.0; 
  fphase[76] = 0.0; 
  fphase[77] = 0.0; 
  fphase[78] = 0.0; 
  fphase[79] = 0.0; 
  fphase[80] = 0.0; 
  fphase[81] = 0.0; 
  fphase[82] = 0.0; 
  fphase[83] = 0.0; 
  fphase[84] = 0.0; 
  fphase[85] = 0.0; 
  fphase[86] = 0.0; 
  fphase[87] = 0.0; 
  fphase[88] = 0.0; 
  fphase[89] = 0.0; 
  fphase[90] = 0.0; 
  fphase[91] = 0.0; 
  fphase[92] = 0.0; 
  fphase[93] = 0.0; 
  fphase[94] = 0.0; 
  fphase[95] = 0.0; 
  fphase[96] = 0.0; 
  fphase[97] = 0.0; 
  fphase[98] = 0.0; 
  fphase[99] = 0.0; 
  fphase[100] = 0.0; 
  fphase[101] = 0.0; 
  fphase[102] = 0.0; 
  fphase[103] = 0.0; 
  fphase[104] = 0.0; 
  fphase[105] = 0.0; 
  fphase[106] = 0.0; 
  fphase[107] = 0.0; 
  fphase[108] = 0.0; 
  fphase[109] = 0.0; 
  fphase[110] = 0.0; 
  fphase[111] = 0.0; 
}
