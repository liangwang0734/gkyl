#include <MaxwellianOnBasisModDecl.h>

void MaxwellianOnBasisGauss2x2vMax_P2_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, double *flowUOrd, double *vtSqOrd, double *fMFacOrd) {

  flowUOrd[0] = 0.5*flowU[0]-0.5590169943749475*(flowU[5]+flowU[4]); 
  flowUOrd[1] = 0.4472135954999581*flowU[5]-0.5590169943749475*flowU[4]-0.6708203932499369*flowU[2]+0.5*flowU[0]; 
  flowUOrd[2] = 0.4472135954999581*flowU[5]-0.5590169943749475*flowU[4]+0.6708203932499369*flowU[2]+0.5*flowU[0]; 
  flowUOrd[3] = (-0.5590169943749475*flowU[5])+0.4472135954999581*flowU[4]-0.6708203932499369*flowU[1]+0.5*flowU[0]; 
  flowUOrd[4] = 0.4472135954999581*(flowU[5]+flowU[4])+0.9*flowU[3]-0.6708203932499369*(flowU[2]+flowU[1])+0.5*flowU[0]; 
  flowUOrd[5] = 0.4472135954999581*(flowU[5]+flowU[4])-0.9*flowU[3]+0.6708203932499369*flowU[2]-0.6708203932499369*flowU[1]+0.5*flowU[0]; 
  flowUOrd[6] = (-0.5590169943749475*flowU[5])+0.4472135954999581*flowU[4]+0.6708203932499369*flowU[1]+0.5*flowU[0]; 
  flowUOrd[7] = 0.4472135954999581*(flowU[5]+flowU[4])-0.9*flowU[3]-0.6708203932499369*flowU[2]+0.6708203932499369*flowU[1]+0.5*flowU[0]; 
  flowUOrd[8] = 0.4472135954999581*(flowU[5]+flowU[4])+0.9*flowU[3]+0.6708203932499369*(flowU[2]+flowU[1])+0.5*flowU[0]; 
  flowUOrd[9] = 0.5*flowU[9]-0.5590169943749475*(flowU[14]+flowU[13]); 
  flowUOrd[10] = 0.4472135954999581*flowU[14]-0.5590169943749475*flowU[13]-0.6708203932499369*flowU[11]+0.5*flowU[9]; 
  flowUOrd[11] = 0.4472135954999581*flowU[14]-0.5590169943749475*flowU[13]+0.6708203932499369*flowU[11]+0.5*flowU[9]; 
  flowUOrd[12] = (-0.5590169943749475*flowU[14])+0.4472135954999581*flowU[13]-0.6708203932499369*flowU[10]+0.5*flowU[9]; 
  flowUOrd[13] = 0.4472135954999581*(flowU[14]+flowU[13])+0.9*flowU[12]-0.6708203932499369*(flowU[11]+flowU[10])+0.5*flowU[9]; 
  flowUOrd[14] = 0.4472135954999581*(flowU[14]+flowU[13])-0.9*flowU[12]+0.6708203932499369*flowU[11]-0.6708203932499369*flowU[10]+0.5*flowU[9]; 
  flowUOrd[15] = (-0.5590169943749475*flowU[14])+0.4472135954999581*flowU[13]+0.6708203932499369*flowU[10]+0.5*flowU[9]; 
  flowUOrd[16] = 0.4472135954999581*(flowU[14]+flowU[13])-0.9*flowU[12]-0.6708203932499369*flowU[11]+0.6708203932499369*flowU[10]+0.5*flowU[9]; 
  flowUOrd[17] = 0.4472135954999581*(flowU[14]+flowU[13])+0.9*flowU[12]+0.6708203932499369*(flowU[11]+flowU[10])+0.5*flowU[9]; 

  vtSqOrd[0] = 0.5*vtSq[0]-0.5590169943749475*(vtSq[5]+vtSq[4]); 
  vtSqOrd[1] = 0.4472135954999581*vtSq[5]-0.5590169943749475*vtSq[4]-0.6708203932499369*vtSq[2]+0.5*vtSq[0]; 
  vtSqOrd[2] = 0.4472135954999581*vtSq[5]-0.5590169943749475*vtSq[4]+0.6708203932499369*vtSq[2]+0.5*vtSq[0]; 
  vtSqOrd[3] = (-0.5590169943749475*vtSq[5])+0.4472135954999581*vtSq[4]-0.6708203932499369*vtSq[1]+0.5*vtSq[0]; 
  vtSqOrd[4] = 0.4472135954999581*(vtSq[5]+vtSq[4])+0.9*vtSq[3]-0.6708203932499369*(vtSq[2]+vtSq[1])+0.5*vtSq[0]; 
  vtSqOrd[5] = 0.4472135954999581*(vtSq[5]+vtSq[4])-0.9*vtSq[3]+0.6708203932499369*vtSq[2]-0.6708203932499369*vtSq[1]+0.5*vtSq[0]; 
  vtSqOrd[6] = (-0.5590169943749475*vtSq[5])+0.4472135954999581*vtSq[4]+0.6708203932499369*vtSq[1]+0.5*vtSq[0]; 
  vtSqOrd[7] = 0.4472135954999581*(vtSq[5]+vtSq[4])-0.9*vtSq[3]-0.6708203932499369*vtSq[2]+0.6708203932499369*vtSq[1]+0.5*vtSq[0]; 
  vtSqOrd[8] = 0.4472135954999581*(vtSq[5]+vtSq[4])+0.9*vtSq[3]+0.6708203932499369*(vtSq[2]+vtSq[1])+0.5*vtSq[0]; 

  if (vtSqOrd[0] <= 0.0)
    fMFacOrd[0] = 0.;
  else
    fMFacOrd[0] = (0.1591549430918953*(0.5*den[0]-0.5590169943749475*(den[5]+den[4])))/vtSqOrd[0]; 
  if (vtSqOrd[1] <= 0.0)
    fMFacOrd[1] = 0.;
  else
    fMFacOrd[1] = (0.1591549430918953*(0.4472135954999581*den[5]-0.5590169943749475*den[4]-0.6708203932499369*den[2]+0.5*den[0]))/vtSqOrd[1]; 
  if (vtSqOrd[2] <= 0.0)
    fMFacOrd[2] = 0.;
  else
    fMFacOrd[2] = (0.1591549430918953*(0.4472135954999581*den[5]-0.5590169943749475*den[4]+0.6708203932499369*den[2]+0.5*den[0]))/vtSqOrd[2]; 
  if (vtSqOrd[3] <= 0.0)
    fMFacOrd[3] = 0.;
  else
    fMFacOrd[3] = (0.1591549430918953*((-0.5590169943749475*den[5])+0.4472135954999581*den[4]-0.6708203932499369*den[1]+0.5*den[0]))/vtSqOrd[3]; 
  if (vtSqOrd[4] <= 0.0)
    fMFacOrd[4] = 0.;
  else
    fMFacOrd[4] = (0.1591549430918953*(0.4472135954999581*(den[5]+den[4])+0.9*den[3]-0.6708203932499369*(den[2]+den[1])+0.5*den[0]))/vtSqOrd[4]; 
  if (vtSqOrd[5] <= 0.0)
    fMFacOrd[5] = 0.;
  else
    fMFacOrd[5] = (0.1591549430918953*(0.4472135954999581*(den[5]+den[4])-0.9*den[3]+0.6708203932499369*den[2]-0.6708203932499369*den[1]+0.5*den[0]))/vtSqOrd[5]; 
  if (vtSqOrd[6] <= 0.0)
    fMFacOrd[6] = 0.;
  else
    fMFacOrd[6] = (0.1591549430918953*((-0.5590169943749475*den[5])+0.4472135954999581*den[4]+0.6708203932499369*den[1]+0.5*den[0]))/vtSqOrd[6]; 
  if (vtSqOrd[7] <= 0.0)
    fMFacOrd[7] = 0.;
  else
    fMFacOrd[7] = (0.1591549430918953*(0.4472135954999581*(den[5]+den[4])-0.9*den[3]-0.6708203932499369*den[2]+0.6708203932499369*den[1]+0.5*den[0]))/vtSqOrd[7]; 
  if (vtSqOrd[8] <= 0.0)
    fMFacOrd[8] = 0.;
  else
    fMFacOrd[8] = (0.1591549430918953*(0.4472135954999581*(den[5]+den[4])+0.9*den[3]+0.6708203932499369*(den[2]+den[1])+0.5*den[0]))/vtSqOrd[8]; 

}
void MaxwellianOnBasisGauss2x2vMax_P2_phaseQuad(const double *flowUOrd, const double *vtSqOrd, const double *fMFacOrd, const double *wc, const double *dxv, double *fMOut) {

  double fMquad[81];
  fMquad[0] = fMFacOrd[0]*exp(-(0.5*(std::pow(wc[3]-1.0*flowUOrd[9],2.0)+std::pow(wc[2]-1.0*flowUOrd[0],2.0)))/vtSqOrd[0]); 
  fMquad[1] = fMFacOrd[0]*exp(-(0.5*(std::pow((-1.0*flowUOrd[9])+wc[3]-0.3872983346207417*dxv[3],2.0)+std::pow(wc[2]-1.0*flowUOrd[0],2.0)))/vtSqOrd[0]); 
  fMquad[2] = fMFacOrd[0]*exp(-(0.5*(std::pow((-1.0*flowUOrd[9])+wc[3]+0.3872983346207417*dxv[3],2.0)+std::pow(wc[2]-1.0*flowUOrd[0],2.0)))/vtSqOrd[0]); 
  fMquad[3] = fMFacOrd[0]*exp(-(0.5*(std::pow(wc[3]-1.0*flowUOrd[9],2.0)+std::pow(wc[2]-0.3872983346207417*dxv[2]-1.0*flowUOrd[0],2.0)))/vtSqOrd[0]); 
  fMquad[4] = fMFacOrd[0]*exp(-(0.5*(std::pow((-1.0*flowUOrd[9])+wc[3]-0.3872983346207417*dxv[3],2.0)+std::pow(wc[2]-0.3872983346207417*dxv[2]-1.0*flowUOrd[0],2.0)))/vtSqOrd[0]); 
  fMquad[5] = fMFacOrd[0]*exp(-(0.5*(std::pow((-1.0*flowUOrd[9])+wc[3]+0.3872983346207417*dxv[3],2.0)+std::pow(wc[2]-0.3872983346207417*dxv[2]-1.0*flowUOrd[0],2.0)))/vtSqOrd[0]); 
  fMquad[6] = fMFacOrd[0]*exp(-(0.5*(std::pow(wc[3]-1.0*flowUOrd[9],2.0)+std::pow(wc[2]+0.3872983346207417*dxv[2]-1.0*flowUOrd[0],2.0)))/vtSqOrd[0]); 
  fMquad[7] = fMFacOrd[0]*exp(-(0.5*(std::pow((-1.0*flowUOrd[9])+wc[3]-0.3872983346207417*dxv[3],2.0)+std::pow(wc[2]+0.3872983346207417*dxv[2]-1.0*flowUOrd[0],2.0)))/vtSqOrd[0]); 
  fMquad[8] = fMFacOrd[0]*exp(-(0.5*(std::pow((-1.0*flowUOrd[9])+wc[3]+0.3872983346207417*dxv[3],2.0)+std::pow(wc[2]+0.3872983346207417*dxv[2]-1.0*flowUOrd[0],2.0)))/vtSqOrd[0]); 
  fMquad[9] = fMFacOrd[1]*exp(-(0.5*(std::pow(wc[3]-1.0*flowUOrd[10],2.0)+std::pow(wc[2]-1.0*flowUOrd[1],2.0)))/vtSqOrd[1]); 
  fMquad[10] = fMFacOrd[1]*exp(-(0.5*(std::pow((-1.0*flowUOrd[10])+wc[3]-0.3872983346207417*dxv[3],2.0)+std::pow(wc[2]-1.0*flowUOrd[1],2.0)))/vtSqOrd[1]); 
  fMquad[11] = fMFacOrd[1]*exp(-(0.5*(std::pow((-1.0*flowUOrd[10])+wc[3]+0.3872983346207417*dxv[3],2.0)+std::pow(wc[2]-1.0*flowUOrd[1],2.0)))/vtSqOrd[1]); 
  fMquad[12] = fMFacOrd[1]*exp(-(0.5*(std::pow(wc[3]-1.0*flowUOrd[10],2.0)+std::pow(wc[2]-0.3872983346207417*dxv[2]-1.0*flowUOrd[1],2.0)))/vtSqOrd[1]); 
  fMquad[13] = fMFacOrd[1]*exp(-(0.5*(std::pow((-1.0*flowUOrd[10])+wc[3]-0.3872983346207417*dxv[3],2.0)+std::pow(wc[2]-0.3872983346207417*dxv[2]-1.0*flowUOrd[1],2.0)))/vtSqOrd[1]); 
  fMquad[14] = fMFacOrd[1]*exp(-(0.5*(std::pow((-1.0*flowUOrd[10])+wc[3]+0.3872983346207417*dxv[3],2.0)+std::pow(wc[2]-0.3872983346207417*dxv[2]-1.0*flowUOrd[1],2.0)))/vtSqOrd[1]); 
  fMquad[15] = fMFacOrd[1]*exp(-(0.5*(std::pow(wc[3]-1.0*flowUOrd[10],2.0)+std::pow(wc[2]+0.3872983346207417*dxv[2]-1.0*flowUOrd[1],2.0)))/vtSqOrd[1]); 
  fMquad[16] = fMFacOrd[1]*exp(-(0.5*(std::pow((-1.0*flowUOrd[10])+wc[3]-0.3872983346207417*dxv[3],2.0)+std::pow(wc[2]+0.3872983346207417*dxv[2]-1.0*flowUOrd[1],2.0)))/vtSqOrd[1]); 
  fMquad[17] = fMFacOrd[1]*exp(-(0.5*(std::pow((-1.0*flowUOrd[10])+wc[3]+0.3872983346207417*dxv[3],2.0)+std::pow(wc[2]+0.3872983346207417*dxv[2]-1.0*flowUOrd[1],2.0)))/vtSqOrd[1]); 
  fMquad[18] = fMFacOrd[2]*exp(-(0.5*(std::pow(wc[3]-1.0*flowUOrd[11],2.0)+std::pow(wc[2]-1.0*flowUOrd[2],2.0)))/vtSqOrd[2]); 
  fMquad[19] = fMFacOrd[2]*exp(-(0.5*(std::pow((-1.0*flowUOrd[11])+wc[3]-0.3872983346207417*dxv[3],2.0)+std::pow(wc[2]-1.0*flowUOrd[2],2.0)))/vtSqOrd[2]); 
  fMquad[20] = fMFacOrd[2]*exp(-(0.5*(std::pow((-1.0*flowUOrd[11])+wc[3]+0.3872983346207417*dxv[3],2.0)+std::pow(wc[2]-1.0*flowUOrd[2],2.0)))/vtSqOrd[2]); 
  fMquad[21] = fMFacOrd[2]*exp(-(0.5*(std::pow(wc[3]-1.0*flowUOrd[11],2.0)+std::pow(wc[2]-1.0*flowUOrd[2]-0.3872983346207417*dxv[2],2.0)))/vtSqOrd[2]); 
  fMquad[22] = fMFacOrd[2]*exp(-(0.5*(std::pow((-1.0*flowUOrd[11])+wc[3]-0.3872983346207417*dxv[3],2.0)+std::pow(wc[2]-1.0*flowUOrd[2]-0.3872983346207417*dxv[2],2.0)))/vtSqOrd[2]); 
  fMquad[23] = fMFacOrd[2]*exp(-(0.5*(std::pow((-1.0*flowUOrd[11])+wc[3]+0.3872983346207417*dxv[3],2.0)+std::pow(wc[2]-1.0*flowUOrd[2]-0.3872983346207417*dxv[2],2.0)))/vtSqOrd[2]); 
  fMquad[24] = fMFacOrd[2]*exp(-(0.5*(std::pow(wc[3]-1.0*flowUOrd[11],2.0)+std::pow(wc[2]-1.0*flowUOrd[2]+0.3872983346207417*dxv[2],2.0)))/vtSqOrd[2]); 
  fMquad[25] = fMFacOrd[2]*exp(-(0.5*(std::pow((-1.0*flowUOrd[11])+wc[3]-0.3872983346207417*dxv[3],2.0)+std::pow(wc[2]-1.0*flowUOrd[2]+0.3872983346207417*dxv[2],2.0)))/vtSqOrd[2]); 
  fMquad[26] = fMFacOrd[2]*exp(-(0.5*(std::pow((-1.0*flowUOrd[11])+wc[3]+0.3872983346207417*dxv[3],2.0)+std::pow(wc[2]-1.0*flowUOrd[2]+0.3872983346207417*dxv[2],2.0)))/vtSqOrd[2]); 
  fMquad[27] = fMFacOrd[3]*exp(-(0.5*(std::pow(wc[3]-1.0*flowUOrd[12],2.0)+std::pow(wc[2]-1.0*flowUOrd[3],2.0)))/vtSqOrd[3]); 
  fMquad[28] = fMFacOrd[3]*exp(-(0.5*(std::pow((-1.0*flowUOrd[12])+wc[3]-0.3872983346207417*dxv[3],2.0)+std::pow(wc[2]-1.0*flowUOrd[3],2.0)))/vtSqOrd[3]); 
  fMquad[29] = fMFacOrd[3]*exp(-(0.5*(std::pow((-1.0*flowUOrd[12])+wc[3]+0.3872983346207417*dxv[3],2.0)+std::pow(wc[2]-1.0*flowUOrd[3],2.0)))/vtSqOrd[3]); 
  fMquad[30] = fMFacOrd[3]*exp(-(0.5*(std::pow(wc[3]-1.0*flowUOrd[12],2.0)+std::pow((-1.0*flowUOrd[3])+wc[2]-0.3872983346207417*dxv[2],2.0)))/vtSqOrd[3]); 
  fMquad[31] = fMFacOrd[3]*exp(-(0.5*(std::pow((-1.0*flowUOrd[12])+wc[3]-0.3872983346207417*dxv[3],2.0)+std::pow((-1.0*flowUOrd[3])+wc[2]-0.3872983346207417*dxv[2],2.0)))/vtSqOrd[3]); 
  fMquad[32] = fMFacOrd[3]*exp(-(0.5*(std::pow((-1.0*flowUOrd[12])+wc[3]+0.3872983346207417*dxv[3],2.0)+std::pow((-1.0*flowUOrd[3])+wc[2]-0.3872983346207417*dxv[2],2.0)))/vtSqOrd[3]); 
  fMquad[33] = fMFacOrd[3]*exp(-(0.5*(std::pow(wc[3]-1.0*flowUOrd[12],2.0)+std::pow((-1.0*flowUOrd[3])+wc[2]+0.3872983346207417*dxv[2],2.0)))/vtSqOrd[3]); 
  fMquad[34] = fMFacOrd[3]*exp(-(0.5*(std::pow((-1.0*flowUOrd[12])+wc[3]-0.3872983346207417*dxv[3],2.0)+std::pow((-1.0*flowUOrd[3])+wc[2]+0.3872983346207417*dxv[2],2.0)))/vtSqOrd[3]); 
  fMquad[35] = fMFacOrd[3]*exp(-(0.5*(std::pow((-1.0*flowUOrd[12])+wc[3]+0.3872983346207417*dxv[3],2.0)+std::pow((-1.0*flowUOrd[3])+wc[2]+0.3872983346207417*dxv[2],2.0)))/vtSqOrd[3]); 
  fMquad[36] = fMFacOrd[4]*exp(-(0.5*(std::pow(wc[3]-1.0*flowUOrd[13],2.0)+std::pow(wc[2]-1.0*flowUOrd[4],2.0)))/vtSqOrd[4]); 
  fMquad[37] = fMFacOrd[4]*exp(-(0.5*(std::pow((-1.0*flowUOrd[13])+wc[3]-0.3872983346207417*dxv[3],2.0)+std::pow(wc[2]-1.0*flowUOrd[4],2.0)))/vtSqOrd[4]); 
  fMquad[38] = fMFacOrd[4]*exp(-(0.5*(std::pow((-1.0*flowUOrd[13])+wc[3]+0.3872983346207417*dxv[3],2.0)+std::pow(wc[2]-1.0*flowUOrd[4],2.0)))/vtSqOrd[4]); 
  fMquad[39] = fMFacOrd[4]*exp(-(0.5*(std::pow(wc[3]-1.0*flowUOrd[13],2.0)+std::pow((-1.0*flowUOrd[4])+wc[2]-0.3872983346207417*dxv[2],2.0)))/vtSqOrd[4]); 
  fMquad[40] = fMFacOrd[4]*exp(-(0.5*(std::pow((-1.0*flowUOrd[13])+wc[3]-0.3872983346207417*dxv[3],2.0)+std::pow((-1.0*flowUOrd[4])+wc[2]-0.3872983346207417*dxv[2],2.0)))/vtSqOrd[4]); 
  fMquad[41] = fMFacOrd[4]*exp(-(0.5*(std::pow((-1.0*flowUOrd[13])+wc[3]+0.3872983346207417*dxv[3],2.0)+std::pow((-1.0*flowUOrd[4])+wc[2]-0.3872983346207417*dxv[2],2.0)))/vtSqOrd[4]); 
  fMquad[42] = fMFacOrd[4]*exp(-(0.5*(std::pow(wc[3]-1.0*flowUOrd[13],2.0)+std::pow((-1.0*flowUOrd[4])+wc[2]+0.3872983346207417*dxv[2],2.0)))/vtSqOrd[4]); 
  fMquad[43] = fMFacOrd[4]*exp(-(0.5*(std::pow((-1.0*flowUOrd[13])+wc[3]-0.3872983346207417*dxv[3],2.0)+std::pow((-1.0*flowUOrd[4])+wc[2]+0.3872983346207417*dxv[2],2.0)))/vtSqOrd[4]); 
  fMquad[44] = fMFacOrd[4]*exp(-(0.5*(std::pow((-1.0*flowUOrd[13])+wc[3]+0.3872983346207417*dxv[3],2.0)+std::pow((-1.0*flowUOrd[4])+wc[2]+0.3872983346207417*dxv[2],2.0)))/vtSqOrd[4]); 
  fMquad[45] = fMFacOrd[5]*exp(-(0.5*(std::pow(wc[3]-1.0*flowUOrd[14],2.0)+std::pow(wc[2]-1.0*flowUOrd[5],2.0)))/vtSqOrd[5]); 
  fMquad[46] = fMFacOrd[5]*exp(-(0.5*(std::pow((-1.0*flowUOrd[14])+wc[3]-0.3872983346207417*dxv[3],2.0)+std::pow(wc[2]-1.0*flowUOrd[5],2.0)))/vtSqOrd[5]); 
  fMquad[47] = fMFacOrd[5]*exp(-(0.5*(std::pow((-1.0*flowUOrd[14])+wc[3]+0.3872983346207417*dxv[3],2.0)+std::pow(wc[2]-1.0*flowUOrd[5],2.0)))/vtSqOrd[5]); 
  fMquad[48] = fMFacOrd[5]*exp(-(0.5*(std::pow(wc[3]-1.0*flowUOrd[14],2.0)+std::pow((-1.0*flowUOrd[5])+wc[2]-0.3872983346207417*dxv[2],2.0)))/vtSqOrd[5]); 
  fMquad[49] = fMFacOrd[5]*exp(-(0.5*(std::pow((-1.0*flowUOrd[14])+wc[3]-0.3872983346207417*dxv[3],2.0)+std::pow((-1.0*flowUOrd[5])+wc[2]-0.3872983346207417*dxv[2],2.0)))/vtSqOrd[5]); 
  fMquad[50] = fMFacOrd[5]*exp(-(0.5*(std::pow((-1.0*flowUOrd[14])+wc[3]+0.3872983346207417*dxv[3],2.0)+std::pow((-1.0*flowUOrd[5])+wc[2]-0.3872983346207417*dxv[2],2.0)))/vtSqOrd[5]); 
  fMquad[51] = fMFacOrd[5]*exp(-(0.5*(std::pow(wc[3]-1.0*flowUOrd[14],2.0)+std::pow((-1.0*flowUOrd[5])+wc[2]+0.3872983346207417*dxv[2],2.0)))/vtSqOrd[5]); 
  fMquad[52] = fMFacOrd[5]*exp(-(0.5*(std::pow((-1.0*flowUOrd[14])+wc[3]-0.3872983346207417*dxv[3],2.0)+std::pow((-1.0*flowUOrd[5])+wc[2]+0.3872983346207417*dxv[2],2.0)))/vtSqOrd[5]); 
  fMquad[53] = fMFacOrd[5]*exp(-(0.5*(std::pow((-1.0*flowUOrd[14])+wc[3]+0.3872983346207417*dxv[3],2.0)+std::pow((-1.0*flowUOrd[5])+wc[2]+0.3872983346207417*dxv[2],2.0)))/vtSqOrd[5]); 
  fMquad[54] = fMFacOrd[6]*exp(-(0.5*(std::pow(wc[3]-1.0*flowUOrd[15],2.0)+std::pow(wc[2]-1.0*flowUOrd[6],2.0)))/vtSqOrd[6]); 
  fMquad[55] = fMFacOrd[6]*exp(-(0.5*(std::pow((-1.0*flowUOrd[15])+wc[3]-0.3872983346207417*dxv[3],2.0)+std::pow(wc[2]-1.0*flowUOrd[6],2.0)))/vtSqOrd[6]); 
  fMquad[56] = fMFacOrd[6]*exp(-(0.5*(std::pow((-1.0*flowUOrd[15])+wc[3]+0.3872983346207417*dxv[3],2.0)+std::pow(wc[2]-1.0*flowUOrd[6],2.0)))/vtSqOrd[6]); 
  fMquad[57] = fMFacOrd[6]*exp(-(0.5*(std::pow(wc[3]-1.0*flowUOrd[15],2.0)+std::pow((-1.0*flowUOrd[6])+wc[2]-0.3872983346207417*dxv[2],2.0)))/vtSqOrd[6]); 
  fMquad[58] = fMFacOrd[6]*exp(-(0.5*(std::pow((-1.0*flowUOrd[15])+wc[3]-0.3872983346207417*dxv[3],2.0)+std::pow((-1.0*flowUOrd[6])+wc[2]-0.3872983346207417*dxv[2],2.0)))/vtSqOrd[6]); 
  fMquad[59] = fMFacOrd[6]*exp(-(0.5*(std::pow((-1.0*flowUOrd[15])+wc[3]+0.3872983346207417*dxv[3],2.0)+std::pow((-1.0*flowUOrd[6])+wc[2]-0.3872983346207417*dxv[2],2.0)))/vtSqOrd[6]); 
  fMquad[60] = fMFacOrd[6]*exp(-(0.5*(std::pow(wc[3]-1.0*flowUOrd[15],2.0)+std::pow((-1.0*flowUOrd[6])+wc[2]+0.3872983346207417*dxv[2],2.0)))/vtSqOrd[6]); 
  fMquad[61] = fMFacOrd[6]*exp(-(0.5*(std::pow((-1.0*flowUOrd[15])+wc[3]-0.3872983346207417*dxv[3],2.0)+std::pow((-1.0*flowUOrd[6])+wc[2]+0.3872983346207417*dxv[2],2.0)))/vtSqOrd[6]); 
  fMquad[62] = fMFacOrd[6]*exp(-(0.5*(std::pow((-1.0*flowUOrd[15])+wc[3]+0.3872983346207417*dxv[3],2.0)+std::pow((-1.0*flowUOrd[6])+wc[2]+0.3872983346207417*dxv[2],2.0)))/vtSqOrd[6]); 
  fMquad[63] = fMFacOrd[7]*exp(-(0.5*(std::pow(wc[3]-1.0*flowUOrd[16],2.0)+std::pow(wc[2]-1.0*flowUOrd[7],2.0)))/vtSqOrd[7]); 
  fMquad[64] = fMFacOrd[7]*exp(-(0.5*(std::pow((-1.0*flowUOrd[16])+wc[3]-0.3872983346207417*dxv[3],2.0)+std::pow(wc[2]-1.0*flowUOrd[7],2.0)))/vtSqOrd[7]); 
  fMquad[65] = fMFacOrd[7]*exp(-(0.5*(std::pow((-1.0*flowUOrd[16])+wc[3]+0.3872983346207417*dxv[3],2.0)+std::pow(wc[2]-1.0*flowUOrd[7],2.0)))/vtSqOrd[7]); 
  fMquad[66] = fMFacOrd[7]*exp(-(0.5*(std::pow(wc[3]-1.0*flowUOrd[16],2.0)+std::pow((-1.0*flowUOrd[7])+wc[2]-0.3872983346207417*dxv[2],2.0)))/vtSqOrd[7]); 
  fMquad[67] = fMFacOrd[7]*exp(-(0.5*(std::pow((-1.0*flowUOrd[16])+wc[3]-0.3872983346207417*dxv[3],2.0)+std::pow((-1.0*flowUOrd[7])+wc[2]-0.3872983346207417*dxv[2],2.0)))/vtSqOrd[7]); 
  fMquad[68] = fMFacOrd[7]*exp(-(0.5*(std::pow((-1.0*flowUOrd[16])+wc[3]+0.3872983346207417*dxv[3],2.0)+std::pow((-1.0*flowUOrd[7])+wc[2]-0.3872983346207417*dxv[2],2.0)))/vtSqOrd[7]); 
  fMquad[69] = fMFacOrd[7]*exp(-(0.5*(std::pow(wc[3]-1.0*flowUOrd[16],2.0)+std::pow((-1.0*flowUOrd[7])+wc[2]+0.3872983346207417*dxv[2],2.0)))/vtSqOrd[7]); 
  fMquad[70] = fMFacOrd[7]*exp(-(0.5*(std::pow((-1.0*flowUOrd[16])+wc[3]-0.3872983346207417*dxv[3],2.0)+std::pow((-1.0*flowUOrd[7])+wc[2]+0.3872983346207417*dxv[2],2.0)))/vtSqOrd[7]); 
  fMquad[71] = fMFacOrd[7]*exp(-(0.5*(std::pow((-1.0*flowUOrd[16])+wc[3]+0.3872983346207417*dxv[3],2.0)+std::pow((-1.0*flowUOrd[7])+wc[2]+0.3872983346207417*dxv[2],2.0)))/vtSqOrd[7]); 
  fMquad[72] = fMFacOrd[8]*exp(-(0.5*(std::pow(wc[3]-1.0*flowUOrd[17],2.0)+std::pow(wc[2]-1.0*flowUOrd[8],2.0)))/vtSqOrd[8]); 
  fMquad[73] = fMFacOrd[8]*exp(-(0.5*(std::pow((-1.0*flowUOrd[17])+wc[3]-0.3872983346207417*dxv[3],2.0)+std::pow(wc[2]-1.0*flowUOrd[8],2.0)))/vtSqOrd[8]); 
  fMquad[74] = fMFacOrd[8]*exp(-(0.5*(std::pow((-1.0*flowUOrd[17])+wc[3]+0.3872983346207417*dxv[3],2.0)+std::pow(wc[2]-1.0*flowUOrd[8],2.0)))/vtSqOrd[8]); 
  fMquad[75] = fMFacOrd[8]*exp(-(0.5*(std::pow(wc[3]-1.0*flowUOrd[17],2.0)+std::pow((-1.0*flowUOrd[8])+wc[2]-0.3872983346207417*dxv[2],2.0)))/vtSqOrd[8]); 
  fMquad[76] = fMFacOrd[8]*exp(-(0.5*(std::pow((-1.0*flowUOrd[17])+wc[3]-0.3872983346207417*dxv[3],2.0)+std::pow((-1.0*flowUOrd[8])+wc[2]-0.3872983346207417*dxv[2],2.0)))/vtSqOrd[8]); 
  fMquad[77] = fMFacOrd[8]*exp(-(0.5*(std::pow((-1.0*flowUOrd[17])+wc[3]+0.3872983346207417*dxv[3],2.0)+std::pow((-1.0*flowUOrd[8])+wc[2]-0.3872983346207417*dxv[2],2.0)))/vtSqOrd[8]); 
  fMquad[78] = fMFacOrd[8]*exp(-(0.5*(std::pow(wc[3]-1.0*flowUOrd[17],2.0)+std::pow((-1.0*flowUOrd[8])+wc[2]+0.3872983346207417*dxv[2],2.0)))/vtSqOrd[8]); 
  fMquad[79] = fMFacOrd[8]*exp(-(0.5*(std::pow((-1.0*flowUOrd[17])+wc[3]-0.3872983346207417*dxv[3],2.0)+std::pow((-1.0*flowUOrd[8])+wc[2]+0.3872983346207417*dxv[2],2.0)))/vtSqOrd[8]); 
  fMquad[80] = fMFacOrd[8]*exp(-(0.5*(std::pow((-1.0*flowUOrd[17])+wc[3]+0.3872983346207417*dxv[3],2.0)+std::pow((-1.0*flowUOrd[8])+wc[2]+0.3872983346207417*dxv[2],2.0)))/vtSqOrd[8]); 

  fMOut[0] = 0.0771604938271605*(fMquad[80]+fMquad[79]+fMquad[78]+fMquad[77]+fMquad[76]+fMquad[75]+fMquad[74]+fMquad[73]+fMquad[72]+fMquad[71]+fMquad[70]+fMquad[69]+fMquad[68]+fMquad[67]+fMquad[66]+fMquad[65]+fMquad[64]+fMquad[63])+0.1234567901234568*(fMquad[62]+fMquad[61]+fMquad[60]+fMquad[59]+fMquad[58]+fMquad[57]+fMquad[56]+fMquad[55]+fMquad[54])+0.0771604938271605*(fMquad[53]+fMquad[52]+fMquad[51]+fMquad[50]+fMquad[49]+fMquad[48]+fMquad[47]+fMquad[46]+fMquad[45]+fMquad[44]+fMquad[43]+fMquad[42]+fMquad[41]+fMquad[40]+fMquad[39]+fMquad[38]+fMquad[37]+fMquad[36])+0.1234567901234568*(fMquad[35]+fMquad[34]+fMquad[33]+fMquad[32]+fMquad[31]+fMquad[30]+fMquad[29]+fMquad[28]+fMquad[27]+fMquad[26]+fMquad[25]+fMquad[24]+fMquad[23]+fMquad[22]+fMquad[21]+fMquad[20]+fMquad[19]+fMquad[18]+fMquad[17]+fMquad[16]+fMquad[15]+fMquad[14]+fMquad[13]+fMquad[12]+fMquad[11]+fMquad[10]+fMquad[9])+0.1975308641975309*(fMquad[8]+fMquad[7]+fMquad[6]+fMquad[5]+fMquad[4]+fMquad[3]+fMquad[2]+fMquad[1]+fMquad[0]); 
  fMOut[1] = 0.1035216656249903*(fMquad[80]+fMquad[79]+fMquad[78]+fMquad[77]+fMquad[76]+fMquad[75]+fMquad[74]+fMquad[73]+fMquad[72]+fMquad[71]+fMquad[70]+fMquad[69]+fMquad[68]+fMquad[67]+fMquad[66]+fMquad[65]+fMquad[64]+fMquad[63])+0.1656346649999844*(fMquad[62]+fMquad[61]+fMquad[60]+fMquad[59]+fMquad[58]+fMquad[57]+fMquad[56]+fMquad[55]+fMquad[54])-0.1035216656249903*(fMquad[53]+fMquad[52]+fMquad[51]+fMquad[50]+fMquad[49]+fMquad[48]+fMquad[47]+fMquad[46]+fMquad[45]+fMquad[44]+fMquad[43]+fMquad[42]+fMquad[41]+fMquad[40]+fMquad[39]+fMquad[38]+fMquad[37]+fMquad[36])-0.1656346649999844*(fMquad[35]+fMquad[34]+fMquad[33]+fMquad[32]+fMquad[31]+fMquad[30]+fMquad[29]+fMquad[28]+fMquad[27]); 
  fMOut[2] = 0.1035216656249903*(fMquad[80]+fMquad[79]+fMquad[78]+fMquad[77]+fMquad[76]+fMquad[75]+fMquad[74]+fMquad[73]+fMquad[72])-0.1035216656249903*(fMquad[71]+fMquad[70]+fMquad[69]+fMquad[68]+fMquad[67]+fMquad[66]+fMquad[65]+fMquad[64]+fMquad[63])+0.1035216656249903*(fMquad[53]+fMquad[52]+fMquad[51]+fMquad[50]+fMquad[49]+fMquad[48]+fMquad[47]+fMquad[46]+fMquad[45])-0.1035216656249903*(fMquad[44]+fMquad[43]+fMquad[42]+fMquad[41]+fMquad[40]+fMquad[39]+fMquad[38]+fMquad[37]+fMquad[36])+0.1656346649999844*(fMquad[26]+fMquad[25]+fMquad[24]+fMquad[23]+fMquad[22]+fMquad[21]+fMquad[20]+fMquad[19]+fMquad[18])-0.1656346649999844*(fMquad[17]+fMquad[16]+fMquad[15]+fMquad[14]+fMquad[13]+fMquad[12]+fMquad[11]+fMquad[10]+fMquad[9]); 
  fMOut[3] = 0.1035216656249903*(fMquad[80]+fMquad[79]+fMquad[78])-0.1035216656249903*(fMquad[77]+fMquad[76]+fMquad[75])+0.1035216656249903*(fMquad[71]+fMquad[70]+fMquad[69])-0.1035216656249903*(fMquad[68]+fMquad[67]+fMquad[66])+0.1656346649999844*(fMquad[62]+fMquad[61]+fMquad[60])-0.1656346649999844*(fMquad[59]+fMquad[58]+fMquad[57])+0.1035216656249903*(fMquad[53]+fMquad[52]+fMquad[51])-0.1035216656249903*(fMquad[50]+fMquad[49]+fMquad[48])+0.1035216656249903*(fMquad[44]+fMquad[43]+fMquad[42])-0.1035216656249903*(fMquad[41]+fMquad[40]+fMquad[39])+0.1656346649999844*(fMquad[35]+fMquad[34]+fMquad[33])-0.1656346649999844*(fMquad[32]+fMquad[31]+fMquad[30])+0.1656346649999844*(fMquad[26]+fMquad[25]+fMquad[24])-0.1656346649999844*(fMquad[23]+fMquad[22]+fMquad[21])+0.1656346649999844*(fMquad[17]+fMquad[16]+fMquad[15])-0.1656346649999844*(fMquad[14]+fMquad[13]+fMquad[12])+0.2650154639999751*(fMquad[8]+fMquad[7]+fMquad[6])-0.2650154639999751*(fMquad[5]+fMquad[4]+fMquad[3]); 
  fMOut[4] = 0.1035216656249903*fMquad[80]-0.1035216656249903*fMquad[79]+0.1035216656249903*fMquad[77]-0.1035216656249903*fMquad[76]+0.1035216656249903*fMquad[74]-0.1035216656249903*fMquad[73]+0.1035216656249903*fMquad[71]-0.1035216656249903*fMquad[70]+0.1035216656249903*fMquad[68]-0.1035216656249903*fMquad[67]+0.1035216656249903*fMquad[65]-0.1035216656249903*fMquad[64]+0.1656346649999844*fMquad[62]-0.1656346649999844*fMquad[61]+0.1656346649999844*fMquad[59]-0.1656346649999844*fMquad[58]+0.1656346649999844*fMquad[56]-0.1656346649999844*fMquad[55]+0.1035216656249903*fMquad[53]-0.1035216656249903*fMquad[52]+0.1035216656249903*fMquad[50]-0.1035216656249903*fMquad[49]+0.1035216656249903*fMquad[47]-0.1035216656249903*fMquad[46]+0.1035216656249903*fMquad[44]-0.1035216656249903*fMquad[43]+0.1035216656249903*fMquad[41]-0.1035216656249903*fMquad[40]+0.1035216656249903*fMquad[38]-0.1035216656249903*fMquad[37]+0.1656346649999844*fMquad[35]-0.1656346649999844*fMquad[34]+0.1656346649999844*fMquad[32]-0.1656346649999844*fMquad[31]+0.1656346649999844*fMquad[29]-0.1656346649999844*fMquad[28]+0.1656346649999844*fMquad[26]-0.1656346649999844*fMquad[25]+0.1656346649999844*fMquad[23]-0.1656346649999844*fMquad[22]+0.1656346649999844*fMquad[20]-0.1656346649999844*fMquad[19]+0.1656346649999844*fMquad[17]-0.1656346649999844*fMquad[16]+0.1656346649999844*fMquad[14]-0.1656346649999844*fMquad[13]+0.1656346649999844*fMquad[11]-0.1656346649999844*fMquad[10]+0.2650154639999751*fMquad[8]-0.2650154639999751*fMquad[7]+0.2650154639999751*fMquad[5]-0.2650154639999751*fMquad[4]+0.2650154639999751*fMquad[2]-0.2650154639999751*fMquad[1]; 
  fMOut[5] = 0.1388888888888889*(fMquad[80]+fMquad[79]+fMquad[78]+fMquad[77]+fMquad[76]+fMquad[75]+fMquad[74]+fMquad[73]+fMquad[72])-0.1388888888888889*(fMquad[71]+fMquad[70]+fMquad[69]+fMquad[68]+fMquad[67]+fMquad[66]+fMquad[65]+fMquad[64]+fMquad[63]+fMquad[53]+fMquad[52]+fMquad[51]+fMquad[50]+fMquad[49]+fMquad[48]+fMquad[47]+fMquad[46]+fMquad[45])+0.1388888888888889*(fMquad[44]+fMquad[43]+fMquad[42]+fMquad[41]+fMquad[40]+fMquad[39]+fMquad[38]+fMquad[37]+fMquad[36]); 
  fMOut[6] = 0.1388888888888889*(fMquad[80]+fMquad[79]+fMquad[78])-0.1388888888888889*(fMquad[77]+fMquad[76]+fMquad[75])+0.1388888888888889*(fMquad[71]+fMquad[70]+fMquad[69])-0.1388888888888889*(fMquad[68]+fMquad[67]+fMquad[66])+0.2222222222222222*(fMquad[62]+fMquad[61]+fMquad[60])-0.2222222222222222*(fMquad[59]+fMquad[58]+fMquad[57])-0.1388888888888889*(fMquad[53]+fMquad[52]+fMquad[51])+0.1388888888888889*(fMquad[50]+fMquad[49]+fMquad[48])-0.1388888888888889*(fMquad[44]+fMquad[43]+fMquad[42])+0.1388888888888889*(fMquad[41]+fMquad[40]+fMquad[39])-0.2222222222222222*(fMquad[35]+fMquad[34]+fMquad[33])+0.2222222222222222*(fMquad[32]+fMquad[31]+fMquad[30]); 
  fMOut[7] = 0.1388888888888889*(fMquad[80]+fMquad[79]+fMquad[78])-0.1388888888888889*(fMquad[77]+fMquad[76]+fMquad[75]+fMquad[71]+fMquad[70]+fMquad[69])+0.1388888888888889*(fMquad[68]+fMquad[67]+fMquad[66]+fMquad[53]+fMquad[52]+fMquad[51])-0.1388888888888889*(fMquad[50]+fMquad[49]+fMquad[48]+fMquad[44]+fMquad[43]+fMquad[42])+0.1388888888888889*(fMquad[41]+fMquad[40]+fMquad[39])+0.2222222222222222*(fMquad[26]+fMquad[25]+fMquad[24])-0.2222222222222222*(fMquad[23]+fMquad[22]+fMquad[21]+fMquad[17]+fMquad[16]+fMquad[15])+0.2222222222222222*(fMquad[14]+fMquad[13]+fMquad[12]); 
  fMOut[8] = 0.1388888888888889*fMquad[80]-0.1388888888888889*fMquad[79]+0.1388888888888889*fMquad[77]-0.1388888888888889*fMquad[76]+0.1388888888888889*fMquad[74]-0.1388888888888889*fMquad[73]+0.1388888888888889*fMquad[71]-0.1388888888888889*fMquad[70]+0.1388888888888889*fMquad[68]-0.1388888888888889*fMquad[67]+0.1388888888888889*fMquad[65]-0.1388888888888889*fMquad[64]+0.2222222222222222*fMquad[62]-0.2222222222222222*fMquad[61]+0.2222222222222222*fMquad[59]-0.2222222222222222*fMquad[58]+0.2222222222222222*fMquad[56]-0.2222222222222222*fMquad[55]-0.1388888888888889*fMquad[53]+0.1388888888888889*fMquad[52]-0.1388888888888889*fMquad[50]+0.1388888888888889*fMquad[49]-0.1388888888888889*fMquad[47]+0.1388888888888889*fMquad[46]-0.1388888888888889*fMquad[44]+0.1388888888888889*fMquad[43]-0.1388888888888889*fMquad[41]+0.1388888888888889*fMquad[40]-0.1388888888888889*fMquad[38]+0.1388888888888889*fMquad[37]-0.2222222222222222*fMquad[35]+0.2222222222222222*fMquad[34]-0.2222222222222222*fMquad[32]+0.2222222222222222*fMquad[31]-0.2222222222222222*fMquad[29]+0.2222222222222222*fMquad[28]; 
  fMOut[9] = 0.1388888888888889*fMquad[80]-0.1388888888888889*fMquad[79]+0.1388888888888889*fMquad[77]-0.1388888888888889*fMquad[76]+0.1388888888888889*fMquad[74]-0.1388888888888889*(fMquad[73]+fMquad[71])+0.1388888888888889*fMquad[70]-0.1388888888888889*fMquad[68]+0.1388888888888889*fMquad[67]-0.1388888888888889*fMquad[65]+0.1388888888888889*(fMquad[64]+fMquad[53])-0.1388888888888889*fMquad[52]+0.1388888888888889*fMquad[50]-0.1388888888888889*fMquad[49]+0.1388888888888889*fMquad[47]-0.1388888888888889*(fMquad[46]+fMquad[44])+0.1388888888888889*fMquad[43]-0.1388888888888889*fMquad[41]+0.1388888888888889*fMquad[40]-0.1388888888888889*fMquad[38]+0.1388888888888889*fMquad[37]+0.2222222222222222*fMquad[26]-0.2222222222222222*fMquad[25]+0.2222222222222222*fMquad[23]-0.2222222222222222*fMquad[22]+0.2222222222222222*fMquad[20]-0.2222222222222222*(fMquad[19]+fMquad[17])+0.2222222222222222*fMquad[16]-0.2222222222222222*fMquad[14]+0.2222222222222222*fMquad[13]-0.2222222222222222*fMquad[11]+0.2222222222222222*fMquad[10]; 
  fMOut[10] = 0.1388888888888889*fMquad[80]-0.1388888888888889*(fMquad[79]+fMquad[77])+0.1388888888888889*(fMquad[76]+fMquad[71])-0.1388888888888889*(fMquad[70]+fMquad[68])+0.1388888888888889*fMquad[67]+0.2222222222222222*fMquad[62]-0.2222222222222222*(fMquad[61]+fMquad[59])+0.2222222222222222*fMquad[58]+0.1388888888888889*fMquad[53]-0.1388888888888889*(fMquad[52]+fMquad[50])+0.1388888888888889*(fMquad[49]+fMquad[44])-0.1388888888888889*(fMquad[43]+fMquad[41])+0.1388888888888889*fMquad[40]+0.2222222222222222*fMquad[35]-0.2222222222222222*(fMquad[34]+fMquad[32])+0.2222222222222222*(fMquad[31]+fMquad[26])-0.2222222222222222*(fMquad[25]+fMquad[23])+0.2222222222222222*(fMquad[22]+fMquad[17])-0.2222222222222222*(fMquad[16]+fMquad[14])+0.2222222222222222*fMquad[13]+0.3555555555555556*fMquad[8]-0.3555555555555556*(fMquad[7]+fMquad[5])+0.3555555555555556*fMquad[4]; 
  fMOut[11] = 0.06901444374999355*(fMquad[80]+fMquad[79]+fMquad[78]+fMquad[77]+fMquad[76]+fMquad[75]+fMquad[74]+fMquad[73]+fMquad[72]+fMquad[71]+fMquad[70]+fMquad[69]+fMquad[68]+fMquad[67]+fMquad[66]+fMquad[65]+fMquad[64]+fMquad[63])+0.1104231099999896*(fMquad[62]+fMquad[61]+fMquad[60]+fMquad[59]+fMquad[58]+fMquad[57]+fMquad[56]+fMquad[55]+fMquad[54])+0.06901444374999355*(fMquad[53]+fMquad[52]+fMquad[51]+fMquad[50]+fMquad[49]+fMquad[48]+fMquad[47]+fMquad[46]+fMquad[45]+fMquad[44]+fMquad[43]+fMquad[42]+fMquad[41]+fMquad[40]+fMquad[39]+fMquad[38]+fMquad[37]+fMquad[36])+0.1104231099999896*(fMquad[35]+fMquad[34]+fMquad[33]+fMquad[32]+fMquad[31]+fMquad[30]+fMquad[29]+fMquad[28]+fMquad[27])-0.138028887499987*(fMquad[26]+fMquad[25]+fMquad[24]+fMquad[23]+fMquad[22]+fMquad[21]+fMquad[20]+fMquad[19]+fMquad[18]+fMquad[17]+fMquad[16]+fMquad[15]+fMquad[14]+fMquad[13]+fMquad[12]+fMquad[11]+fMquad[10]+fMquad[9])-0.2208462199999792*(fMquad[8]+fMquad[7]+fMquad[6]+fMquad[5]+fMquad[4]+fMquad[3]+fMquad[2]+fMquad[1]+fMquad[0]); 
  fMOut[12] = 0.06901444374999355*(fMquad[80]+fMquad[79]+fMquad[78]+fMquad[77]+fMquad[76]+fMquad[75]+fMquad[74]+fMquad[73]+fMquad[72]+fMquad[71]+fMquad[70]+fMquad[69]+fMquad[68]+fMquad[67]+fMquad[66]+fMquad[65]+fMquad[64]+fMquad[63])-0.138028887499987*(fMquad[62]+fMquad[61]+fMquad[60]+fMquad[59]+fMquad[58]+fMquad[57]+fMquad[56]+fMquad[55]+fMquad[54])+0.06901444374999355*(fMquad[53]+fMquad[52]+fMquad[51]+fMquad[50]+fMquad[49]+fMquad[48]+fMquad[47]+fMquad[46]+fMquad[45]+fMquad[44]+fMquad[43]+fMquad[42]+fMquad[41]+fMquad[40]+fMquad[39]+fMquad[38]+fMquad[37]+fMquad[36])-0.138028887499987*(fMquad[35]+fMquad[34]+fMquad[33]+fMquad[32]+fMquad[31]+fMquad[30]+fMquad[29]+fMquad[28]+fMquad[27])+0.1104231099999896*(fMquad[26]+fMquad[25]+fMquad[24]+fMquad[23]+fMquad[22]+fMquad[21]+fMquad[20]+fMquad[19]+fMquad[18]+fMquad[17]+fMquad[16]+fMquad[15]+fMquad[14]+fMquad[13]+fMquad[12]+fMquad[11]+fMquad[10]+fMquad[9])-0.2208462199999792*(fMquad[8]+fMquad[7]+fMquad[6]+fMquad[5]+fMquad[4]+fMquad[3]+fMquad[2]+fMquad[1]+fMquad[0]); 
  fMOut[13] = 0.06901444374999355*(fMquad[80]+fMquad[79]+fMquad[78]+fMquad[77]+fMquad[76]+fMquad[75])-0.08626805468749191*(fMquad[74]+fMquad[73]+fMquad[72])+0.06901444374999355*(fMquad[71]+fMquad[70]+fMquad[69]+fMquad[68]+fMquad[67]+fMquad[66])-0.08626805468749191*(fMquad[65]+fMquad[64]+fMquad[63])+0.1104231099999896*(fMquad[62]+fMquad[61]+fMquad[60]+fMquad[59]+fMquad[58]+fMquad[57])-0.138028887499987*(fMquad[56]+fMquad[55]+fMquad[54])+0.06901444374999355*(fMquad[53]+fMquad[52]+fMquad[51]+fMquad[50]+fMquad[49]+fMquad[48])-0.08626805468749191*(fMquad[47]+fMquad[46]+fMquad[45])+0.06901444374999355*(fMquad[44]+fMquad[43]+fMquad[42]+fMquad[41]+fMquad[40]+fMquad[39])-0.08626805468749191*(fMquad[38]+fMquad[37]+fMquad[36])+0.1104231099999896*(fMquad[35]+fMquad[34]+fMquad[33]+fMquad[32]+fMquad[31]+fMquad[30])-0.138028887499987*(fMquad[29]+fMquad[28]+fMquad[27])+0.1104231099999896*(fMquad[26]+fMquad[25]+fMquad[24]+fMquad[23]+fMquad[22]+fMquad[21])-0.138028887499987*(fMquad[20]+fMquad[19]+fMquad[18])+0.1104231099999896*(fMquad[17]+fMquad[16]+fMquad[15]+fMquad[14]+fMquad[13]+fMquad[12])-0.138028887499987*(fMquad[11]+fMquad[10]+fMquad[9])+0.1766769759999834*(fMquad[8]+fMquad[7]+fMquad[6]+fMquad[5]+fMquad[4]+fMquad[3])-0.2208462199999792*(fMquad[2]+fMquad[1]+fMquad[0]); 
  fMOut[14] = 0.06901444374999355*(fMquad[80]+fMquad[79])-0.08626805468749191*fMquad[78]+0.06901444374999355*(fMquad[77]+fMquad[76])-0.08626805468749191*fMquad[75]+0.06901444374999355*(fMquad[74]+fMquad[73])-0.08626805468749191*fMquad[72]+0.06901444374999355*(fMquad[71]+fMquad[70])-0.08626805468749191*fMquad[69]+0.06901444374999355*(fMquad[68]+fMquad[67])-0.08626805468749191*fMquad[66]+0.06901444374999355*(fMquad[65]+fMquad[64])-0.08626805468749191*fMquad[63]+0.1104231099999896*(fMquad[62]+fMquad[61])-0.138028887499987*fMquad[60]+0.1104231099999896*(fMquad[59]+fMquad[58])-0.138028887499987*fMquad[57]+0.1104231099999896*(fMquad[56]+fMquad[55])-0.138028887499987*fMquad[54]+0.06901444374999355*(fMquad[53]+fMquad[52])-0.08626805468749191*fMquad[51]+0.06901444374999355*(fMquad[50]+fMquad[49])-0.08626805468749191*fMquad[48]+0.06901444374999355*(fMquad[47]+fMquad[46])-0.08626805468749191*fMquad[45]+0.06901444374999355*(fMquad[44]+fMquad[43])-0.08626805468749191*fMquad[42]+0.06901444374999355*(fMquad[41]+fMquad[40])-0.08626805468749191*fMquad[39]+0.06901444374999355*(fMquad[38]+fMquad[37])-0.08626805468749191*fMquad[36]+0.1104231099999896*(fMquad[35]+fMquad[34])-0.138028887499987*fMquad[33]+0.1104231099999896*(fMquad[32]+fMquad[31])-0.138028887499987*fMquad[30]+0.1104231099999896*(fMquad[29]+fMquad[28])-0.138028887499987*fMquad[27]+0.1104231099999896*(fMquad[26]+fMquad[25])-0.138028887499987*fMquad[24]+0.1104231099999896*(fMquad[23]+fMquad[22])-0.138028887499987*fMquad[21]+0.1104231099999896*(fMquad[20]+fMquad[19])-0.138028887499987*fMquad[18]+0.1104231099999896*(fMquad[17]+fMquad[16])-0.138028887499987*fMquad[15]+0.1104231099999896*(fMquad[14]+fMquad[13])-0.138028887499987*fMquad[12]+0.1104231099999896*(fMquad[11]+fMquad[10])-0.138028887499987*fMquad[9]+0.1766769759999834*(fMquad[8]+fMquad[7])-0.2208462199999792*fMquad[6]+0.1766769759999834*(fMquad[5]+fMquad[4])-0.2208462199999792*fMquad[3]+0.1766769759999834*(fMquad[2]+fMquad[1])-0.2208462199999792*fMquad[0]; 

}
