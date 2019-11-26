#include <GkLBOModDecl.h> 
double GkLBOconstNuSurfPositivity1x2vSer_Vpar_P1(const double m_, const double *cflRateByDirL, const double *cflRateByDirR, const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double dtApprox, const double *BmagInv, const double nuSum, const double vMuMidMax, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // m_:              species mass. 
  // cflRateByDir[3]: CFL rate in each direction. 
  // w[3]:            Cell-center coordinates. 
  // dxv[3]:          Cell spacing. 
  // nuSum:           collisionalities added (self and cross species collisionalities). 
  // vMuMidMax:       maximum midpoint value of v-u. 
  // nuUSum[2]:       sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[2]:    sum of thermal speeds squared time their respective collisionalities. 
  // fl/fr:           Distribution function in left/right cells 
  // outl/outr:       Incremented distribution function in left/right cells 
  double rdv = 1.0/dxvl[1]; 
  double rdv2L = 2.0/dxvl[1]; 
  double rdv2R = 2.0/dxvr[1]; 
  double rdvSq4L = rdv2L*rdv2L; 
  double rdvSq4R = rdv2R*rdv2R; 

  double alphaDrSurf[4]; 
  alphaDrSurf[0] = (2.0*wl[1]+dxvl[1])*nuSum-1.414213562373095*nuUSum[0]; 
  alphaDrSurf[1] = -1.414213562373095*nuUSum[1]; 

  double rCtrlL[4], rCtrlR[4];  // rCtrl=f1/f0 at each control node in dimensions other than vx 
  rCtrlL[0] = (1.414213562373095*(3.0*fl[7]-5.196152422706631*(fl[6]+fl[4])+9.0*fl[2]))/(20.78460969082652*EPSILON+1.414213562373095*(1.732050807568877*fl[5]-3.0*(fl[3]+fl[1])+5.196152422706631*fl[0])); 
  rCtrlL[1] = -(1.414213562373095*(3.0*fl[7]+5.196152422706631*fl[6]-1.0*(5.196152422706631*fl[4]+9.0*fl[2])))/(20.78460969082652*EPSILON+1.414213562373095*((-1.732050807568877*fl[5])+3.0*(fl[1]-1.0*fl[3])+5.196152422706631*fl[0])); 
  rCtrlL[2] = -(1.414213562373095*(3.0*fl[7]+5.196152422706631*(fl[4]-1.0*fl[6])-9.0*fl[2]))/(20.78460969082652*EPSILON+1.414213562373095*((-1.732050807568877*fl[5])+3.0*(fl[3]-1.0*fl[1])+5.196152422706631*fl[0])); 
  rCtrlL[3] = (1.414213562373095*(3.0*fl[7]+5.196152422706631*(fl[6]+fl[4])+9.0*fl[2]))/(20.78460969082652*EPSILON+1.414213562373095*(1.732050807568877*fl[5]+3.0*(fl[3]+fl[1])+5.196152422706631*fl[0])); 
  rCtrlR[0] = (1.414213562373095*(3.0*fr[7]-5.196152422706631*(fr[6]+fr[4])+9.0*fr[2]))/(20.78460969082652*EPSILON+1.414213562373095*(1.732050807568877*fr[5]-3.0*(fr[3]+fr[1])+5.196152422706631*fr[0])); 
  rCtrlR[1] = -(1.414213562373095*(3.0*fr[7]+5.196152422706631*fr[6]-1.0*(5.196152422706631*fr[4]+9.0*fr[2])))/(20.78460969082652*EPSILON+1.414213562373095*((-1.732050807568877*fr[5])+3.0*(fr[1]-1.0*fr[3])+5.196152422706631*fr[0])); 
  rCtrlR[2] = -(1.414213562373095*(3.0*fr[7]+5.196152422706631*(fr[4]-1.0*fr[6])-9.0*fr[2]))/(20.78460969082652*EPSILON+1.414213562373095*((-1.732050807568877*fr[5])+3.0*(fr[3]-1.0*fr[1])+5.196152422706631*fr[0])); 
  rCtrlR[3] = (1.414213562373095*(3.0*fr[7]+5.196152422706631*(fr[6]+fr[4])+9.0*fr[2]))/(20.78460969082652*EPSILON+1.414213562373095*(1.732050807568877*fr[5]+3.0*(fr[3]+fr[1])+5.196152422706631*fr[0])); 
  double fCtrlL[4], fCtrlR[4];  // fCtrl = anti-limited f evaluated at each control node on vx surface 
  // control node [x,vy] = [-1/3,-1/3] 
  fCtrlL[0] = 0.04811252243246882*(2.449489742783178*fl[5]-4.242640687119286*(fl[3]+fl[1])+7.348469228349534*fl[0])*limTheta(rCtrlL[0],1.0); 
  fCtrlR[0] = 0.04811252243246882*(2.449489742783178*fr[5]-4.242640687119286*(fr[3]+fr[1])+7.348469228349534*fr[0])*limTheta(rCtrlR[0],-1.0); 
  // control node [x,vy] = [1/3,-1/3] 
  fCtrlL[1] = -0.04811252243246882*(2.449489742783178*fl[5]+4.242640687119286*fl[3]-4.242640687119286*fl[1]-7.348469228349534*fl[0])*limTheta(rCtrlL[1],1.0); 
  fCtrlR[1] = -0.04811252243246882*(2.449489742783178*fr[5]+4.242640687119286*fr[3]-4.242640687119286*fr[1]-7.348469228349534*fr[0])*limTheta(rCtrlR[1],-1.0); 
  // control node [x,vy] = [-1/3,1/3] 
  fCtrlL[2] = -0.04811252243246882*(2.449489742783178*fl[5]-4.242640687119286*fl[3]+4.242640687119286*fl[1]-7.348469228349534*fl[0])*limTheta(rCtrlL[2],1.0); 
  fCtrlR[2] = -0.04811252243246882*(2.449489742783178*fr[5]-4.242640687119286*fr[3]+4.242640687119286*fr[1]-7.348469228349534*fr[0])*limTheta(rCtrlR[2],-1.0); 
  // control node [x,vy] = [1/3,1/3] 
  fCtrlL[3] = 0.04811252243246882*(2.449489742783178*fl[5]+4.242640687119286*(fl[3]+fl[1])+7.348469228349534*fl[0])*limTheta(rCtrlL[3],1.0); 
  fCtrlR[3] = 0.04811252243246882*(2.449489742783178*fr[5]+4.242640687119286*(fr[3]+fr[1])+7.348469228349534*fr[0])*limTheta(rCtrlR[3],-1.0); 
  double fL_AL[4], fR_AL[4];  // f_AL = mode coefficients of anti-limited f on surface 
  fL_AL[0] = 0.5*(fCtrlL[3]+fCtrlL[2]+fCtrlL[1]+fCtrlL[0]); 
  fL_AL[1] = 0.8660254037844386*(fCtrlL[3]-1.0*fCtrlL[2]+fCtrlL[1]-1.0*fCtrlL[0]); 
  fL_AL[2] = 0.8660254037844386*(fCtrlL[3]+fCtrlL[2]-1.0*(fCtrlL[1]+fCtrlL[0])); 
  fL_AL[3] = 1.5*(fCtrlL[3]-1.0*(fCtrlL[2]+fCtrlL[1])+fCtrlL[0]); 
  fR_AL[0] = 0.5*(fCtrlR[3]+fCtrlR[2]+fCtrlR[1]+fCtrlR[0]); 
  fR_AL[1] = 0.8660254037844386*(fCtrlR[3]-1.0*fCtrlR[2]+fCtrlR[1]-1.0*fCtrlR[0]); 
  fR_AL[2] = 0.8660254037844386*(fCtrlR[3]+fCtrlR[2]-1.0*(fCtrlR[1]+fCtrlR[0])); 
  fR_AL[3] = 1.5*(fCtrlR[3]-1.0*(fCtrlR[2]+fCtrlR[1])+fCtrlR[0]); 
  double alphaQuad; 
  // determine upwinding and enforce limiters at each surface quadrature node 
  double fhatALQuad[4], fhatAL[4]; 
  alphaQuad = 0.5*alphaDrSurf[1]-0.5*alphaDrSurf[0]; 
  if(alphaQuad > 0) {
  fhatALQuad[0] = 0.5*(fL_AL[3]-1.0*(fL_AL[2]+fL_AL[1])+fL_AL[0]); 
  } else {
  fhatALQuad[0] = 0.5*(fR_AL[3]-1.0*(fR_AL[2]+fR_AL[1])+fR_AL[0]); 
  } 
  alphaQuad = -0.5*(alphaDrSurf[1]+alphaDrSurf[0]); 
  if(alphaQuad > 0) {
  fhatALQuad[1] = -0.5*(fL_AL[3]+fL_AL[2]-1.0*(fL_AL[1]+fL_AL[0])); 
  } else {
  fhatALQuad[1] = -0.5*(fR_AL[3]+fR_AL[2]-1.0*(fR_AL[1]+fR_AL[0])); 
  } 
  alphaQuad = 0.5*alphaDrSurf[1]-0.5*alphaDrSurf[0]; 
  if(alphaQuad > 0) {
  fhatALQuad[2] = -0.5*(fL_AL[3]-1.0*fL_AL[2]+fL_AL[1]-1.0*fL_AL[0]); 
  } else {
  fhatALQuad[2] = -0.5*(fR_AL[3]-1.0*fR_AL[2]+fR_AL[1]-1.0*fR_AL[0]); 
  } 
  alphaQuad = -0.5*(alphaDrSurf[1]+alphaDrSurf[0]); 
  if(alphaQuad > 0) {
  fhatALQuad[3] = 0.5*(fL_AL[3]+fL_AL[2]+fL_AL[1]+fL_AL[0]); 
  } else {
  fhatALQuad[3] = 0.5*(fR_AL[3]+fR_AL[2]+fR_AL[1]+fR_AL[0]); 
  } 
  fhatAL[0] = 0.5*(fhatALQuad[3]+fhatALQuad[2]+fhatALQuad[1]+fhatALQuad[0]); 
  fhatAL[1] = 0.5*(fhatALQuad[3]-1.0*fhatALQuad[2]+fhatALQuad[1]-1.0*fhatALQuad[0]); 
  fhatAL[2] = 0.5*(fhatALQuad[3]+fhatALQuad[2]-1.0*(fhatALQuad[1]+fhatALQuad[0])); 
  fhatAL[3] = 0.5*(fhatALQuad[3]-1.0*(fhatALQuad[2]+fhatALQuad[1])+fhatALQuad[0]); 

  // begin surface update 
 
  double fluxFracL, fluxFracR, flim = 0.;
  double GhatDragCtrl[4];
  fluxFracL = cflRateByDirL[0] == 0. ? 0.5 : cflRateByDirL[1]/cflRateByDirL[0]; 
  fluxFracR = cflRateByDirR[0] == 0. ? 0.5 : cflRateByDirR[1]/cflRateByDirR[0]; 
  // Control node [x,mu] = [-1/3,-1/3]. 
  GhatDragCtrl[0] = (-0.1443375672974065*(alphaDrSurf[1]*fhatAL[3]+alphaDrSurf[0]*fhatAL[2]))+0.08333333333333333*(alphaDrSurf[0]*fhatAL[3]+alphaDrSurf[1]*fhatAL[2])+0.25*(alphaDrSurf[1]*fhatAL[1]+alphaDrSurf[0]*fhatAL[0])-0.1443375672974065*(alphaDrSurf[0]*fhatAL[1]+fhatAL[0]*alphaDrSurf[1]); 
  if(GhatDragCtrl[0] > EPSILON) {
    flim = std::max(0., 0.02777777777777778*(2.449489742783178*fl[7]-4.242640687119286*fl[6]+4.242640687119286*fl[5]-4.242640687119286*fl[4]-7.348469228349534*fl[3]+7.348469228349534*fl[2]-7.348469228349534*fl[1]+12.72792206135786*fl[0])); 
    GhatDragCtrl[0] = std::min(GhatDragCtrl[0], std::abs(fluxFracL*flim/dtApprox/rdv2L)); 
  } else if(GhatDragCtrl[0] < -EPSILON) {
    flim = std::max(0., -0.02777777777777778*(2.449489742783178*fr[7]-4.242640687119286*(fr[6]+fr[5]+fr[4])+7.348469228349534*(fr[3]+fr[2]+fr[1])-12.72792206135786*fr[0])); 
    GhatDragCtrl[0] = -std::min(-GhatDragCtrl[0], std::abs(fluxFracR*flim/dtApprox/rdv2R)); 
  } else GhatDragCtrl[0] = 0.; 
  // Control node [x,mu] = [1/3,-1/3]. 
  GhatDragCtrl[1] = (-0.1443375672974065*(alphaDrSurf[1]*fhatAL[3]+alphaDrSurf[0]*fhatAL[2]))-0.08333333333333333*(alphaDrSurf[0]*fhatAL[3]+alphaDrSurf[1]*fhatAL[2])+0.25*(alphaDrSurf[1]*fhatAL[1]+alphaDrSurf[0]*fhatAL[0])+0.1443375672974065*(alphaDrSurf[0]*fhatAL[1]+fhatAL[0]*alphaDrSurf[1]); 
  if(GhatDragCtrl[1] > EPSILON) {
    flim = std::max(0., -0.02777777777777778*(2.449489742783178*fl[7]+4.242640687119286*(fl[6]+fl[5])-4.242640687119286*fl[4]+7.348469228349534*fl[3]-7.348469228349534*(fl[2]+fl[1])-12.72792206135786*fl[0])); 
    GhatDragCtrl[1] = std::min(GhatDragCtrl[1], std::abs(fluxFracL*flim/dtApprox/rdv2L)); 
  } else if(GhatDragCtrl[1] < -EPSILON) {
    flim = std::max(0., 0.02777777777777778*(2.449489742783178*fr[7]+4.242640687119286*fr[6]-4.242640687119286*(fr[5]+fr[4])-7.348469228349534*(fr[3]+fr[2])+7.348469228349534*fr[1]+12.72792206135786*fr[0])); 
    GhatDragCtrl[1] = -std::min(-GhatDragCtrl[1], std::abs(fluxFracR*flim/dtApprox/rdv2R)); 
  } else GhatDragCtrl[1] = 0.; 
  // Control node [x,mu] = [-1/3,1/3]. 
  GhatDragCtrl[2] = 0.1443375672974065*(alphaDrSurf[1]*fhatAL[3]+alphaDrSurf[0]*fhatAL[2])-0.08333333333333333*(alphaDrSurf[0]*fhatAL[3]+alphaDrSurf[1]*fhatAL[2])+0.25*(alphaDrSurf[1]*fhatAL[1]+alphaDrSurf[0]*fhatAL[0])-0.1443375672974065*(alphaDrSurf[0]*fhatAL[1]+fhatAL[0]*alphaDrSurf[1]); 
  if(GhatDragCtrl[2] > EPSILON) {
    flim = std::max(0., -0.02777777777777778*(2.449489742783178*fl[7]-4.242640687119286*fl[6]+4.242640687119286*(fl[5]+fl[4])-7.348469228349534*(fl[3]+fl[2])+7.348469228349534*fl[1]-12.72792206135786*fl[0])); 
    GhatDragCtrl[2] = std::min(GhatDragCtrl[2], std::abs(fluxFracL*flim/dtApprox/rdv2L)); 
  } else if(GhatDragCtrl[2] < -EPSILON) {
    flim = std::max(0., 0.02777777777777778*(2.449489742783178*fr[7]-4.242640687119286*(fr[6]+fr[5])+4.242640687119286*fr[4]+7.348469228349534*fr[3]-7.348469228349534*(fr[2]+fr[1])+12.72792206135786*fr[0])); 
    GhatDragCtrl[2] = -std::min(-GhatDragCtrl[2], std::abs(fluxFracR*flim/dtApprox/rdv2R)); 
  } else GhatDragCtrl[2] = 0.; 
  // Control node [x,mu] = [1/3,1/3]. 
  GhatDragCtrl[3] = 0.1443375672974065*(alphaDrSurf[1]*fhatAL[3]+alphaDrSurf[0]*fhatAL[2])+0.08333333333333333*(alphaDrSurf[0]*fhatAL[3]+alphaDrSurf[1]*fhatAL[2])+0.25*(alphaDrSurf[1]*fhatAL[1]+alphaDrSurf[0]*fhatAL[0])+0.1443375672974065*(alphaDrSurf[0]*fhatAL[1]+fhatAL[0]*alphaDrSurf[1]); 
  if(GhatDragCtrl[3] > EPSILON) {
    flim = std::max(0., 0.02777777777777778*(2.449489742783178*fl[7]+4.242640687119286*(fl[6]+fl[5]+fl[4])+7.348469228349534*(fl[3]+fl[2]+fl[1])+12.72792206135786*fl[0])); 
    GhatDragCtrl[3] = std::min(GhatDragCtrl[3], std::abs(fluxFracL*flim/dtApprox/rdv2L)); 
  } else if(GhatDragCtrl[3] < -EPSILON) {
    flim = std::max(0., -0.02777777777777778*(2.449489742783178*fr[7]+4.242640687119286*fr[6]-4.242640687119286*fr[5]+4.242640687119286*fr[4]-7.348469228349534*fr[3]+7.348469228349534*fr[2]-7.348469228349534*fr[1]-12.72792206135786*fr[0])); 
    GhatDragCtrl[3] = -std::min(-GhatDragCtrl[3], std::abs(fluxFracR*flim/dtApprox/rdv2R)); 
  } else GhatDragCtrl[3] = 0.; 

  double Gdiff[8]; 
  double Ghat[8]; 
  double incr2[8]; 

  if ( ((-0.06804138174397717*fr[7])+0.06804138174397717*fl[7]+0.1178511301977579*fr[6]-0.1178511301977579*fl[6]+0.05892556509887893*fr[5]+0.05892556509887893*fl[5]+0.1178511301977579*fr[4]-0.1178511301977579*fl[4]-0.1020620726159657*fr[3]-0.1020620726159657*fl[3]-0.2041241452319315*fr[2]+0.2041241452319315*fl[2]-0.1020620726159657*fr[1]-0.1020620726159657*fl[1]+0.1767766952966368*fr[0]+0.1767766952966368*fl[0]>=0.0) && (0.06804138174397717*fr[7]-0.06804138174397717*fl[7]-0.1178511301977579*fr[6]+0.1178511301977579*fl[6]-0.05892556509887893*fr[5]-0.05892556509887893*fl[5]+0.1178511301977579*fr[4]-0.1178511301977579*fl[4]+0.1020620726159657*fr[3]+0.1020620726159657*fl[3]-0.2041241452319315*fr[2]+0.2041241452319315*fl[2]-0.1020620726159657*fr[1]-0.1020620726159657*fl[1]+0.1767766952966368*fr[0]+0.1767766952966368*fl[0]>=0.0) && (0.06804138174397717*fr[7]-0.06804138174397717*fl[7]+0.1178511301977579*fr[6]-0.1178511301977579*fl[6]-0.05892556509887893*fr[5]-0.05892556509887893*fl[5]-0.1178511301977579*fr[4]+0.1178511301977579*fl[4]-0.1020620726159657*fr[3]-0.1020620726159657*fl[3]-0.2041241452319315*fr[2]+0.2041241452319315*fl[2]+0.1020620726159657*fr[1]+0.1020620726159657*fl[1]+0.1767766952966368*fr[0]+0.1767766952966368*fl[0]>=0.0) && ((-0.06804138174397717*fr[7])+0.06804138174397717*fl[7]-0.1178511301977579*fr[6]+0.1178511301977579*fl[6]+0.05892556509887893*fr[5]+0.05892556509887893*fl[5]-0.1178511301977579*fr[4]+0.1178511301977579*fl[4]+0.1020620726159657*fr[3]+0.1020620726159657*fl[3]-0.2041241452319315*fr[2]+0.2041241452319315*fl[2]+0.1020620726159657*fr[1]+0.1020620726159657*fl[1]+0.1767766952966368*fr[0]+0.1767766952966368*fl[0]>=0.0) ) {
    incr2[2] = nuVtSqSum[1]*((-0.3535533905932737*fr[4])+0.3535533905932737*fl[4]+0.3061862178478971*(fr[1]+fl[1]))+nuVtSqSum[0]*((-0.3535533905932737*fr[2])+0.3535533905932737*fl[2]+0.3061862178478971*(fr[0]+fl[0])); 
    incr2[4] = nuVtSqSum[0]*((-0.3535533905932737*fr[4])+0.3535533905932737*fl[4]+0.3061862178478971*(fr[1]+fl[1]))+nuVtSqSum[1]*((-0.3535533905932737*fr[2])+0.3535533905932737*fl[2]+0.3061862178478971*(fr[0]+fl[0])); 
    incr2[6] = nuVtSqSum[1]*((-0.3535533905932737*fr[7])+0.3535533905932737*fl[7]+0.3061862178478971*(fr[5]+fl[5]))+nuVtSqSum[0]*((-0.3535533905932737*fr[6])+0.3535533905932737*fl[6]+0.3061862178478971*(fr[3]+fl[3])); 
    incr2[7] = nuVtSqSum[0]*((-0.3535533905932737*fr[7])+0.3535533905932737*fl[7]+0.3061862178478971*(fr[5]+fl[5]))+nuVtSqSum[1]*((-0.3535533905932737*fr[6])+0.3535533905932737*fl[6]+0.3061862178478971*(fr[3]+fl[3])); 

    Gdiff[0] = nuVtSqSum[1]*((-1.530931089239486*(fr[4]+fl[4]))+1.590990257669731*fr[1]-1.590990257669731*fl[1])+nuVtSqSum[0]*((-1.530931089239486*(fr[2]+fl[2]))+1.590990257669731*fr[0]-1.590990257669731*fl[0]); 
    Gdiff[1] = nuVtSqSum[0]*((-1.530931089239486*(fr[4]+fl[4]))+1.590990257669731*fr[1]-1.590990257669731*fl[1])+nuVtSqSum[1]*((-1.530931089239486*(fr[2]+fl[2]))+1.590990257669731*fr[0]-1.590990257669731*fl[0]); 
    Gdiff[3] = nuVtSqSum[1]*((-1.530931089239486*(fr[7]+fl[7]))+1.590990257669731*fr[5]-1.590990257669731*fl[5])+nuVtSqSum[0]*((-1.530931089239486*(fr[6]+fl[6]))+1.590990257669731*fr[3]-1.590990257669731*fl[3]); 
    Gdiff[5] = nuVtSqSum[0]*((-1.530931089239486*(fr[7]+fl[7]))+1.590990257669731*fr[5]-1.590990257669731*fl[5])+nuVtSqSum[1]*((-1.530931089239486*(fr[6]+fl[6]))+1.590990257669731*fr[3]-1.590990257669731*fl[3]); 

    Ghat[0] = Gdiff[0]*rdv+0.7071067811865475*(GhatDragCtrl[3]+GhatDragCtrl[2]+GhatDragCtrl[1]+GhatDragCtrl[0]); 
    Ghat[1] = Gdiff[1]*rdv+1.224744871391589*GhatDragCtrl[3]-1.224744871391589*GhatDragCtrl[2]+1.224744871391589*GhatDragCtrl[1]-1.224744871391589*GhatDragCtrl[0]; 
    Ghat[3] = Gdiff[3]*rdv+1.224744871391589*(GhatDragCtrl[3]+GhatDragCtrl[2])-1.224744871391589*(GhatDragCtrl[1]+GhatDragCtrl[0]); 
    Ghat[5] = Gdiff[5]*rdv+2.121320343559642*GhatDragCtrl[3]-2.121320343559642*(GhatDragCtrl[2]+GhatDragCtrl[1])+2.121320343559642*GhatDragCtrl[0]; 
  } else {

    double xBar[4];
    xBar[0] = (0.05103103630798284*fr[7]+0.05103103630798284*fl[7]-0.08838834764831849*fr[6]-0.08838834764831849*fl[6]+0.08838834764831843*fr[5]-0.08838834764831843*fl[5]-0.08838834764831849*fr[4]-0.08838834764831849*fl[4]-0.1530931089239485*fr[3]+0.1530931089239485*fl[3]+0.1530931089239485*fr[2]+0.1530931089239485*fl[2]-0.1530931089239485*fr[1]+0.1530931089239485*fl[1]+0.2651650429449551*fr[0]-0.2651650429449551*fl[0])/(0.5*(1.224744871391589*((-0.3333333333333333*(1.732050807568877*fr[4]-1.0*fr[7]))-0.5773502691896258*fr[6]+fr[2])-1.224744871391589*((-0.3333333333333333*(1.732050807568877*fl[4]-1.0*fl[7]))-0.5773502691896258*fl[6]+fl[2]))-0.25*(2.449489742783178*((-0.3333333333333333*(1.732050807568877*fr[4]-1.0*fr[7]))-0.5773502691896258*fr[6]+fr[2])-2.449489742783178*((-0.3333333333333333*(1.732050807568877*fl[4]-1.0*fl[7]))-0.5773502691896258*fl[6]+fl[2])-2.121320343559642*((-0.3333333333333333*(1.732050807568877*fr[1]-1.0*fr[5]))-0.5773502691896258*fr[3]+fr[0])-2.121320343559642*((-0.3333333333333333*(1.732050807568877*fl[1]-1.0*fl[5]))-0.5773502691896258*fl[3]+fl[0]))); 
    xBar[1] = ((-0.05103103630798284*fr[7])-0.05103103630798284*fl[7]+0.08838834764831849*fr[6]+0.08838834764831849*fl[6]-0.08838834764831843*fr[5]+0.08838834764831843*fl[5]-0.08838834764831849*fr[4]-0.08838834764831849*fl[4]+0.1530931089239485*fr[3]-0.1530931089239485*fl[3]+0.1530931089239485*fr[2]+0.1530931089239485*fl[2]-0.1530931089239485*fr[1]+0.1530931089239485*fl[1]+0.2651650429449551*fr[0]-0.2651650429449551*fl[0])/(0.5*(1.224744871391589*((-0.3333333333333333*(fr[7]+1.732050807568877*fr[4]))+0.5773502691896258*fr[6]+fr[2])-1.224744871391589*((-0.3333333333333333*(fl[7]+1.732050807568877*fl[4]))+0.5773502691896258*fl[6]+fl[2]))-0.25*(2.449489742783178*((-0.3333333333333333*(fr[7]+1.732050807568877*fr[4]))+0.5773502691896258*fr[6]+fr[2])-2.449489742783178*((-0.3333333333333333*(fl[7]+1.732050807568877*fl[4]))+0.5773502691896258*fl[6]+fl[2])-2.121320343559642*((-0.3333333333333333*(fr[5]+1.732050807568877*fr[1]))+0.5773502691896258*fr[3]+fr[0])-2.121320343559642*((-0.3333333333333333*(fl[5]+1.732050807568877*fl[1]))+0.5773502691896258*fl[3]+fl[0]))); 
    xBar[2] = ((-0.05103103630798284*fr[7])-0.05103103630798284*fl[7]-0.08838834764831849*fr[6]-0.08838834764831849*fl[6]-0.08838834764831843*fr[5]+0.08838834764831843*fl[5]+0.08838834764831849*fr[4]+0.08838834764831849*fl[4]-0.1530931089239485*fr[3]+0.1530931089239485*fl[3]+0.1530931089239485*fr[2]+0.1530931089239485*fl[2]+0.1530931089239485*fr[1]-0.1530931089239485*fl[1]+0.2651650429449551*fr[0]-0.2651650429449551*fl[0])/(0.5*(1.224744871391589*(0.3333333333333333*(1.732050807568877*fr[4]-1.0*fr[7])-0.5773502691896258*fr[6]+fr[2])-1.224744871391589*(0.3333333333333333*(1.732050807568877*fl[4]-1.0*fl[7])-0.5773502691896258*fl[6]+fl[2]))-0.25*(2.449489742783178*(0.3333333333333333*(1.732050807568877*fr[4]-1.0*fr[7])-0.5773502691896258*fr[6]+fr[2])-2.449489742783178*(0.3333333333333333*(1.732050807568877*fl[4]-1.0*fl[7])-0.5773502691896258*fl[6]+fl[2])-2.121320343559642*(0.3333333333333333*(1.732050807568877*fr[1]-1.0*fr[5])-0.5773502691896258*fr[3]+fr[0])-2.121320343559642*(0.3333333333333333*(1.732050807568877*fl[1]-1.0*fl[5])-0.5773502691896258*fl[3]+fl[0]))); 
    xBar[3] = (0.05103103630798284*fr[7]+0.05103103630798284*fl[7]+0.08838834764831849*fr[6]+0.08838834764831849*fl[6]+0.08838834764831843*fr[5]-0.08838834764831843*fl[5]+0.08838834764831849*fr[4]+0.08838834764831849*fl[4]+0.1530931089239485*fr[3]-0.1530931089239485*fl[3]+0.1530931089239485*fr[2]+0.1530931089239485*fl[2]+0.1530931089239485*fr[1]-0.1530931089239485*fl[1]+0.2651650429449551*fr[0]-0.2651650429449551*fl[0])/(0.5*(1.224744871391589*(0.3333333333333333*(fr[7]+1.732050807568877*fr[4])+0.5773502691896258*fr[6]+fr[2])-1.224744871391589*(0.3333333333333333*(fl[7]+1.732050807568877*fl[4])+0.5773502691896258*fl[6]+fl[2]))-0.25*(2.449489742783178*(0.3333333333333333*(fr[7]+1.732050807568877*fr[4])+0.5773502691896258*fr[6]+fr[2])-2.449489742783178*(0.3333333333333333*(fl[7]+1.732050807568877*fl[4])+0.5773502691896258*fl[6]+fl[2])-2.121320343559642*(0.3333333333333333*(fr[5]+1.732050807568877*fr[1])+0.5773502691896258*fr[3]+fr[0])-2.121320343559642*(0.3333333333333333*(fl[5]+1.732050807568877*fl[1])+0.5773502691896258*fl[3]+fl[0]))); 

    double xBarSq[4];
    xBarSq[0] = xBar[0]*xBar[0]; 
    xBarSq[1] = xBar[1]*xBar[1]; 
    xBarSq[2] = xBar[2]*xBar[2]; 
    xBarSq[3] = xBar[3]*xBar[3]; 

    double g1[4];
    g1[0] = (3.0*xBar[0])/(1.0-1.0*xBarSq[0])-(1.0*xBar[0]*xBarSq[0])/(1.0-1.0*xBarSq[0]); 
    g1[1] = (3.0*xBar[1])/(1.0-1.0*xBarSq[1])-(1.0*xBar[1]*xBarSq[1])/(1.0-1.0*xBarSq[1]); 
    g1[2] = (3.0*xBar[2])/(1.0-1.0*xBarSq[2])-(1.0*xBar[2]*xBarSq[2])/(1.0-1.0*xBarSq[2]); 
    g1[3] = (3.0*xBar[3])/(1.0-1.0*xBarSq[3])-(1.0*xBar[3]*xBarSq[3])/(1.0-1.0*xBarSq[3]); 

    double gBound[4];
    double gBoundP[4];

    if (std::abs(g1[0]) > 1.0e-15) {
      double g1Sq = g1[0]*g1[0];
      gBound[0] = (-(1.387778780781446e-17*g1[0]*fr[7])/std::sinh(g1[0]))+(1.387778780781446e-17*g1[0]*fl[7])/std::sinh(g1[0])+(2.775557561562891e-17*g1[0]*fr[6])/std::sinh(g1[0])-(2.775557561562891e-17*g1[0]*fl[6])/std::sinh(g1[0])+(0.05892556509887895*g1[0]*fr[5])/std::sinh(g1[0])+(0.05892556509887895*g1[0]*fl[5])/std::sinh(g1[0])+(1.387778780781446e-17*g1[0]*fr[4])/std::sinh(g1[0])-(1.387778780781446e-17*g1[0]*fl[4])/std::sinh(g1[0])-(0.1020620726159658*g1[0]*fr[3])/std::sinh(g1[0])-(0.1020620726159658*g1[0]*fl[3])/std::sinh(g1[0])-(2.775557561562891e-17*g1[0]*fr[2])/std::sinh(g1[0])+(2.775557561562891e-17*g1[0]*fl[2])/std::sinh(g1[0])-(0.1020620726159657*g1[0]*fr[1])/std::sinh(g1[0])-(0.1020620726159657*g1[0]*fl[1])/std::sinh(g1[0])+(0.1767766952966369*fr[0]*g1[0])/std::sinh(g1[0])+(0.1767766952966369*fl[0]*g1[0])/std::sinh(g1[0]); 
      gBoundP[0] = (-(1.387778780781446e-17*g1Sq*fr[7])/std::sinh(g1[0]))+(1.387778780781446e-17*g1Sq*fl[7])/std::sinh(g1[0])+(2.775557561562891e-17*g1Sq*fr[6])/std::sinh(g1[0])-(2.775557561562891e-17*g1Sq*fl[6])/std::sinh(g1[0])+(0.05892556509887895*g1Sq*fr[5])/std::sinh(g1[0])+(0.05892556509887895*g1Sq*fl[5])/std::sinh(g1[0])+(1.387778780781446e-17*g1Sq*fr[4])/std::sinh(g1[0])-(1.387778780781446e-17*g1Sq*fl[4])/std::sinh(g1[0])-(0.1020620726159658*g1Sq*fr[3])/std::sinh(g1[0])-(0.1020620726159658*g1Sq*fl[3])/std::sinh(g1[0])-(2.775557561562891e-17*g1Sq*fr[2])/std::sinh(g1[0])+(2.775557561562891e-17*g1Sq*fl[2])/std::sinh(g1[0])-(0.1020620726159657*g1Sq*fr[1])/std::sinh(g1[0])-(0.1020620726159657*g1Sq*fl[1])/std::sinh(g1[0])+(0.1767766952966369*fr[0]*g1Sq)/std::sinh(g1[0])+(0.1767766952966369*fl[0]*g1Sq)/std::sinh(g1[0]); 
    } else {
      gBound[0] = (-1.387778780781446e-17*fr[7])+1.387778780781446e-17*fl[7]+2.775557561562891e-17*fr[6]-2.775557561562891e-17*fl[6]+0.05892556509887895*fr[5]+0.05892556509887895*fl[5]+1.387778780781446e-17*fr[4]-1.387778780781446e-17*fl[4]-0.1020620726159658*fr[3]-0.1020620726159658*fl[3]-2.775557561562891e-17*fr[2]+2.775557561562891e-17*fl[2]-0.1020620726159657*fr[1]-0.1020620726159657*fl[1]+0.1767766952966369*fr[0]+0.1767766952966369*fl[0]; 
    };

    if (std::abs(g1[1]) > 1.0e-15) {
      double g1Sq = g1[1]*g1[1];
      gBound[1] = (1.387778780781446e-17*g1[1]*fr[7])/std::sinh(g1[1])-(1.387778780781446e-17*g1[1]*fl[7])/std::sinh(g1[1])-(2.775557561562891e-17*g1[1]*fr[6])/std::sinh(g1[1])+(2.775557561562891e-17*g1[1]*fl[6])/std::sinh(g1[1])-(0.05892556509887895*g1[1]*fr[5])/std::sinh(g1[1])-(0.05892556509887895*g1[1]*fl[5])/std::sinh(g1[1])+(1.387778780781446e-17*g1[1]*fr[4])/std::sinh(g1[1])-(1.387778780781446e-17*g1[1]*fl[4])/std::sinh(g1[1])+(0.1020620726159658*g1[1]*fr[3])/std::sinh(g1[1])+(0.1020620726159658*g1[1]*fl[3])/std::sinh(g1[1])-(2.775557561562891e-17*g1[1]*fr[2])/std::sinh(g1[1])+(2.775557561562891e-17*g1[1]*fl[2])/std::sinh(g1[1])-(0.1020620726159657*fr[1]*g1[1])/std::sinh(g1[1])-(0.1020620726159657*fl[1]*g1[1])/std::sinh(g1[1])+(0.1767766952966369*fr[0]*g1[1])/std::sinh(g1[1])+(0.1767766952966369*fl[0]*g1[1])/std::sinh(g1[1]); 
      gBoundP[1] = (1.387778780781446e-17*g1Sq*fr[7])/std::sinh(g1[1])-(1.387778780781446e-17*g1Sq*fl[7])/std::sinh(g1[1])-(2.775557561562891e-17*g1Sq*fr[6])/std::sinh(g1[1])+(2.775557561562891e-17*g1Sq*fl[6])/std::sinh(g1[1])-(0.05892556509887895*g1Sq*fr[5])/std::sinh(g1[1])-(0.05892556509887895*g1Sq*fl[5])/std::sinh(g1[1])+(1.387778780781446e-17*g1Sq*fr[4])/std::sinh(g1[1])-(1.387778780781446e-17*g1Sq*fl[4])/std::sinh(g1[1])+(0.1020620726159658*g1Sq*fr[3])/std::sinh(g1[1])+(0.1020620726159658*g1Sq*fl[3])/std::sinh(g1[1])-(2.775557561562891e-17*g1Sq*fr[2])/std::sinh(g1[1])+(2.775557561562891e-17*g1Sq*fl[2])/std::sinh(g1[1])-(0.1020620726159657*fr[1]*g1Sq)/std::sinh(g1[1])-(0.1020620726159657*fl[1]*g1Sq)/std::sinh(g1[1])+(0.1767766952966369*fr[0]*g1Sq)/std::sinh(g1[1])+(0.1767766952966369*fl[0]*g1Sq)/std::sinh(g1[1]); 
    } else {
      gBound[1] = 1.387778780781446e-17*fr[7]-1.387778780781446e-17*fl[7]-2.775557561562891e-17*fr[6]+2.775557561562891e-17*fl[6]-0.05892556509887895*fr[5]-0.05892556509887895*fl[5]+1.387778780781446e-17*fr[4]-1.387778780781446e-17*fl[4]+0.1020620726159658*fr[3]+0.1020620726159658*fl[3]-2.775557561562891e-17*fr[2]+2.775557561562891e-17*fl[2]-0.1020620726159657*fr[1]-0.1020620726159657*fl[1]+0.1767766952966369*fr[0]+0.1767766952966369*fl[0]; 
    };

    if (std::abs(g1[2]) > 1.0e-15) {
      double g1Sq = g1[2]*g1[2];
      gBound[2] = (1.387778780781446e-17*g1[2]*fr[7])/std::sinh(g1[2])-(1.387778780781446e-17*g1[2]*fl[7])/std::sinh(g1[2])+(2.775557561562891e-17*g1[2]*fr[6])/std::sinh(g1[2])-(2.775557561562891e-17*g1[2]*fl[6])/std::sinh(g1[2])-(0.05892556509887895*g1[2]*fr[5])/std::sinh(g1[2])-(0.05892556509887895*g1[2]*fl[5])/std::sinh(g1[2])-(1.387778780781446e-17*g1[2]*fr[4])/std::sinh(g1[2])+(1.387778780781446e-17*g1[2]*fl[4])/std::sinh(g1[2])-(0.1020620726159658*g1[2]*fr[3])/std::sinh(g1[2])-(0.1020620726159658*g1[2]*fl[3])/std::sinh(g1[2])-(2.775557561562891e-17*fr[2]*g1[2])/std::sinh(g1[2])+(2.775557561562891e-17*fl[2]*g1[2])/std::sinh(g1[2])+(0.1020620726159657*fr[1]*g1[2])/std::sinh(g1[2])+(0.1020620726159657*fl[1]*g1[2])/std::sinh(g1[2])+(0.1767766952966369*fr[0]*g1[2])/std::sinh(g1[2])+(0.1767766952966369*fl[0]*g1[2])/std::sinh(g1[2]); 
      gBoundP[2] = (1.387778780781446e-17*g1Sq*fr[7])/std::sinh(g1[2])-(1.387778780781446e-17*g1Sq*fl[7])/std::sinh(g1[2])+(2.775557561562891e-17*g1Sq*fr[6])/std::sinh(g1[2])-(2.775557561562891e-17*g1Sq*fl[6])/std::sinh(g1[2])-(0.05892556509887895*g1Sq*fr[5])/std::sinh(g1[2])-(0.05892556509887895*g1Sq*fl[5])/std::sinh(g1[2])-(1.387778780781446e-17*g1Sq*fr[4])/std::sinh(g1[2])+(1.387778780781446e-17*g1Sq*fl[4])/std::sinh(g1[2])-(0.1020620726159658*g1Sq*fr[3])/std::sinh(g1[2])-(0.1020620726159658*g1Sq*fl[3])/std::sinh(g1[2])-(2.775557561562891e-17*fr[2]*g1Sq)/std::sinh(g1[2])+(2.775557561562891e-17*fl[2]*g1Sq)/std::sinh(g1[2])+(0.1020620726159657*fr[1]*g1Sq)/std::sinh(g1[2])+(0.1020620726159657*fl[1]*g1Sq)/std::sinh(g1[2])+(0.1767766952966369*fr[0]*g1Sq)/std::sinh(g1[2])+(0.1767766952966369*fl[0]*g1Sq)/std::sinh(g1[2]); 
    } else {
      gBound[2] = 1.387778780781446e-17*fr[7]-1.387778780781446e-17*fl[7]+2.775557561562891e-17*fr[6]-2.775557561562891e-17*fl[6]-0.05892556509887895*fr[5]-0.05892556509887895*fl[5]-1.387778780781446e-17*fr[4]+1.387778780781446e-17*fl[4]-0.1020620726159658*fr[3]-0.1020620726159658*fl[3]-2.775557561562891e-17*fr[2]+2.775557561562891e-17*fl[2]+0.1020620726159657*fr[1]+0.1020620726159657*fl[1]+0.1767766952966369*fr[0]+0.1767766952966369*fl[0]; 
    };

    if (std::abs(g1[3]) > 1.0e-15) {
      double g1Sq = g1[3]*g1[3];
      gBound[3] = (-(1.387778780781446e-17*g1[3]*fr[7])/std::sinh(g1[3]))+(1.387778780781446e-17*g1[3]*fl[7])/std::sinh(g1[3])-(2.775557561562891e-17*g1[3]*fr[6])/std::sinh(g1[3])+(2.775557561562891e-17*g1[3]*fl[6])/std::sinh(g1[3])+(0.05892556509887895*g1[3]*fr[5])/std::sinh(g1[3])+(0.05892556509887895*g1[3]*fl[5])/std::sinh(g1[3])-(1.387778780781446e-17*g1[3]*fr[4])/std::sinh(g1[3])+(1.387778780781446e-17*g1[3]*fl[4])/std::sinh(g1[3])+(0.1020620726159658*fr[3]*g1[3])/std::sinh(g1[3])+(0.1020620726159658*fl[3]*g1[3])/std::sinh(g1[3])-(2.775557561562891e-17*fr[2]*g1[3])/std::sinh(g1[3])+(2.775557561562891e-17*fl[2]*g1[3])/std::sinh(g1[3])+(0.1020620726159657*fr[1]*g1[3])/std::sinh(g1[3])+(0.1020620726159657*fl[1]*g1[3])/std::sinh(g1[3])+(0.1767766952966369*fr[0]*g1[3])/std::sinh(g1[3])+(0.1767766952966369*fl[0]*g1[3])/std::sinh(g1[3]); 
      gBoundP[3] = (-(1.387778780781446e-17*g1Sq*fr[7])/std::sinh(g1[3]))+(1.387778780781446e-17*g1Sq*fl[7])/std::sinh(g1[3])-(2.775557561562891e-17*g1Sq*fr[6])/std::sinh(g1[3])+(2.775557561562891e-17*g1Sq*fl[6])/std::sinh(g1[3])+(0.05892556509887895*g1Sq*fr[5])/std::sinh(g1[3])+(0.05892556509887895*g1Sq*fl[5])/std::sinh(g1[3])-(1.387778780781446e-17*g1Sq*fr[4])/std::sinh(g1[3])+(1.387778780781446e-17*g1Sq*fl[4])/std::sinh(g1[3])+(0.1020620726159658*fr[3]*g1Sq)/std::sinh(g1[3])+(0.1020620726159658*fl[3]*g1Sq)/std::sinh(g1[3])-(2.775557561562891e-17*fr[2]*g1Sq)/std::sinh(g1[3])+(2.775557561562891e-17*fl[2]*g1Sq)/std::sinh(g1[3])+(0.1020620726159657*fr[1]*g1Sq)/std::sinh(g1[3])+(0.1020620726159657*fl[1]*g1Sq)/std::sinh(g1[3])+(0.1767766952966369*fr[0]*g1Sq)/std::sinh(g1[3])+(0.1767766952966369*fl[0]*g1Sq)/std::sinh(g1[3]); 
    } else {
      gBound[3] = (-1.387778780781446e-17*fr[7])+1.387778780781446e-17*fl[7]-2.775557561562891e-17*fr[6]+2.775557561562891e-17*fl[6]+0.05892556509887895*fr[5]+0.05892556509887895*fl[5]-1.387778780781446e-17*fr[4]+1.387778780781446e-17*fl[4]+0.1020620726159658*fr[3]+0.1020620726159658*fl[3]-2.775557561562891e-17*fr[2]+2.775557561562891e-17*fl[2]+0.1020620726159657*fr[1]+0.1020620726159657*fl[1]+0.1767766952966369*fr[0]+0.1767766952966369*fl[0]; 
    };

    incr2[2] = nuVtSqSum[1]*(0.75*(gBound[3]+gBound[2])-0.75*(gBound[1]+gBound[0]))+0.4330127018922193*nuVtSqSum[0]*(gBound[3]+gBound[2]+gBound[1]+gBound[0]); 
    incr2[4] = nuVtSqSum[0]*(0.75*(gBound[3]+gBound[2])-0.75*(gBound[1]+gBound[0]))+0.4330127018922193*nuVtSqSum[1]*(gBound[3]+gBound[2]+gBound[1]+gBound[0]); 
    incr2[6] = nuVtSqSum[1]*(1.299038105676658*gBound[3]-1.299038105676658*(gBound[2]+gBound[1])+1.299038105676658*gBound[0])+nuVtSqSum[0]*(0.75*gBound[3]-0.75*gBound[2]+0.75*gBound[1]-0.75*gBound[0]); 
    incr2[7] = nuVtSqSum[0]*(1.299038105676658*gBound[3]-1.299038105676658*(gBound[2]+gBound[1])+1.299038105676658*gBound[0])+nuVtSqSum[1]*(0.75*gBound[3]-0.75*gBound[2]+0.75*gBound[1]-0.75*gBound[0]); 

    Gdiff[0] = nuVtSqSum[1]*(0.8660254037844386*(gBoundP[3]+gBoundP[2])-0.8660254037844386*(gBoundP[1]+gBoundP[0]))+0.5*nuVtSqSum[0]*(gBoundP[3]+gBoundP[2]+gBoundP[1]+gBoundP[0]); 
    Gdiff[1] = nuVtSqSum[0]*(0.8660254037844386*(gBoundP[3]+gBoundP[2])-0.8660254037844386*(gBoundP[1]+gBoundP[0]))+0.5*nuVtSqSum[1]*(gBoundP[3]+gBoundP[2]+gBoundP[1]+gBoundP[0]); 
    Gdiff[3] = nuVtSqSum[1]*(1.5*gBoundP[3]-1.5*(gBoundP[2]+gBoundP[1])+1.5*gBoundP[0])+nuVtSqSum[0]*(0.8660254037844386*gBoundP[3]-0.8660254037844386*gBoundP[2]+0.8660254037844386*gBoundP[1]-0.8660254037844386*gBoundP[0]); 
    Gdiff[5] = nuVtSqSum[0]*(1.5*gBoundP[3]-1.5*(gBoundP[2]+gBoundP[1])+1.5*gBoundP[0])+nuVtSqSum[1]*(0.8660254037844386*gBoundP[3]-0.8660254037844386*gBoundP[2]+0.8660254037844386*gBoundP[1]-0.8660254037844386*gBoundP[0]); 

    Ghat[0] = Gdiff[0]*rdv+0.7071067811865475*(GhatDragCtrl[3]+GhatDragCtrl[2]+GhatDragCtrl[1]+GhatDragCtrl[0]); 
    Ghat[1] = Gdiff[1]*rdv+1.224744871391589*GhatDragCtrl[3]-1.224744871391589*GhatDragCtrl[2]+1.224744871391589*GhatDragCtrl[1]-1.224744871391589*GhatDragCtrl[0]; 
    Ghat[3] = Gdiff[3]*rdv+1.224744871391589*(GhatDragCtrl[3]+GhatDragCtrl[2])-1.224744871391589*(GhatDragCtrl[1]+GhatDragCtrl[0]); 
    Ghat[5] = Gdiff[5]*rdv+2.121320343559642*GhatDragCtrl[3]-2.121320343559642*(GhatDragCtrl[2]+GhatDragCtrl[1])+2.121320343559642*GhatDragCtrl[0]; 
  };

  double incr1[8]; 
  incr1[0] = -0.5*Ghat[0]; 
  incr1[1] = -0.5*Ghat[1]; 
  incr1[2] = 0.8660254037844386*Ghat[0]; 
  incr1[3] = -0.5*Ghat[3]; 
  incr1[4] = 0.8660254037844386*Ghat[1]; 
  incr1[5] = -0.5*Ghat[5]; 
  incr1[6] = 0.8660254037844386*Ghat[3]; 
  incr1[7] = 0.8660254037844386*Ghat[5]; 

  outr[0] += incr1[0]*rdv2R; 
  outr[1] += incr1[1]*rdv2R; 
  outr[2] += incr2[2]*rdvSq4R+incr1[2]*rdv2R; 
  outr[3] += incr1[3]*rdv2R; 
  outr[4] += incr2[4]*rdvSq4R+incr1[4]*rdv2R; 
  outr[5] += incr1[5]*rdv2R; 
  outr[6] += incr2[6]*rdvSq4R+incr1[6]*rdv2R; 
  outr[7] += incr2[7]*rdvSq4R+incr1[7]*rdv2R; 

  outl[0] += -1.0*incr1[0]*rdv2L; 
  outl[1] += -1.0*incr1[1]*rdv2L; 
  outl[2] += incr1[2]*rdv2L-1.0*incr2[2]*rdvSq4L; 
  outl[3] += -1.0*incr1[3]*rdv2L; 
  outl[4] += incr1[4]*rdv2L-1.0*incr2[4]*rdvSq4L; 
  outl[5] += -1.0*incr1[5]*rdv2L; 
  outl[6] += incr1[6]*rdv2L-1.0*incr2[6]*rdvSq4L; 
  outl[7] += incr1[7]*rdv2L-1.0*incr2[7]*rdvSq4L; 

  return std::abs(wl[1]-(0.7071067811865475*nuUSum[0])/nuSum); 
} 
double GkLBOconstNuSurfPositivity1x2vSer_Mu_P1(const double m_, const double *cflRateByDirL, const double *cflRateByDirR, const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double dtApprox, const double *BmagInv, const double nuSum, const double vMuMidMax, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // m_:              species mass. 
  // cflRateByDir[3]: CFL rate in each direction. 
  // w[3]:            Cell-center coordinates. 
  // dxv[3]:          Cell spacing. 
  // nuSum:           collisionalities added (self and cross species collisionalities). 
  // vMuMidMax:       maximum midpoint value of v-u. 
  // nuUSum[2]:       sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[2]:    sum of thermal speeds squared time their respective collisionalities. 
  // fl/fr:           Distribution function in left/right cells 
  // outl/outr:       Incremented distribution function in left/right cells 
  double rdv = 1.0/dxvl[2]; 
  double rdv2L = 2.0/dxvl[2]; 
  double rdv2R = 2.0/dxvr[2]; 
  double rdvSq4L = rdv2L*rdv2L; 
  double rdvSq4R = rdv2R*rdv2R; 

  double alphaDrSurf[4]; 
  alphaDrSurf[0] = (4.0*wl[2]+2.0*dxvl[2])*nuSum; 

  double rCtrlL[4], rCtrlR[4];  // rCtrl=f1/f0 at each control node in dimensions other than vy 
  rCtrlL[0] = (1.414213562373095*(3.0*fl[7]-5.196152422706631*(fl[6]+fl[5])+9.0*fl[3]))/(20.78460969082652*EPSILON+1.414213562373095*(1.732050807568877*fl[4]-3.0*(fl[2]+fl[1])+5.196152422706631*fl[0])); 
  rCtrlL[1] = -(1.414213562373095*(3.0*fl[7]+5.196152422706631*fl[6]-1.0*(5.196152422706631*fl[5]+9.0*fl[3])))/(20.78460969082652*EPSILON+1.414213562373095*((-1.732050807568877*fl[4])+3.0*(fl[1]-1.0*fl[2])+5.196152422706631*fl[0])); 
  rCtrlL[2] = -(1.414213562373095*(3.0*fl[7]+5.196152422706631*(fl[5]-1.0*fl[6])-9.0*fl[3]))/(20.78460969082652*EPSILON+1.414213562373095*((-1.732050807568877*fl[4])+3.0*(fl[2]-1.0*fl[1])+5.196152422706631*fl[0])); 
  rCtrlL[3] = (1.414213562373095*(3.0*fl[7]+5.196152422706631*(fl[6]+fl[5])+9.0*fl[3]))/(20.78460969082652*EPSILON+1.414213562373095*(1.732050807568877*fl[4]+3.0*(fl[2]+fl[1])+5.196152422706631*fl[0])); 
  rCtrlR[0] = (1.414213562373095*(3.0*fr[7]-5.196152422706631*(fr[6]+fr[5])+9.0*fr[3]))/(20.78460969082652*EPSILON+1.414213562373095*(1.732050807568877*fr[4]-3.0*(fr[2]+fr[1])+5.196152422706631*fr[0])); 
  rCtrlR[1] = -(1.414213562373095*(3.0*fr[7]+5.196152422706631*fr[6]-1.0*(5.196152422706631*fr[5]+9.0*fr[3])))/(20.78460969082652*EPSILON+1.414213562373095*((-1.732050807568877*fr[4])+3.0*(fr[1]-1.0*fr[2])+5.196152422706631*fr[0])); 
  rCtrlR[2] = -(1.414213562373095*(3.0*fr[7]+5.196152422706631*(fr[5]-1.0*fr[6])-9.0*fr[3]))/(20.78460969082652*EPSILON+1.414213562373095*((-1.732050807568877*fr[4])+3.0*(fr[2]-1.0*fr[1])+5.196152422706631*fr[0])); 
  rCtrlR[3] = (1.414213562373095*(3.0*fr[7]+5.196152422706631*(fr[6]+fr[5])+9.0*fr[3]))/(20.78460969082652*EPSILON+1.414213562373095*(1.732050807568877*fr[4]+3.0*(fr[2]+fr[1])+5.196152422706631*fr[0])); 
  double fCtrlL[4], fCtrlR[4];  // fCtrl = anti-limited f evaluated at each control node on vy surface 
  // control node [x,vx] = [-1/3,-1/3] 
  fCtrlL[0] = 0.04811252243246882*(2.449489742783178*fl[4]-4.242640687119286*(fl[2]+fl[1])+7.348469228349534*fl[0])*limTheta(rCtrlL[0],1.0); 
  fCtrlR[0] = 0.04811252243246882*(2.449489742783178*fr[4]-4.242640687119286*(fr[2]+fr[1])+7.348469228349534*fr[0])*limTheta(rCtrlR[0],-1.0); 
  // control node [x,vx] = [1/3,-1/3] 
  fCtrlL[1] = -0.04811252243246882*(2.449489742783178*fl[4]+4.242640687119286*fl[2]-4.242640687119286*fl[1]-7.348469228349534*fl[0])*limTheta(rCtrlL[1],1.0); 
  fCtrlR[1] = -0.04811252243246882*(2.449489742783178*fr[4]+4.242640687119286*fr[2]-4.242640687119286*fr[1]-7.348469228349534*fr[0])*limTheta(rCtrlR[1],-1.0); 
  // control node [x,vx] = [-1/3,1/3] 
  fCtrlL[2] = -0.04811252243246882*(2.449489742783178*fl[4]-4.242640687119286*fl[2]+4.242640687119286*fl[1]-7.348469228349534*fl[0])*limTheta(rCtrlL[2],1.0); 
  fCtrlR[2] = -0.04811252243246882*(2.449489742783178*fr[4]-4.242640687119286*fr[2]+4.242640687119286*fr[1]-7.348469228349534*fr[0])*limTheta(rCtrlR[2],-1.0); 
  // control node [x,vx] = [1/3,1/3] 
  fCtrlL[3] = 0.04811252243246882*(2.449489742783178*fl[4]+4.242640687119286*(fl[2]+fl[1])+7.348469228349534*fl[0])*limTheta(rCtrlL[3],1.0); 
  fCtrlR[3] = 0.04811252243246882*(2.449489742783178*fr[4]+4.242640687119286*(fr[2]+fr[1])+7.348469228349534*fr[0])*limTheta(rCtrlR[3],-1.0); 
  double fL_AL[4], fR_AL[4];  // f_AL = mode coefficients of anti-limited f on surface 
  fL_AL[0] = 0.5*(fCtrlL[3]+fCtrlL[2]+fCtrlL[1]+fCtrlL[0]); 
  fL_AL[1] = 0.8660254037844386*(fCtrlL[3]-1.0*fCtrlL[2]+fCtrlL[1]-1.0*fCtrlL[0]); 
  fL_AL[2] = 0.8660254037844386*(fCtrlL[3]+fCtrlL[2]-1.0*(fCtrlL[1]+fCtrlL[0])); 
  fL_AL[3] = 1.5*(fCtrlL[3]-1.0*(fCtrlL[2]+fCtrlL[1])+fCtrlL[0]); 
  fR_AL[0] = 0.5*(fCtrlR[3]+fCtrlR[2]+fCtrlR[1]+fCtrlR[0]); 
  fR_AL[1] = 0.8660254037844386*(fCtrlR[3]-1.0*fCtrlR[2]+fCtrlR[1]-1.0*fCtrlR[0]); 
  fR_AL[2] = 0.8660254037844386*(fCtrlR[3]+fCtrlR[2]-1.0*(fCtrlR[1]+fCtrlR[0])); 
  fR_AL[3] = 1.5*(fCtrlR[3]-1.0*(fCtrlR[2]+fCtrlR[1])+fCtrlR[0]); 
  double alphaQuad; 
  // determine upwinding and enforce limiters at each surface quadrature node 
  double fhatALQuad[4], fhatAL[4]; 
  alphaQuad = -0.5*alphaDrSurf[0]; 
  if(alphaQuad > 0) {
  fhatALQuad[0] = 0.5*(fL_AL[3]-1.0*(fL_AL[2]+fL_AL[1])+fL_AL[0]); 
  } else {
  fhatALQuad[0] = 0.5*(fR_AL[3]-1.0*(fR_AL[2]+fR_AL[1])+fR_AL[0]); 
  } 
  alphaQuad = -0.5*alphaDrSurf[0]; 
  if(alphaQuad > 0) {
  fhatALQuad[1] = -0.5*(fL_AL[3]+fL_AL[2]-1.0*(fL_AL[1]+fL_AL[0])); 
  } else {
  fhatALQuad[1] = -0.5*(fR_AL[3]+fR_AL[2]-1.0*(fR_AL[1]+fR_AL[0])); 
  } 
  alphaQuad = -0.5*alphaDrSurf[0]; 
  if(alphaQuad > 0) {
  fhatALQuad[2] = -0.5*(fL_AL[3]-1.0*fL_AL[2]+fL_AL[1]-1.0*fL_AL[0]); 
  } else {
  fhatALQuad[2] = -0.5*(fR_AL[3]-1.0*fR_AL[2]+fR_AL[1]-1.0*fR_AL[0]); 
  } 
  alphaQuad = -0.5*alphaDrSurf[0]; 
  if(alphaQuad > 0) {
  fhatALQuad[3] = 0.5*(fL_AL[3]+fL_AL[2]+fL_AL[1]+fL_AL[0]); 
  } else {
  fhatALQuad[3] = 0.5*(fR_AL[3]+fR_AL[2]+fR_AL[1]+fR_AL[0]); 
  } 
  fhatAL[0] = 0.5*(fhatALQuad[3]+fhatALQuad[2]+fhatALQuad[1]+fhatALQuad[0]); 
  fhatAL[1] = 0.5*(fhatALQuad[3]-1.0*fhatALQuad[2]+fhatALQuad[1]-1.0*fhatALQuad[0]); 
  fhatAL[2] = 0.5*(fhatALQuad[3]+fhatALQuad[2]-1.0*(fhatALQuad[1]+fhatALQuad[0])); 
  fhatAL[3] = 0.5*(fhatALQuad[3]-1.0*(fhatALQuad[2]+fhatALQuad[1])+fhatALQuad[0]); 

  // begin surface update 
 
  double fluxFracL, fluxFracR, flim = 0.;
  double GhatDragCtrl[4];
  fluxFracL = cflRateByDirL[0] == 0. ? 0.5 : cflRateByDirL[2]/cflRateByDirL[0]; 
  fluxFracR = cflRateByDirR[0] == 0. ? 0.5 : cflRateByDirR[2]/cflRateByDirR[0]; 
  // Control node [x,vpar] = [-1/3,-1/3]. 
  GhatDragCtrl[0] = alphaDrSurf[0]*(0.08333333333333333*fhatAL[3]-0.1443375672974065*(fhatAL[2]+fhatAL[1])+0.25*fhatAL[0]); 
  if(GhatDragCtrl[0] > EPSILON) {
    flim = std::max(0., 0.02777777777777778*(2.449489742783178*fl[7]-4.242640687119286*(fl[6]+fl[5])+4.242640687119286*fl[4]+7.348469228349534*fl[3]-7.348469228349534*(fl[2]+fl[1])+12.72792206135786*fl[0])); 
    GhatDragCtrl[0] = std::min(GhatDragCtrl[0], std::abs(fluxFracL*flim/dtApprox/rdv2L)); 
  } else if(GhatDragCtrl[0] < -EPSILON) {
    flim = std::max(0., -0.02777777777777778*(2.449489742783178*fr[7]-4.242640687119286*(fr[6]+fr[5]+fr[4])+7.348469228349534*(fr[3]+fr[2]+fr[1])-12.72792206135786*fr[0])); 
    GhatDragCtrl[0] = -std::min(-GhatDragCtrl[0], std::abs(fluxFracR*flim/dtApprox/rdv2R)); 
  } else GhatDragCtrl[0] = 0.; 
  // Control node [x,vpar] = [1/3,-1/3]. 
  GhatDragCtrl[1] = alphaDrSurf[0]*((-0.08333333333333333*fhatAL[3])-0.1443375672974065*fhatAL[2]+0.1443375672974065*fhatAL[1]+0.25*fhatAL[0]); 
  if(GhatDragCtrl[1] > EPSILON) {
    flim = std::max(0., -0.02777777777777778*(2.449489742783178*fl[7]+4.242640687119286*fl[6]-4.242640687119286*fl[5]+4.242640687119286*fl[4]-7.348469228349534*fl[3]+7.348469228349534*fl[2]-7.348469228349534*fl[1]-12.72792206135786*fl[0])); 
    GhatDragCtrl[1] = std::min(GhatDragCtrl[1], std::abs(fluxFracL*flim/dtApprox/rdv2L)); 
  } else if(GhatDragCtrl[1] < -EPSILON) {
    flim = std::max(0., 0.02777777777777778*(2.449489742783178*fr[7]+4.242640687119286*fr[6]-4.242640687119286*(fr[5]+fr[4])-7.348469228349534*(fr[3]+fr[2])+7.348469228349534*fr[1]+12.72792206135786*fr[0])); 
    GhatDragCtrl[1] = -std::min(-GhatDragCtrl[1], std::abs(fluxFracR*flim/dtApprox/rdv2R)); 
  } else GhatDragCtrl[1] = 0.; 
  // Control node [x,vpar] = [-1/3,1/3]. 
  GhatDragCtrl[2] = alphaDrSurf[0]*((-0.08333333333333333*fhatAL[3])+0.1443375672974065*fhatAL[2]-0.1443375672974065*fhatAL[1]+0.25*fhatAL[0]); 
  if(GhatDragCtrl[2] > EPSILON) {
    flim = std::max(0., -0.02777777777777778*(2.449489742783178*fl[7]-4.242640687119286*fl[6]+4.242640687119286*(fl[5]+fl[4])-7.348469228349534*(fl[3]+fl[2])+7.348469228349534*fl[1]-12.72792206135786*fl[0])); 
    GhatDragCtrl[2] = std::min(GhatDragCtrl[2], std::abs(fluxFracL*flim/dtApprox/rdv2L)); 
  } else if(GhatDragCtrl[2] < -EPSILON) {
    flim = std::max(0., 0.02777777777777778*(2.449489742783178*fr[7]-4.242640687119286*fr[6]+4.242640687119286*fr[5]-4.242640687119286*fr[4]-7.348469228349534*fr[3]+7.348469228349534*fr[2]-7.348469228349534*fr[1]+12.72792206135786*fr[0])); 
    GhatDragCtrl[2] = -std::min(-GhatDragCtrl[2], std::abs(fluxFracR*flim/dtApprox/rdv2R)); 
  } else GhatDragCtrl[2] = 0.; 
  // Control node [x,vpar] = [1/3,1/3]. 
  GhatDragCtrl[3] = alphaDrSurf[0]*(0.08333333333333333*fhatAL[3]+0.1443375672974065*(fhatAL[2]+fhatAL[1])+0.25*fhatAL[0]); 
  if(GhatDragCtrl[3] > EPSILON) {
    flim = std::max(0., 0.02777777777777778*(2.449489742783178*fl[7]+4.242640687119286*(fl[6]+fl[5]+fl[4])+7.348469228349534*(fl[3]+fl[2]+fl[1])+12.72792206135786*fl[0])); 
    GhatDragCtrl[3] = std::min(GhatDragCtrl[3], std::abs(fluxFracL*flim/dtApprox/rdv2L)); 
  } else if(GhatDragCtrl[3] < -EPSILON) {
    flim = std::max(0., -0.02777777777777778*(2.449489742783178*fr[7]+4.242640687119286*(fr[6]+fr[5])-4.242640687119286*fr[4]+7.348469228349534*fr[3]-7.348469228349534*(fr[2]+fr[1])-12.72792206135786*fr[0])); 
    GhatDragCtrl[3] = -std::min(-GhatDragCtrl[3], std::abs(fluxFracR*flim/dtApprox/rdv2R)); 
  } else GhatDragCtrl[3] = 0.; 

  double diffFac[2]; 
  diffFac[0] = (BmagInv[1]*nuVtSqSum[1]+BmagInv[0]*nuVtSqSum[0])*(1.414213562373095*wl[2]+0.7071067811865475*dxvl[2])*m_; 
  diffFac[1] = (BmagInv[0]*nuVtSqSum[1]+nuVtSqSum[0]*BmagInv[1])*(1.414213562373095*wl[2]+0.7071067811865475*dxvl[2])*m_; 

  double Gdiff[8]; 
  double Ghat[8]; 
  double incr2[8]; 

  if ( ((-0.06804138174397717*fr[7])+0.06804138174397717*fl[7]+0.1178511301977579*fr[6]-0.1178511301977579*fl[6]+0.1178511301977579*fr[5]-0.1178511301977579*fl[5]+0.05892556509887893*fr[4]+0.05892556509887893*fl[4]-0.2041241452319315*fr[3]+0.2041241452319315*fl[3]-0.1020620726159657*fr[2]-0.1020620726159657*fl[2]-0.1020620726159657*fr[1]-0.1020620726159657*fl[1]+0.1767766952966368*fr[0]+0.1767766952966368*fl[0]>=0.0) && (0.06804138174397717*fr[7]-0.06804138174397717*fl[7]-0.1178511301977579*fr[6]+0.1178511301977579*fl[6]+0.1178511301977579*fr[5]-0.1178511301977579*fl[5]-0.05892556509887893*fr[4]-0.05892556509887893*fl[4]-0.2041241452319315*fr[3]+0.2041241452319315*fl[3]+0.1020620726159657*fr[2]+0.1020620726159657*fl[2]-0.1020620726159657*fr[1]-0.1020620726159657*fl[1]+0.1767766952966368*fr[0]+0.1767766952966368*fl[0]>=0.0) && (0.06804138174397717*fr[7]-0.06804138174397717*fl[7]+0.1178511301977579*fr[6]-0.1178511301977579*fl[6]-0.1178511301977579*fr[5]+0.1178511301977579*fl[5]-0.05892556509887893*fr[4]-0.05892556509887893*fl[4]-0.2041241452319315*fr[3]+0.2041241452319315*fl[3]-0.1020620726159657*fr[2]-0.1020620726159657*fl[2]+0.1020620726159657*fr[1]+0.1020620726159657*fl[1]+0.1767766952966368*fr[0]+0.1767766952966368*fl[0]>=0.0) && ((-0.06804138174397717*fr[7])+0.06804138174397717*fl[7]-0.1178511301977579*fr[6]+0.1178511301977579*fl[6]-0.1178511301977579*fr[5]+0.1178511301977579*fl[5]+0.05892556509887893*fr[4]+0.05892556509887893*fl[4]-0.2041241452319315*fr[3]+0.2041241452319315*fl[3]+0.1020620726159657*fr[2]+0.1020620726159657*fl[2]+0.1020620726159657*fr[1]+0.1020620726159657*fl[1]+0.1767766952966368*fr[0]+0.1767766952966368*fl[0]>=0.0) ) {
    incr2[3] = diffFac[1]*(0.3535533905932737*fl[5]-0.3535533905932737*fr[5])+diffFac[0]*(0.3535533905932737*fl[3]-0.3535533905932737*fr[3])+0.3061862178478971*(diffFac[1]*(fr[1]+fl[1])+diffFac[0]*(fr[0]+fl[0])); 
    incr2[5] = diffFac[0]*(0.3535533905932737*fl[5]-0.3535533905932737*fr[5])+diffFac[1]*(0.3535533905932737*fl[3]-0.3535533905932737*fr[3])+0.3061862178478971*(diffFac[0]*(fr[1]+fl[1])+(fr[0]+fl[0])*diffFac[1]); 
    incr2[6] = diffFac[1]*(0.3535533905932737*fl[7]-0.3535533905932737*fr[7])+diffFac[0]*(0.3535533905932737*fl[6]-0.3535533905932737*fr[6])+0.3061862178478971*(diffFac[1]*(fr[4]+fl[4])+diffFac[0]*(fr[2]+fl[2])); 
    incr2[7] = diffFac[0]*(0.3535533905932737*fl[7]-0.3535533905932737*fr[7])+diffFac[1]*(0.3535533905932737*fl[6]-0.3535533905932737*fr[6])+0.3061862178478971*(diffFac[0]*(fr[4]+fl[4])+diffFac[1]*(fr[2]+fl[2])); 

    Gdiff[0] = (-1.530931089239486*(diffFac[1]*(fr[5]+fl[5])+diffFac[0]*(fr[3]+fl[3])))+diffFac[1]*(1.590990257669731*fr[1]-1.590990257669731*fl[1])+diffFac[0]*(1.590990257669731*fr[0]-1.590990257669731*fl[0]); 
    Gdiff[1] = (-1.530931089239486*(diffFac[0]*(fr[5]+fl[5])+diffFac[1]*(fr[3]+fl[3])))+diffFac[0]*(1.590990257669731*fr[1]-1.590990257669731*fl[1])+(1.590990257669731*fr[0]-1.590990257669731*fl[0])*diffFac[1]; 
    Gdiff[2] = (-1.530931089239486*(diffFac[1]*(fr[7]+fl[7])+diffFac[0]*(fr[6]+fl[6])))+diffFac[1]*(1.590990257669731*fr[4]-1.590990257669731*fl[4])+diffFac[0]*(1.590990257669731*fr[2]-1.590990257669731*fl[2]); 
    Gdiff[4] = (-1.530931089239486*(diffFac[0]*(fr[7]+fl[7])+diffFac[1]*(fr[6]+fl[6])))+diffFac[0]*(1.590990257669731*fr[4]-1.590990257669731*fl[4])+diffFac[1]*(1.590990257669731*fr[2]-1.590990257669731*fl[2]); 

    Ghat[0] = Gdiff[0]*rdv+0.7071067811865475*(GhatDragCtrl[3]+GhatDragCtrl[2]+GhatDragCtrl[1]+GhatDragCtrl[0]); 
    Ghat[1] = Gdiff[1]*rdv+1.224744871391589*GhatDragCtrl[3]-1.224744871391589*GhatDragCtrl[2]+1.224744871391589*GhatDragCtrl[1]-1.224744871391589*GhatDragCtrl[0]; 
    Ghat[2] = Gdiff[2]*rdv+1.224744871391589*(GhatDragCtrl[3]+GhatDragCtrl[2])-1.224744871391589*(GhatDragCtrl[1]+GhatDragCtrl[0]); 
    Ghat[4] = Gdiff[4]*rdv+2.121320343559642*GhatDragCtrl[3]-2.121320343559642*(GhatDragCtrl[2]+GhatDragCtrl[1])+2.121320343559642*GhatDragCtrl[0]; 
  } else {

    double xBar[4];
    xBar[0] = (0.05103103630798284*fr[7]+0.05103103630798284*fl[7]-0.08838834764831849*fr[6]-0.08838834764831849*fl[6]-0.08838834764831849*fr[5]-0.08838834764831849*fl[5]+0.08838834764831843*fr[4]-0.08838834764831843*fl[4]+0.1530931089239485*fr[3]+0.1530931089239485*fl[3]-0.1530931089239485*fr[2]+0.1530931089239485*fl[2]-0.1530931089239485*fr[1]+0.1530931089239485*fl[1]+0.2651650429449551*fr[0]-0.2651650429449551*fl[0])/(0.5*(1.224744871391589*((-0.3333333333333333*(1.732050807568877*fr[5]-1.0*fr[7]))-0.5773502691896258*fr[6]+fr[3])-1.224744871391589*((-0.3333333333333333*(1.732050807568877*fl[5]-1.0*fl[7]))-0.5773502691896258*fl[6]+fl[3]))-0.25*(2.449489742783178*((-0.3333333333333333*(1.732050807568877*fr[5]-1.0*fr[7]))-0.5773502691896258*fr[6]+fr[3])-2.449489742783178*((-0.3333333333333333*(1.732050807568877*fl[5]-1.0*fl[7]))-0.5773502691896258*fl[6]+fl[3])-2.121320343559642*((-0.3333333333333333*(1.732050807568877*fr[1]-1.0*fr[4]))-0.5773502691896258*fr[2]+fr[0])-2.121320343559642*((-0.3333333333333333*(1.732050807568877*fl[1]-1.0*fl[4]))-0.5773502691896258*fl[2]+fl[0]))); 
    xBar[1] = ((-0.05103103630798284*fr[7])-0.05103103630798284*fl[7]+0.08838834764831849*fr[6]+0.08838834764831849*fl[6]-0.08838834764831849*fr[5]-0.08838834764831849*fl[5]-0.08838834764831843*fr[4]+0.08838834764831843*fl[4]+0.1530931089239485*fr[3]+0.1530931089239485*fl[3]+0.1530931089239485*fr[2]-0.1530931089239485*fl[2]-0.1530931089239485*fr[1]+0.1530931089239485*fl[1]+0.2651650429449551*fr[0]-0.2651650429449551*fl[0])/(0.5*(1.224744871391589*((-0.3333333333333333*(fr[7]+1.732050807568877*fr[5]))+0.5773502691896258*fr[6]+fr[3])-1.224744871391589*((-0.3333333333333333*(fl[7]+1.732050807568877*fl[5]))+0.5773502691896258*fl[6]+fl[3]))-0.25*(2.449489742783178*((-0.3333333333333333*(fr[7]+1.732050807568877*fr[5]))+0.5773502691896258*fr[6]+fr[3])-2.449489742783178*((-0.3333333333333333*(fl[7]+1.732050807568877*fl[5]))+0.5773502691896258*fl[6]+fl[3])-2.121320343559642*((-0.3333333333333333*(fr[4]+1.732050807568877*fr[1]))+0.5773502691896258*fr[2]+fr[0])-2.121320343559642*((-0.3333333333333333*(fl[4]+1.732050807568877*fl[1]))+0.5773502691896258*fl[2]+fl[0]))); 
    xBar[2] = ((-0.05103103630798284*fr[7])-0.05103103630798284*fl[7]-0.08838834764831849*fr[6]-0.08838834764831849*fl[6]+0.08838834764831849*fr[5]+0.08838834764831849*fl[5]-0.08838834764831843*fr[4]+0.08838834764831843*fl[4]+0.1530931089239485*fr[3]+0.1530931089239485*fl[3]-0.1530931089239485*fr[2]+0.1530931089239485*fl[2]+0.1530931089239485*fr[1]-0.1530931089239485*fl[1]+0.2651650429449551*fr[0]-0.2651650429449551*fl[0])/(0.5*(1.224744871391589*(0.3333333333333333*(1.732050807568877*fr[5]-1.0*fr[7])-0.5773502691896258*fr[6]+fr[3])-1.224744871391589*(0.3333333333333333*(1.732050807568877*fl[5]-1.0*fl[7])-0.5773502691896258*fl[6]+fl[3]))-0.25*(2.449489742783178*(0.3333333333333333*(1.732050807568877*fr[5]-1.0*fr[7])-0.5773502691896258*fr[6]+fr[3])-2.449489742783178*(0.3333333333333333*(1.732050807568877*fl[5]-1.0*fl[7])-0.5773502691896258*fl[6]+fl[3])-2.121320343559642*(0.3333333333333333*(1.732050807568877*fr[1]-1.0*fr[4])-0.5773502691896258*fr[2]+fr[0])-2.121320343559642*(0.3333333333333333*(1.732050807568877*fl[1]-1.0*fl[4])-0.5773502691896258*fl[2]+fl[0]))); 
    xBar[3] = (0.05103103630798284*fr[7]+0.05103103630798284*fl[7]+0.08838834764831849*fr[6]+0.08838834764831849*fl[6]+0.08838834764831849*fr[5]+0.08838834764831849*fl[5]+0.08838834764831843*fr[4]-0.08838834764831843*fl[4]+0.1530931089239485*fr[3]+0.1530931089239485*fl[3]+0.1530931089239485*fr[2]-0.1530931089239485*fl[2]+0.1530931089239485*fr[1]-0.1530931089239485*fl[1]+0.2651650429449551*fr[0]-0.2651650429449551*fl[0])/(0.5*(1.224744871391589*(0.3333333333333333*(fr[7]+1.732050807568877*fr[5])+0.5773502691896258*fr[6]+fr[3])-1.224744871391589*(0.3333333333333333*(fl[7]+1.732050807568877*fl[5])+0.5773502691896258*fl[6]+fl[3]))-0.25*(2.449489742783178*(0.3333333333333333*(fr[7]+1.732050807568877*fr[5])+0.5773502691896258*fr[6]+fr[3])-2.449489742783178*(0.3333333333333333*(fl[7]+1.732050807568877*fl[5])+0.5773502691896258*fl[6]+fl[3])-2.121320343559642*(0.3333333333333333*(fr[4]+1.732050807568877*fr[1])+0.5773502691896258*fr[2]+fr[0])-2.121320343559642*(0.3333333333333333*(fl[4]+1.732050807568877*fl[1])+0.5773502691896258*fl[2]+fl[0]))); 

    double xBarSq[4];
    xBarSq[0] = xBar[0]*xBar[0]; 
    xBarSq[1] = xBar[1]*xBar[1]; 
    xBarSq[2] = xBar[2]*xBar[2]; 
    xBarSq[3] = xBar[3]*xBar[3]; 

    double g1[4];
    g1[0] = (3.0*xBar[0])/(1.0-1.0*xBarSq[0])-(1.0*xBar[0]*xBarSq[0])/(1.0-1.0*xBarSq[0]); 
    g1[1] = (3.0*xBar[1])/(1.0-1.0*xBarSq[1])-(1.0*xBar[1]*xBarSq[1])/(1.0-1.0*xBarSq[1]); 
    g1[2] = (3.0*xBar[2])/(1.0-1.0*xBarSq[2])-(1.0*xBar[2]*xBarSq[2])/(1.0-1.0*xBarSq[2]); 
    g1[3] = (3.0*xBar[3])/(1.0-1.0*xBarSq[3])-(1.0*xBar[3]*xBarSq[3])/(1.0-1.0*xBarSq[3]); 

    double gBound[4];
    double gBoundP[4];

    if (std::abs(g1[0]) > 1.0e-15) {
      double g1Sq = g1[0]*g1[0];
      gBound[0] = (-(1.387778780781446e-17*g1[0]*fr[7])/std::sinh(g1[0]))+(1.387778780781446e-17*g1[0]*fl[7])/std::sinh(g1[0])+(2.775557561562891e-17*g1[0]*fr[6])/std::sinh(g1[0])-(2.775557561562891e-17*g1[0]*fl[6])/std::sinh(g1[0])+(1.387778780781446e-17*g1[0]*fr[5])/std::sinh(g1[0])-(1.387778780781446e-17*g1[0]*fl[5])/std::sinh(g1[0])+(0.05892556509887895*g1[0]*fr[4])/std::sinh(g1[0])+(0.05892556509887895*g1[0]*fl[4])/std::sinh(g1[0])-(2.775557561562891e-17*g1[0]*fr[3])/std::sinh(g1[0])+(2.775557561562891e-17*g1[0]*fl[3])/std::sinh(g1[0])-(0.1020620726159658*g1[0]*fr[2])/std::sinh(g1[0])-(0.1020620726159658*g1[0]*fl[2])/std::sinh(g1[0])-(0.1020620726159657*g1[0]*fr[1])/std::sinh(g1[0])-(0.1020620726159657*g1[0]*fl[1])/std::sinh(g1[0])+(0.1767766952966369*fr[0]*g1[0])/std::sinh(g1[0])+(0.1767766952966369*fl[0]*g1[0])/std::sinh(g1[0]); 
      gBoundP[0] = (-(1.387778780781446e-17*g1Sq*fr[7])/std::sinh(g1[0]))+(1.387778780781446e-17*g1Sq*fl[7])/std::sinh(g1[0])+(2.775557561562891e-17*g1Sq*fr[6])/std::sinh(g1[0])-(2.775557561562891e-17*g1Sq*fl[6])/std::sinh(g1[0])+(1.387778780781446e-17*g1Sq*fr[5])/std::sinh(g1[0])-(1.387778780781446e-17*g1Sq*fl[5])/std::sinh(g1[0])+(0.05892556509887895*g1Sq*fr[4])/std::sinh(g1[0])+(0.05892556509887895*g1Sq*fl[4])/std::sinh(g1[0])-(2.775557561562891e-17*g1Sq*fr[3])/std::sinh(g1[0])+(2.775557561562891e-17*g1Sq*fl[3])/std::sinh(g1[0])-(0.1020620726159658*g1Sq*fr[2])/std::sinh(g1[0])-(0.1020620726159658*g1Sq*fl[2])/std::sinh(g1[0])-(0.1020620726159657*g1Sq*fr[1])/std::sinh(g1[0])-(0.1020620726159657*g1Sq*fl[1])/std::sinh(g1[0])+(0.1767766952966369*fr[0]*g1Sq)/std::sinh(g1[0])+(0.1767766952966369*fl[0]*g1Sq)/std::sinh(g1[0]); 
    } else {
      gBound[0] = (-1.387778780781446e-17*fr[7])+1.387778780781446e-17*fl[7]+2.775557561562891e-17*fr[6]-2.775557561562891e-17*fl[6]+1.387778780781446e-17*fr[5]-1.387778780781446e-17*fl[5]+0.05892556509887895*fr[4]+0.05892556509887895*fl[4]-2.775557561562891e-17*fr[3]+2.775557561562891e-17*fl[3]-0.1020620726159658*fr[2]-0.1020620726159658*fl[2]-0.1020620726159657*fr[1]-0.1020620726159657*fl[1]+0.1767766952966369*fr[0]+0.1767766952966369*fl[0]; 
    };

    if (std::abs(g1[1]) > 1.0e-15) {
      double g1Sq = g1[1]*g1[1];
      gBound[1] = (1.387778780781446e-17*g1[1]*fr[7])/std::sinh(g1[1])-(1.387778780781446e-17*g1[1]*fl[7])/std::sinh(g1[1])-(2.775557561562891e-17*g1[1]*fr[6])/std::sinh(g1[1])+(2.775557561562891e-17*g1[1]*fl[6])/std::sinh(g1[1])+(1.387778780781446e-17*g1[1]*fr[5])/std::sinh(g1[1])-(1.387778780781446e-17*g1[1]*fl[5])/std::sinh(g1[1])-(0.05892556509887895*g1[1]*fr[4])/std::sinh(g1[1])-(0.05892556509887895*g1[1]*fl[4])/std::sinh(g1[1])-(2.775557561562891e-17*g1[1]*fr[3])/std::sinh(g1[1])+(2.775557561562891e-17*g1[1]*fl[3])/std::sinh(g1[1])+(0.1020620726159658*g1[1]*fr[2])/std::sinh(g1[1])+(0.1020620726159658*g1[1]*fl[2])/std::sinh(g1[1])-(0.1020620726159657*fr[1]*g1[1])/std::sinh(g1[1])-(0.1020620726159657*fl[1]*g1[1])/std::sinh(g1[1])+(0.1767766952966369*fr[0]*g1[1])/std::sinh(g1[1])+(0.1767766952966369*fl[0]*g1[1])/std::sinh(g1[1]); 
      gBoundP[1] = (1.387778780781446e-17*g1Sq*fr[7])/std::sinh(g1[1])-(1.387778780781446e-17*g1Sq*fl[7])/std::sinh(g1[1])-(2.775557561562891e-17*g1Sq*fr[6])/std::sinh(g1[1])+(2.775557561562891e-17*g1Sq*fl[6])/std::sinh(g1[1])+(1.387778780781446e-17*g1Sq*fr[5])/std::sinh(g1[1])-(1.387778780781446e-17*g1Sq*fl[5])/std::sinh(g1[1])-(0.05892556509887895*g1Sq*fr[4])/std::sinh(g1[1])-(0.05892556509887895*g1Sq*fl[4])/std::sinh(g1[1])-(2.775557561562891e-17*g1Sq*fr[3])/std::sinh(g1[1])+(2.775557561562891e-17*g1Sq*fl[3])/std::sinh(g1[1])+(0.1020620726159658*g1Sq*fr[2])/std::sinh(g1[1])+(0.1020620726159658*g1Sq*fl[2])/std::sinh(g1[1])-(0.1020620726159657*fr[1]*g1Sq)/std::sinh(g1[1])-(0.1020620726159657*fl[1]*g1Sq)/std::sinh(g1[1])+(0.1767766952966369*fr[0]*g1Sq)/std::sinh(g1[1])+(0.1767766952966369*fl[0]*g1Sq)/std::sinh(g1[1]); 
    } else {
      gBound[1] = 1.387778780781446e-17*fr[7]-1.387778780781446e-17*fl[7]-2.775557561562891e-17*fr[6]+2.775557561562891e-17*fl[6]+1.387778780781446e-17*fr[5]-1.387778780781446e-17*fl[5]-0.05892556509887895*fr[4]-0.05892556509887895*fl[4]-2.775557561562891e-17*fr[3]+2.775557561562891e-17*fl[3]+0.1020620726159658*fr[2]+0.1020620726159658*fl[2]-0.1020620726159657*fr[1]-0.1020620726159657*fl[1]+0.1767766952966369*fr[0]+0.1767766952966369*fl[0]; 
    };

    if (std::abs(g1[2]) > 1.0e-15) {
      double g1Sq = g1[2]*g1[2];
      gBound[2] = (1.387778780781446e-17*g1[2]*fr[7])/std::sinh(g1[2])-(1.387778780781446e-17*g1[2]*fl[7])/std::sinh(g1[2])+(2.775557561562891e-17*g1[2]*fr[6])/std::sinh(g1[2])-(2.775557561562891e-17*g1[2]*fl[6])/std::sinh(g1[2])-(1.387778780781446e-17*g1[2]*fr[5])/std::sinh(g1[2])+(1.387778780781446e-17*g1[2]*fl[5])/std::sinh(g1[2])-(0.05892556509887895*g1[2]*fr[4])/std::sinh(g1[2])-(0.05892556509887895*g1[2]*fl[4])/std::sinh(g1[2])-(2.775557561562891e-17*g1[2]*fr[3])/std::sinh(g1[2])+(2.775557561562891e-17*g1[2]*fl[3])/std::sinh(g1[2])-(0.1020620726159658*fr[2]*g1[2])/std::sinh(g1[2])-(0.1020620726159658*fl[2]*g1[2])/std::sinh(g1[2])+(0.1020620726159657*fr[1]*g1[2])/std::sinh(g1[2])+(0.1020620726159657*fl[1]*g1[2])/std::sinh(g1[2])+(0.1767766952966369*fr[0]*g1[2])/std::sinh(g1[2])+(0.1767766952966369*fl[0]*g1[2])/std::sinh(g1[2]); 
      gBoundP[2] = (1.387778780781446e-17*g1Sq*fr[7])/std::sinh(g1[2])-(1.387778780781446e-17*g1Sq*fl[7])/std::sinh(g1[2])+(2.775557561562891e-17*g1Sq*fr[6])/std::sinh(g1[2])-(2.775557561562891e-17*g1Sq*fl[6])/std::sinh(g1[2])-(1.387778780781446e-17*g1Sq*fr[5])/std::sinh(g1[2])+(1.387778780781446e-17*g1Sq*fl[5])/std::sinh(g1[2])-(0.05892556509887895*g1Sq*fr[4])/std::sinh(g1[2])-(0.05892556509887895*g1Sq*fl[4])/std::sinh(g1[2])-(2.775557561562891e-17*g1Sq*fr[3])/std::sinh(g1[2])+(2.775557561562891e-17*g1Sq*fl[3])/std::sinh(g1[2])-(0.1020620726159658*fr[2]*g1Sq)/std::sinh(g1[2])-(0.1020620726159658*fl[2]*g1Sq)/std::sinh(g1[2])+(0.1020620726159657*fr[1]*g1Sq)/std::sinh(g1[2])+(0.1020620726159657*fl[1]*g1Sq)/std::sinh(g1[2])+(0.1767766952966369*fr[0]*g1Sq)/std::sinh(g1[2])+(0.1767766952966369*fl[0]*g1Sq)/std::sinh(g1[2]); 
    } else {
      gBound[2] = 1.387778780781446e-17*fr[7]-1.387778780781446e-17*fl[7]+2.775557561562891e-17*fr[6]-2.775557561562891e-17*fl[6]-1.387778780781446e-17*fr[5]+1.387778780781446e-17*fl[5]-0.05892556509887895*fr[4]-0.05892556509887895*fl[4]-2.775557561562891e-17*fr[3]+2.775557561562891e-17*fl[3]-0.1020620726159658*fr[2]-0.1020620726159658*fl[2]+0.1020620726159657*fr[1]+0.1020620726159657*fl[1]+0.1767766952966369*fr[0]+0.1767766952966369*fl[0]; 
    };

    if (std::abs(g1[3]) > 1.0e-15) {
      double g1Sq = g1[3]*g1[3];
      gBound[3] = (-(1.387778780781446e-17*g1[3]*fr[7])/std::sinh(g1[3]))+(1.387778780781446e-17*g1[3]*fl[7])/std::sinh(g1[3])-(2.775557561562891e-17*g1[3]*fr[6])/std::sinh(g1[3])+(2.775557561562891e-17*g1[3]*fl[6])/std::sinh(g1[3])-(1.387778780781446e-17*g1[3]*fr[5])/std::sinh(g1[3])+(1.387778780781446e-17*g1[3]*fl[5])/std::sinh(g1[3])+(0.05892556509887895*g1[3]*fr[4])/std::sinh(g1[3])+(0.05892556509887895*g1[3]*fl[4])/std::sinh(g1[3])-(2.775557561562891e-17*fr[3]*g1[3])/std::sinh(g1[3])+(2.775557561562891e-17*fl[3]*g1[3])/std::sinh(g1[3])+(0.1020620726159658*fr[2]*g1[3])/std::sinh(g1[3])+(0.1020620726159658*fl[2]*g1[3])/std::sinh(g1[3])+(0.1020620726159657*fr[1]*g1[3])/std::sinh(g1[3])+(0.1020620726159657*fl[1]*g1[3])/std::sinh(g1[3])+(0.1767766952966369*fr[0]*g1[3])/std::sinh(g1[3])+(0.1767766952966369*fl[0]*g1[3])/std::sinh(g1[3]); 
      gBoundP[3] = (-(1.387778780781446e-17*g1Sq*fr[7])/std::sinh(g1[3]))+(1.387778780781446e-17*g1Sq*fl[7])/std::sinh(g1[3])-(2.775557561562891e-17*g1Sq*fr[6])/std::sinh(g1[3])+(2.775557561562891e-17*g1Sq*fl[6])/std::sinh(g1[3])-(1.387778780781446e-17*g1Sq*fr[5])/std::sinh(g1[3])+(1.387778780781446e-17*g1Sq*fl[5])/std::sinh(g1[3])+(0.05892556509887895*g1Sq*fr[4])/std::sinh(g1[3])+(0.05892556509887895*g1Sq*fl[4])/std::sinh(g1[3])-(2.775557561562891e-17*fr[3]*g1Sq)/std::sinh(g1[3])+(2.775557561562891e-17*fl[3]*g1Sq)/std::sinh(g1[3])+(0.1020620726159658*fr[2]*g1Sq)/std::sinh(g1[3])+(0.1020620726159658*fl[2]*g1Sq)/std::sinh(g1[3])+(0.1020620726159657*fr[1]*g1Sq)/std::sinh(g1[3])+(0.1020620726159657*fl[1]*g1Sq)/std::sinh(g1[3])+(0.1767766952966369*fr[0]*g1Sq)/std::sinh(g1[3])+(0.1767766952966369*fl[0]*g1Sq)/std::sinh(g1[3]); 
    } else {
      gBound[3] = (-1.387778780781446e-17*fr[7])+1.387778780781446e-17*fl[7]-2.775557561562891e-17*fr[6]+2.775557561562891e-17*fl[6]-1.387778780781446e-17*fr[5]+1.387778780781446e-17*fl[5]+0.05892556509887895*fr[4]+0.05892556509887895*fl[4]-2.775557561562891e-17*fr[3]+2.775557561562891e-17*fl[3]+0.1020620726159658*fr[2]+0.1020620726159658*fl[2]+0.1020620726159657*fr[1]+0.1020620726159657*fl[1]+0.1767766952966369*fr[0]+0.1767766952966369*fl[0]; 
    };

    incr2[3] = (0.75*diffFac[1]+0.4330127018922193*diffFac[0])*gBound[3]+(0.75*diffFac[1]+0.4330127018922193*diffFac[0])*gBound[2]+(0.4330127018922193*diffFac[0]-0.75*diffFac[1])*gBound[1]+gBound[0]*(0.4330127018922193*diffFac[0]-0.75*diffFac[1]); 
    incr2[5] = (0.4330127018922193*diffFac[1]+0.75*diffFac[0])*gBound[3]+(0.4330127018922193*diffFac[1]+0.75*diffFac[0])*gBound[2]+(0.4330127018922193*diffFac[1]-0.75*diffFac[0])*gBound[1]+gBound[0]*(0.4330127018922193*diffFac[1]-0.75*diffFac[0]); 
    incr2[6] = (1.299038105676658*diffFac[1]+0.75*diffFac[0])*gBound[3]+((-1.299038105676658*diffFac[1])-0.75*diffFac[0])*gBound[2]+(0.75*diffFac[0]-1.299038105676658*diffFac[1])*gBound[1]+gBound[0]*(1.299038105676658*diffFac[1]-0.75*diffFac[0]); 
    incr2[7] = (0.75*diffFac[1]+1.299038105676658*diffFac[0])*gBound[3]+((-0.75*diffFac[1])-1.299038105676658*diffFac[0])*gBound[2]+(0.75*diffFac[1]-1.299038105676658*diffFac[0])*gBound[1]+gBound[0]*(1.299038105676658*diffFac[0]-0.75*diffFac[1]); 

    Gdiff[0] = (0.8660254037844386*diffFac[1]+0.5*diffFac[0])*gBoundP[3]+(0.8660254037844386*diffFac[1]+0.5*diffFac[0])*gBoundP[2]+(0.5*diffFac[0]-0.8660254037844386*diffFac[1])*gBoundP[1]+gBoundP[0]*(0.5*diffFac[0]-0.8660254037844386*diffFac[1]); 
    Gdiff[1] = (0.5*diffFac[1]+0.8660254037844386*diffFac[0])*gBoundP[3]+(0.5*diffFac[1]+0.8660254037844386*diffFac[0])*gBoundP[2]+(0.5*diffFac[1]-0.8660254037844386*diffFac[0])*gBoundP[1]+gBoundP[0]*(0.5*diffFac[1]-0.8660254037844386*diffFac[0]); 
    Gdiff[2] = (1.5*diffFac[1]+0.8660254037844386*diffFac[0])*gBoundP[3]+((-1.5*diffFac[1])-0.8660254037844386*diffFac[0])*gBoundP[2]+(0.8660254037844386*diffFac[0]-1.5*diffFac[1])*gBoundP[1]+gBoundP[0]*(1.5*diffFac[1]-0.8660254037844386*diffFac[0]); 
    Gdiff[4] = (0.8660254037844386*diffFac[1]+1.5*diffFac[0])*gBoundP[3]+((-0.8660254037844386*diffFac[1])-1.5*diffFac[0])*gBoundP[2]+(0.8660254037844386*diffFac[1]-1.5*diffFac[0])*gBoundP[1]+gBoundP[0]*(1.5*diffFac[0]-0.8660254037844386*diffFac[1]); 

    Ghat[0] = Gdiff[0]*rdv+0.7071067811865475*(GhatDragCtrl[3]+GhatDragCtrl[2]+GhatDragCtrl[1]+GhatDragCtrl[0]); 
    Ghat[1] = Gdiff[1]*rdv+1.224744871391589*GhatDragCtrl[3]-1.224744871391589*GhatDragCtrl[2]+1.224744871391589*GhatDragCtrl[1]-1.224744871391589*GhatDragCtrl[0]; 
    Ghat[2] = Gdiff[2]*rdv+1.224744871391589*(GhatDragCtrl[3]+GhatDragCtrl[2])-1.224744871391589*(GhatDragCtrl[1]+GhatDragCtrl[0]); 
    Ghat[4] = Gdiff[4]*rdv+2.121320343559642*GhatDragCtrl[3]-2.121320343559642*(GhatDragCtrl[2]+GhatDragCtrl[1])+2.121320343559642*GhatDragCtrl[0]; 
  };

  double incr1[8]; 
  incr1[0] = -0.5*Ghat[0]; 
  incr1[1] = -0.5*Ghat[1]; 
  incr1[2] = -0.5*Ghat[2]; 
  incr1[3] = 0.8660254037844386*Ghat[0]; 
  incr1[4] = -0.5*Ghat[4]; 
  incr1[5] = 0.8660254037844386*Ghat[1]; 
  incr1[6] = 0.8660254037844386*Ghat[2]; 
  incr1[7] = 0.8660254037844386*Ghat[4]; 

  outr[0] += incr1[0]*rdv2R; 
  outr[1] += incr1[1]*rdv2R; 
  outr[2] += incr1[2]*rdv2R; 
  outr[3] += incr2[3]*rdvSq4R+incr1[3]*rdv2R; 
  outr[4] += incr1[4]*rdv2R; 
  outr[5] += incr2[5]*rdvSq4R+incr1[5]*rdv2R; 
  outr[6] += incr2[6]*rdvSq4R+incr1[6]*rdv2R; 
  outr[7] += incr2[7]*rdvSq4R+incr1[7]*rdv2R; 

  outl[0] += -1.0*incr1[0]*rdv2L; 
  outl[1] += -1.0*incr1[1]*rdv2L; 
  outl[2] += -1.0*incr1[2]*rdv2L; 
  outl[3] += incr1[3]*rdv2L-1.0*incr2[3]*rdvSq4L; 
  outl[4] += -1.0*incr1[4]*rdv2L; 
  outl[5] += incr1[5]*rdv2L-1.0*incr2[5]*rdvSq4L; 
  outl[6] += incr1[6]*rdv2L-1.0*incr2[6]*rdvSq4L; 
  outl[7] += incr1[7]*rdv2L-1.0*incr2[7]*rdvSq4L; 

  return std::abs(2.0*wl[2]); 
} 
