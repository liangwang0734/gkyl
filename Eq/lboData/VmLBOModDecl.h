#ifndef VMLBO_MOD_DELC_H 
#define VMLBO_MOD_DELC_H 
#include <cmath> 
extern "C" { 
double VmLBOconstNuVol1x1vMaxP1(const double *w, const double *dxv, const double nu, const double *u, const double *vtSq, const double *f, double *out); 
double VmLBOconstNuSurf1x1vMax_VX_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr); 

double VmLBOconstNuVol1x1vMaxP2(const double *w, const double *dxv, const double nu, const double *u, const double *vtSq, const double *f, double *out); 
double VmLBOconstNuSurf1x1vMax_VX_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr); 

double VmLBOconstNuVol1x1vMaxP3(const double *w, const double *dxv, const double nu, const double *u, const double *vtSq, const double *f, double *out); 
double VmLBOconstNuSurf1x1vMax_VX_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr); 

double VmLBOconstNuVol1x1vMaxP4(const double *w, const double *dxv, const double nu, const double *u, const double *vtSq, const double *f, double *out); 
double VmLBOconstNuSurf1x1vMax_VX_P4(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr); 


double VmLBOconstNuVol1x2vMaxP1(const double *w, const double *dxv, const double nu, const double *u, const double *vtSq, const double *f, double *out); 
double VmLBOconstNuSurf1x2vMax_VX_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr); 
double VmLBOconstNuSurf1x2vMax_VY_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr); 

double VmLBOconstNuVol1x2vMaxP2(const double *w, const double *dxv, const double nu, const double *u, const double *vtSq, const double *f, double *out); 
double VmLBOconstNuSurf1x2vMax_VX_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr); 
double VmLBOconstNuSurf1x2vMax_VY_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr); 

double VmLBOconstNuVol1x2vMaxP3(const double *w, const double *dxv, const double nu, const double *u, const double *vtSq, const double *f, double *out); 
double VmLBOconstNuSurf1x2vMax_VX_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr); 
double VmLBOconstNuSurf1x2vMax_VY_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr); 

double VmLBOconstNuVol1x2vMaxP4(const double *w, const double *dxv, const double nu, const double *u, const double *vtSq, const double *f, double *out); 
double VmLBOconstNuSurf1x2vMax_VX_P4(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr); 
double VmLBOconstNuSurf1x2vMax_VY_P4(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr); 


double VmLBOconstNuVol1x3vMaxP1(const double *w, const double *dxv, const double nu, const double *u, const double *vtSq, const double *f, double *out); 
double VmLBOconstNuSurf1x3vMax_VX_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr); 
double VmLBOconstNuSurf1x3vMax_VY_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr); 
double VmLBOconstNuSurf1x3vMax_VZ_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr); 

double VmLBOconstNuVol1x3vMaxP2(const double *w, const double *dxv, const double nu, const double *u, const double *vtSq, const double *f, double *out); 
double VmLBOconstNuSurf1x3vMax_VX_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr); 
double VmLBOconstNuSurf1x3vMax_VY_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr); 
double VmLBOconstNuSurf1x3vMax_VZ_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr); 

double VmLBOconstNuVol1x3vMaxP3(const double *w, const double *dxv, const double nu, const double *u, const double *vtSq, const double *f, double *out); 
double VmLBOconstNuSurf1x3vMax_VX_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr); 
double VmLBOconstNuSurf1x3vMax_VY_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr); 
double VmLBOconstNuSurf1x3vMax_VZ_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr); 

double VmLBOconstNuVol1x3vMaxP4(const double *w, const double *dxv, const double nu, const double *u, const double *vtSq, const double *f, double *out); 
double VmLBOconstNuSurf1x3vMax_VX_P4(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr); 
double VmLBOconstNuSurf1x3vMax_VY_P4(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr); 
double VmLBOconstNuSurf1x3vMax_VZ_P4(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr); 


double VmLBOconstNuVol2x2vMaxP1(const double *w, const double *dxv, const double nu, const double *u, const double *vtSq, const double *f, double *out); 
double VmLBOconstNuSurf2x2vMax_VX_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr); 
double VmLBOconstNuSurf2x2vMax_VY_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr); 

double VmLBOconstNuVol2x2vMaxP2(const double *w, const double *dxv, const double nu, const double *u, const double *vtSq, const double *f, double *out); 
double VmLBOconstNuSurf2x2vMax_VX_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr); 
double VmLBOconstNuSurf2x2vMax_VY_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr); 

double VmLBOconstNuVol2x2vMaxP3(const double *w, const double *dxv, const double nu, const double *u, const double *vtSq, const double *f, double *out); 
double VmLBOconstNuSurf2x2vMax_VX_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr); 
double VmLBOconstNuSurf2x2vMax_VY_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr); 

double VmLBOconstNuVol2x2vMaxP4(const double *w, const double *dxv, const double nu, const double *u, const double *vtSq, const double *f, double *out); 
double VmLBOconstNuSurf2x2vMax_VX_P4(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr); 
double VmLBOconstNuSurf2x2vMax_VY_P4(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr); 


double VmLBOconstNuVol2x3vMaxP1(const double *w, const double *dxv, const double nu, const double *u, const double *vtSq, const double *f, double *out); 
double VmLBOconstNuSurf2x3vMax_VX_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr); 
double VmLBOconstNuSurf2x3vMax_VY_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr); 
double VmLBOconstNuSurf2x3vMax_VZ_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr); 

double VmLBOconstNuVol2x3vMaxP2(const double *w, const double *dxv, const double nu, const double *u, const double *vtSq, const double *f, double *out); 
double VmLBOconstNuSurf2x3vMax_VX_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr); 
double VmLBOconstNuSurf2x3vMax_VY_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr); 
double VmLBOconstNuSurf2x3vMax_VZ_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr); 

double VmLBOconstNuVol2x3vMaxP3(const double *w, const double *dxv, const double nu, const double *u, const double *vtSq, const double *f, double *out); 
double VmLBOconstNuSurf2x3vMax_VX_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr); 
double VmLBOconstNuSurf2x3vMax_VY_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr); 
double VmLBOconstNuSurf2x3vMax_VZ_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr); 

double VmLBOconstNuVol2x3vMaxP4(const double *w, const double *dxv, const double nu, const double *u, const double *vtSq, const double *f, double *out); 
double VmLBOconstNuSurf2x3vMax_VX_P4(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr); 
double VmLBOconstNuSurf2x3vMax_VY_P4(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr); 
double VmLBOconstNuSurf2x3vMax_VZ_P4(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr); 


double VmLBOconstNuVol3x3vMaxP1(const double *w, const double *dxv, const double nu, const double *u, const double *vtSq, const double *f, double *out); 
double VmLBOconstNuSurf3x3vMax_VX_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr); 
double VmLBOconstNuSurf3x3vMax_VY_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr); 
double VmLBOconstNuSurf3x3vMax_VZ_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr); 

double VmLBOconstNuVol3x3vMaxP2(const double *w, const double *dxv, const double nu, const double *u, const double *vtSq, const double *f, double *out); 
double VmLBOconstNuSurf3x3vMax_VX_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr); 
double VmLBOconstNuSurf3x3vMax_VY_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr); 
double VmLBOconstNuSurf3x3vMax_VZ_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr); 

double VmLBOconstNuVol3x3vMaxP3(const double *w, const double *dxv, const double nu, const double *u, const double *vtSq, const double *f, double *out); 
double VmLBOconstNuSurf3x3vMax_VX_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr); 
double VmLBOconstNuSurf3x3vMax_VY_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr); 
double VmLBOconstNuSurf3x3vMax_VZ_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr); 

double VmLBOconstNuVol3x3vMaxP4(const double *w, const double *dxv, const double nu, const double *u, const double *vtSq, const double *f, double *out); 
double VmLBOconstNuSurf3x3vMax_VX_P4(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr); 
double VmLBOconstNuSurf3x3vMax_VY_P4(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr); 
double VmLBOconstNuSurf3x3vMax_VZ_P4(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr); 



 
double VmLBOconstNuVol1x1vSerP1(const double *w, const double *dxv, const double nu, const double *u, const double *vtSq, const double *f, double *out); 
double VmLBOconstNuSurf1x1vSer_VX_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr); 

double VmLBOconstNuVol1x1vSerP2(const double *w, const double *dxv, const double nu, const double *u, const double *vtSq, const double *f, double *out); 
double VmLBOconstNuSurf1x1vSer_VX_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr); 

double VmLBOconstNuVol1x1vSerP3(const double *w, const double *dxv, const double nu, const double *u, const double *vtSq, const double *f, double *out); 
double VmLBOconstNuSurf1x1vSer_VX_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr); 

double VmLBOconstNuVol1x1vSerP4(const double *w, const double *dxv, const double nu, const double *u, const double *vtSq, const double *f, double *out); 
double VmLBOconstNuSurf1x1vSer_VX_P4(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr); 


double VmLBOconstNuVol1x2vSerP1(const double *w, const double *dxv, const double nu, const double *u, const double *vtSq, const double *f, double *out); 
double VmLBOconstNuSurf1x2vSer_VX_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr); 
double VmLBOconstNuSurf1x2vSer_VY_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr); 

double VmLBOconstNuVol1x2vSerP2(const double *w, const double *dxv, const double nu, const double *u, const double *vtSq, const double *f, double *out); 
double VmLBOconstNuSurf1x2vSer_VX_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr); 
double VmLBOconstNuSurf1x2vSer_VY_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr); 

double VmLBOconstNuVol1x2vSerP3(const double *w, const double *dxv, const double nu, const double *u, const double *vtSq, const double *f, double *out); 
double VmLBOconstNuSurf1x2vSer_VX_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr); 
double VmLBOconstNuSurf1x2vSer_VY_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr); 

double VmLBOconstNuVol1x2vSerP4(const double *w, const double *dxv, const double nu, const double *u, const double *vtSq, const double *f, double *out); 
double VmLBOconstNuSurf1x2vSer_VX_P4(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr); 
double VmLBOconstNuSurf1x2vSer_VY_P4(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr); 


double VmLBOconstNuVol1x3vSerP1(const double *w, const double *dxv, const double nu, const double *u, const double *vtSq, const double *f, double *out); 
double VmLBOconstNuSurf1x3vSer_VX_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr); 
double VmLBOconstNuSurf1x3vSer_VY_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr); 
double VmLBOconstNuSurf1x3vSer_VZ_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr); 

double VmLBOconstNuVol1x3vSerP2(const double *w, const double *dxv, const double nu, const double *u, const double *vtSq, const double *f, double *out); 
double VmLBOconstNuSurf1x3vSer_VX_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr); 
double VmLBOconstNuSurf1x3vSer_VY_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr); 
double VmLBOconstNuSurf1x3vSer_VZ_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr); 

double VmLBOconstNuVol1x3vSerP3(const double *w, const double *dxv, const double nu, const double *u, const double *vtSq, const double *f, double *out); 
double VmLBOconstNuSurf1x3vSer_VX_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr); 
double VmLBOconstNuSurf1x3vSer_VY_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr); 
double VmLBOconstNuSurf1x3vSer_VZ_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr); 

double VmLBOconstNuVol1x3vSerP4(const double *w, const double *dxv, const double nu, const double *u, const double *vtSq, const double *f, double *out); 
double VmLBOconstNuSurf1x3vSer_VX_P4(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr); 
double VmLBOconstNuSurf1x3vSer_VY_P4(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr); 
double VmLBOconstNuSurf1x3vSer_VZ_P4(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr); 


double VmLBOconstNuVol2x2vSerP1(const double *w, const double *dxv, const double nu, const double *u, const double *vtSq, const double *f, double *out); 
double VmLBOconstNuSurf2x2vSer_VX_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr); 
double VmLBOconstNuSurf2x2vSer_VY_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr); 

double VmLBOconstNuVol2x2vSerP2(const double *w, const double *dxv, const double nu, const double *u, const double *vtSq, const double *f, double *out); 
double VmLBOconstNuSurf2x2vSer_VX_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr); 
double VmLBOconstNuSurf2x2vSer_VY_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr); 

double VmLBOconstNuVol2x2vSerP3(const double *w, const double *dxv, const double nu, const double *u, const double *vtSq, const double *f, double *out); 
double VmLBOconstNuSurf2x2vSer_VX_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr); 
double VmLBOconstNuSurf2x2vSer_VY_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr); 

double VmLBOconstNuVol2x2vSerP4(const double *w, const double *dxv, const double nu, const double *u, const double *vtSq, const double *f, double *out); 
double VmLBOconstNuSurf2x2vSer_VX_P4(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr); 
double VmLBOconstNuSurf2x2vSer_VY_P4(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr); 


double VmLBOconstNuVol2x3vSerP1(const double *w, const double *dxv, const double nu, const double *u, const double *vtSq, const double *f, double *out); 
double VmLBOconstNuSurf2x3vSer_VX_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr); 
double VmLBOconstNuSurf2x3vSer_VY_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr); 
double VmLBOconstNuSurf2x3vSer_VZ_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr); 

double VmLBOconstNuVol2x3vSerP2(const double *w, const double *dxv, const double nu, const double *u, const double *vtSq, const double *f, double *out); 
double VmLBOconstNuSurf2x3vSer_VX_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr); 
double VmLBOconstNuSurf2x3vSer_VY_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr); 
double VmLBOconstNuSurf2x3vSer_VZ_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr); 

double VmLBOconstNuVol2x3vSerP3(const double *w, const double *dxv, const double nu, const double *u, const double *vtSq, const double *f, double *out); 
double VmLBOconstNuSurf2x3vSer_VX_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr); 
double VmLBOconstNuSurf2x3vSer_VY_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr); 
double VmLBOconstNuSurf2x3vSer_VZ_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr); 

double VmLBOconstNuVol2x3vSerP4(const double *w, const double *dxv, const double nu, const double *u, const double *vtSq, const double *f, double *out); 
double VmLBOconstNuSurf2x3vSer_VX_P4(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr); 
double VmLBOconstNuSurf2x3vSer_VY_P4(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr); 
double VmLBOconstNuSurf2x3vSer_VZ_P4(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr); 


double VmLBOconstNuVol3x3vSerP1(const double *w, const double *dxv, const double nu, const double *u, const double *vtSq, const double *f, double *out); 
double VmLBOconstNuSurf3x3vSer_VX_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr); 
double VmLBOconstNuSurf3x3vSer_VY_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr); 
double VmLBOconstNuSurf3x3vSer_VZ_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr); 

double VmLBOconstNuVol3x3vSerP2(const double *w, const double *dxv, const double nu, const double *u, const double *vtSq, const double *f, double *out); 
double VmLBOconstNuSurf3x3vSer_VX_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr); 
double VmLBOconstNuSurf3x3vSer_VY_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr); 
double VmLBOconstNuSurf3x3vSer_VZ_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr); 

double VmLBOconstNuVol3x3vSerP3(const double *w, const double *dxv, const double nu, const double *u, const double *vtSq, const double *f, double *out); 
double VmLBOconstNuSurf3x3vSer_VX_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr); 
double VmLBOconstNuSurf3x3vSer_VY_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr); 
double VmLBOconstNuSurf3x3vSer_VZ_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr); 

double VmLBOconstNuVol3x3vSerP4(const double *w, const double *dxv, const double nu, const double *u, const double *vtSq, const double *f, double *out); 
double VmLBOconstNuSurf3x3vSer_VX_P4(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr); 
double VmLBOconstNuSurf3x3vSer_VY_P4(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr); 
double VmLBOconstNuSurf3x3vSer_VZ_P4(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr); 



 
} 
#endif 
