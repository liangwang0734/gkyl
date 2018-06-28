#include <GkLBOModDecl.h> 
double GkLBOconstNuBoundarySurf1x2vSer_Vpar_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double mufac, const double nu, const double vMuMidMax, const double *nuU, const double *nuVtSq, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[NDIM]:    Cell-center coordinates. dxv[NDIM]: Cell spacing. idx[NDIM]: current grid index.
  // nu:         constant collisionality. 
  // vMuMidMax:  maximum midpoint value of v-u. 
  // nuU[2*2]: bulk velocity (in 2 directions) times nu. 
  // nuVtSq[2]: thermal speed squared times nu. 
  // fl/fr:      Distribution function in left/right cells 
  // outl/outr:  Incremented distribution function in left/right cells 
  double rdvSq4l = 4.0/(dxvl[1]*dxvl[1]); 
  double rdvSq4r = 4.0/(dxvr[1]*dxvr[1]); 

  double incr2[8]; 
  if (idxr[1] == 1) {

    incr2[2] = nuVtSq[1]*((-1.060660171779821*fr[4])-0.6123724356957944*fr[1])+nuVtSq[0]*((-1.060660171779821*fr[2])-0.6123724356957944*fr[0]); 
    incr2[4] = nuVtSq[0]*((-1.060660171779821*fr[4])-0.6123724356957944*fr[1])+nuVtSq[1]*((-1.060660171779821*fr[2])-0.6123724356957944*fr[0]); 
    incr2[6] = nuVtSq[1]*((-1.060660171779821*fr[7])-0.6123724356957944*fr[5])+nuVtSq[0]*((-1.060660171779821*fr[6])-0.6123724356957944*fr[3]); 
    incr2[7] = nuVtSq[0]*((-1.060660171779821*fr[7])-0.6123724356957944*fr[5])+nuVtSq[1]*((-1.060660171779821*fr[6])-0.6123724356957944*fr[3]); 

    outr[2] += incr2[2]*rdvSq4r; 
    outr[4] += incr2[4]*rdvSq4r; 
    outr[6] += incr2[6]*rdvSq4r; 
    outr[7] += incr2[7]*rdvSq4r; 

  } else {

    incr2[2] = nuVtSq[1]*(0.6123724356957944*fl[1]-1.060660171779821*fl[4])+nuVtSq[0]*(0.6123724356957944*fl[0]-1.060660171779821*fl[2]); 
    incr2[4] = nuVtSq[0]*(0.6123724356957944*fl[1]-1.060660171779821*fl[4])+nuVtSq[1]*(0.6123724356957944*fl[0]-1.060660171779821*fl[2]); 
    incr2[6] = nuVtSq[1]*(0.6123724356957944*fl[5]-1.060660171779821*fl[7])+nuVtSq[0]*(0.6123724356957944*fl[3]-1.060660171779821*fl[6]); 
    incr2[7] = nuVtSq[0]*(0.6123724356957944*fl[5]-1.060660171779821*fl[7])+nuVtSq[1]*(0.6123724356957944*fl[3]-1.060660171779821*fl[6]); 

    outl[2] += incr2[2]*rdvSq4l; 
    outl[4] += incr2[4]*rdvSq4l; 
    outl[6] += incr2[6]*rdvSq4l; 
    outl[7] += incr2[7]*rdvSq4l; 

  }
  return 0.0; 
} 
double GkLBOconstNuBoundarySurf1x2vSer_Vpar_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double mufac, const double nu, const double vMuMidMax, const double *nuU, const double *nuVtSq, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[NDIM]:    Cell-center coordinates. dxv[NDIM]: Cell spacing. idx[NDIM]: current grid index.
  // nu:         constant collisionality. 
  // vMuMidMax:  maximum midpoint value of v-u. 
  // nuU[2*3]: bulk velocity (in 2 directions) times nu. 
  // nuVtSq[3]: thermal speed squared times nu. 
  // fl/fr:      Distribution function in left/right cells 
  // outl/outr:  Incremented distribution function in left/right cells 
  double rdvSq4l = 4.0/(dxvl[1]*dxvl[1]); 
  double rdvSq4r = 4.0/(dxvr[1]*dxvr[1]); 

  double incr2[20]; 
  if (idxr[1] == 1) {

    incr2[2] = nuVtSq[1]*((-1.369306393762915*fr[12])-1.060660171779821*fr[4]-0.6123724356957944*fr[1])+nuVtSq[2]*((-1.060660171779821*fr[11])-0.6123724356957944*fr[7])+nuVtSq[0]*((-1.369306393762915*fr[8])-1.060660171779821*fr[2]-0.6123724356957944*fr[0]); 
    incr2[4] = nuVtSq[2]*((-1.224744871391589*fr[12])-0.9486832980505137*fr[4]-0.5477225575051661*fr[1])+nuVtSq[0]*((-1.369306393762915*fr[12])-1.060660171779821*fr[4]-0.6123724356957944*fr[1])+nuVtSq[1]*((-0.9486832980505138*fr[11])-1.369306393762915*fr[8]-0.5477225575051661*fr[7]-1.060660171779821*fr[2]-0.6123724356957944*fr[0]); 
    incr2[6] = nuVtSq[1]*((-1.369306393762915*fr[18])-1.060660171779821*fr[10]-0.6123724356957944*fr[5])+nuVtSq[2]*((-1.060660171779821*fr[17])-0.6123724356957944*fr[13])+nuVtSq[0]*((-1.369306393762915*fr[14])-1.060660171779821*fr[6]-0.6123724356957944*fr[3]); 
    incr2[8] = nuVtSq[1]*((-5.303300858899106*fr[12])-4.107919181288745*fr[4]-2.371708245126284*fr[1])+nuVtSq[2]*((-4.107919181288745*fr[11])-2.371708245126284*fr[7])+nuVtSq[0]*((-5.303300858899105*fr[8])-4.107919181288745*fr[2]-2.371708245126284*fr[0]); 
    incr2[10] = nuVtSq[2]*((-1.224744871391589*fr[18])-0.9486832980505137*fr[10]-0.5477225575051661*fr[5])+nuVtSq[0]*((-1.369306393762915*fr[18])-1.060660171779821*fr[10]-0.6123724356957944*fr[5])+nuVtSq[1]*((-0.9486832980505137*fr[17])-1.369306393762915*fr[14]-0.5477225575051661*fr[13]-1.060660171779821*fr[6]-0.6123724356957944*fr[3]); 
    incr2[11] = nuVtSq[1]*((-1.224744871391589*fr[12])-0.9486832980505138*fr[4]-0.5477225575051661*fr[1])+nuVtSq[2]*((-0.6776309271789384*fr[11])-1.369306393762915*fr[8]-0.3912303982179757*fr[7]-1.060660171779821*fr[2]-0.6123724356957944*fr[0])+nuVtSq[0]*((-1.060660171779821*fr[11])-0.6123724356957944*fr[7]); 
    incr2[12] = nuVtSq[2]*((-4.743416490252569*fr[12])-3.674234614174767*fr[4]-2.121320343559642*fr[1])+nuVtSq[0]*((-5.303300858899105*fr[12])-4.107919181288746*fr[4]-2.371708245126284*fr[1])+nuVtSq[1]*((-3.674234614174766*fr[11])-5.303300858899106*fr[8]-2.121320343559642*fr[7]-4.107919181288746*fr[2]-2.371708245126284*fr[0]); 
    incr2[14] = nuVtSq[1]*((-5.303300858899106*fr[18])-4.107919181288745*fr[10]-2.371708245126284*fr[5])+nuVtSq[2]*((-4.107919181288745*fr[17])-2.371708245126284*fr[13])+nuVtSq[0]*((-5.303300858899105*fr[14])-4.107919181288745*fr[6]-2.371708245126284*fr[3]); 
    incr2[16] = nuVtSq[1]*((-1.060660171779821*fr[19])-0.6123724356957944*fr[15])+nuVtSq[0]*((-1.060660171779821*fr[16])-0.6123724356957944*fr[9]); 
    incr2[17] = nuVtSq[1]*((-1.224744871391589*fr[18])-0.9486832980505137*fr[10]-0.5477225575051661*fr[5])+nuVtSq[2]*((-0.6776309271789384*fr[17])-1.369306393762915*fr[14]-0.3912303982179757*fr[13]-1.060660171779821*fr[6]-0.6123724356957944*fr[3])+nuVtSq[0]*((-1.060660171779821*fr[17])-0.6123724356957944*fr[13]); 
    incr2[18] = nuVtSq[2]*((-4.743416490252569*fr[18])-3.674234614174766*fr[10]-2.121320343559642*fr[5])+nuVtSq[0]*((-5.303300858899105*fr[18])-4.107919181288745*fr[10]-2.371708245126284*fr[5])+nuVtSq[1]*((-3.674234614174766*fr[17])-5.303300858899106*fr[14]-2.121320343559642*fr[13]-4.107919181288745*fr[6]-2.371708245126284*fr[3]); 
    incr2[19] = nuVtSq[2]*((-0.9486832980505137*fr[19])-0.5477225575051661*fr[15])+nuVtSq[0]*((-1.060660171779821*fr[19])-0.6123724356957944*fr[15])+nuVtSq[1]*((-1.060660171779821*fr[16])-0.6123724356957944*fr[9]); 

    outr[2] += incr2[2]*rdvSq4r; 
    outr[4] += incr2[4]*rdvSq4r; 
    outr[6] += incr2[6]*rdvSq4r; 
    outr[8] += incr2[8]*rdvSq4r; 
    outr[10] += incr2[10]*rdvSq4r; 
    outr[11] += incr2[11]*rdvSq4r; 
    outr[12] += incr2[12]*rdvSq4r; 
    outr[14] += incr2[14]*rdvSq4r; 
    outr[16] += incr2[16]*rdvSq4r; 
    outr[17] += incr2[17]*rdvSq4r; 
    outr[18] += incr2[18]*rdvSq4r; 
    outr[19] += incr2[19]*rdvSq4r; 

  } else {

    incr2[2] = nuVtSq[1]*(1.369306393762915*fl[12]-1.060660171779821*fl[4]+0.6123724356957944*fl[1])+nuVtSq[2]*(0.6123724356957944*fl[7]-1.060660171779821*fl[11])+nuVtSq[0]*(1.369306393762915*fl[8]-1.060660171779821*fl[2]+0.6123724356957944*fl[0]); 
    incr2[4] = nuVtSq[0]*(1.369306393762915*fl[12]-1.060660171779821*fl[4]+0.6123724356957944*fl[1])+nuVtSq[2]*(1.224744871391589*fl[12]-0.9486832980505137*fl[4]+0.5477225575051661*fl[1])+nuVtSq[1]*((-0.9486832980505138*fl[11])+1.369306393762915*fl[8]+0.5477225575051661*fl[7]-1.060660171779821*fl[2]+0.6123724356957944*fl[0]); 
    incr2[6] = nuVtSq[1]*(1.369306393762915*fl[18]-1.060660171779821*fl[10]+0.6123724356957944*fl[5])+nuVtSq[2]*(0.6123724356957944*fl[13]-1.060660171779821*fl[17])+nuVtSq[0]*(1.369306393762915*fl[14]-1.060660171779821*fl[6]+0.6123724356957944*fl[3]); 
    incr2[8] = nuVtSq[1]*((-5.303300858899106*fl[12])+4.107919181288745*fl[4]-2.371708245126284*fl[1])+nuVtSq[2]*(4.107919181288745*fl[11]-2.371708245126284*fl[7])+nuVtSq[0]*((-5.303300858899105*fl[8])+4.107919181288745*fl[2]-2.371708245126284*fl[0]); 
    incr2[10] = nuVtSq[0]*(1.369306393762915*fl[18]-1.060660171779821*fl[10]+0.6123724356957944*fl[5])+nuVtSq[2]*(1.224744871391589*fl[18]-0.9486832980505137*fl[10]+0.5477225575051661*fl[5])+nuVtSq[1]*((-0.9486832980505137*fl[17])+1.369306393762915*fl[14]+0.5477225575051661*fl[13]-1.060660171779821*fl[6]+0.6123724356957944*fl[3]); 
    incr2[11] = nuVtSq[1]*(1.224744871391589*fl[12]-0.9486832980505138*fl[4]+0.5477225575051661*fl[1])+nuVtSq[2]*((-0.6776309271789384*fl[11])+1.369306393762915*fl[8]+0.3912303982179757*fl[7]-1.060660171779821*fl[2]+0.6123724356957944*fl[0])+nuVtSq[0]*(0.6123724356957944*fl[7]-1.060660171779821*fl[11]); 
    incr2[12] = nuVtSq[2]*((-4.743416490252569*fl[12])+3.674234614174767*fl[4]-2.121320343559642*fl[1])+nuVtSq[0]*((-5.303300858899105*fl[12])+4.107919181288746*fl[4]-2.371708245126284*fl[1])+nuVtSq[1]*(3.674234614174766*fl[11]-5.303300858899106*fl[8]-2.121320343559642*fl[7]+4.107919181288746*fl[2]-2.371708245126284*fl[0]); 
    incr2[14] = nuVtSq[1]*((-5.303300858899106*fl[18])+4.107919181288745*fl[10]-2.371708245126284*fl[5])+nuVtSq[2]*(4.107919181288745*fl[17]-2.371708245126284*fl[13])+nuVtSq[0]*((-5.303300858899105*fl[14])+4.107919181288745*fl[6]-2.371708245126284*fl[3]); 
    incr2[16] = nuVtSq[1]*(0.6123724356957944*fl[15]-1.060660171779821*fl[19])+nuVtSq[0]*(0.6123724356957944*fl[9]-1.060660171779821*fl[16]); 
    incr2[17] = nuVtSq[1]*(1.224744871391589*fl[18]-0.9486832980505137*fl[10]+0.5477225575051661*fl[5])+nuVtSq[2]*((-0.6776309271789384*fl[17])+1.369306393762915*fl[14]+0.3912303982179757*fl[13]-1.060660171779821*fl[6]+0.6123724356957944*fl[3])+nuVtSq[0]*(0.6123724356957944*fl[13]-1.060660171779821*fl[17]); 
    incr2[18] = nuVtSq[2]*((-4.743416490252569*fl[18])+3.674234614174766*fl[10]-2.121320343559642*fl[5])+nuVtSq[0]*((-5.303300858899105*fl[18])+4.107919181288745*fl[10]-2.371708245126284*fl[5])+nuVtSq[1]*(3.674234614174766*fl[17]-5.303300858899106*fl[14]-2.121320343559642*fl[13]+4.107919181288745*fl[6]-2.371708245126284*fl[3]); 
    incr2[19] = nuVtSq[2]*(0.5477225575051661*fl[15]-0.9486832980505137*fl[19])+nuVtSq[0]*(0.6123724356957944*fl[15]-1.060660171779821*fl[19])+nuVtSq[1]*(0.6123724356957944*fl[9]-1.060660171779821*fl[16]); 

    outl[2] += incr2[2]*rdvSq4l; 
    outl[4] += incr2[4]*rdvSq4l; 
    outl[6] += incr2[6]*rdvSq4l; 
    outl[8] += incr2[8]*rdvSq4l; 
    outl[10] += incr2[10]*rdvSq4l; 
    outl[11] += incr2[11]*rdvSq4l; 
    outl[12] += incr2[12]*rdvSq4l; 
    outl[14] += incr2[14]*rdvSq4l; 
    outl[16] += incr2[16]*rdvSq4l; 
    outl[17] += incr2[17]*rdvSq4l; 
    outl[18] += incr2[18]*rdvSq4l; 
    outl[19] += incr2[19]*rdvSq4l; 

  }
  return 0.0; 
} 
double GkLBOconstNuBoundarySurf1x2vSer_Mu_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double mufac, const double nu, const double vMuMidMax, const double *nuU, const double *nuVtSq, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[NDIM]:    Cell-center coordinates. dxv[NDIM]: Cell spacing. idx[NDIM]: current grid index.
  // nu:         constant collisionality. 
  // vMuMidMax:  maximum midpoint value of v-u. 
  // nuU[2*2]: bulk velocity (in 2 directions) times nu. 
  // nuVtSq[2]: thermal speed squared times nu. 
  // fl/fr:      Distribution function in left/right cells 
  // outl/outr:  Incremented distribution function in left/right cells 
  double rdvSq4l = 4.0/(dxvl[2]*dxvl[2]); 
  double rdvSq4r = 4.0/(dxvr[2]*dxvr[2]); 

  double incr2[8]; 
  if (idxr[2] == 1) {

    incr2[3] = (nuVtSq[1]*(((-2.121320343559642*wl[2])-1.060660171779821*dxvl[2])*fr[5]+fr[1]*((-1.224744871391589*wl[2])-0.6123724356957944*dxvl[2]))+nuVtSq[0]*(((-2.121320343559642*wl[2])-1.060660171779821*dxvl[2])*fr[3]+fr[0]*((-1.224744871391589*wl[2])-0.6123724356957944*dxvl[2])))*mufac; 
    incr2[5] = (nuVtSq[0]*(((-2.121320343559642*wl[2])-1.060660171779821*dxvl[2])*fr[5]+fr[1]*((-1.224744871391589*wl[2])-0.6123724356957944*dxvl[2]))+nuVtSq[1]*(((-2.121320343559642*wl[2])-1.060660171779821*dxvl[2])*fr[3]+fr[0]*((-1.224744871391589*wl[2])-0.6123724356957944*dxvl[2])))*mufac; 
    incr2[6] = (nuVtSq[1]*(((-2.121320343559642*wl[2])-1.060660171779821*dxvl[2])*fr[7]+((-1.224744871391589*wl[2])-0.6123724356957944*dxvl[2])*fr[4])+nuVtSq[0]*(((-2.121320343559642*wl[2])-1.060660171779821*dxvl[2])*fr[6]+fr[2]*((-1.224744871391589*wl[2])-0.6123724356957944*dxvl[2])))*mufac; 
    incr2[7] = (nuVtSq[0]*(((-2.121320343559642*wl[2])-1.060660171779821*dxvl[2])*fr[7]+((-1.224744871391589*wl[2])-0.6123724356957944*dxvl[2])*fr[4])+nuVtSq[1]*(((-2.121320343559642*wl[2])-1.060660171779821*dxvl[2])*fr[6]+fr[2]*((-1.224744871391589*wl[2])-0.6123724356957944*dxvl[2])))*mufac; 

    outr[3] += incr2[3]*rdvSq4r; 
    outr[5] += incr2[5]*rdvSq4r; 
    outr[6] += incr2[6]*rdvSq4r; 
    outr[7] += incr2[7]*rdvSq4r; 

  } else {

    incr2[3] = (nuVtSq[1]*((1.060660171779821*dxvl[2]-2.121320343559642*wl[2])*fl[5]+fl[1]*(1.224744871391589*wl[2]-0.6123724356957944*dxvl[2]))+nuVtSq[0]*((1.060660171779821*dxvl[2]-2.121320343559642*wl[2])*fl[3]+fl[0]*(1.224744871391589*wl[2]-0.6123724356957944*dxvl[2])))*mufac; 
    incr2[5] = (nuVtSq[0]*((1.060660171779821*dxvl[2]-2.121320343559642*wl[2])*fl[5]+fl[1]*(1.224744871391589*wl[2]-0.6123724356957944*dxvl[2]))+nuVtSq[1]*((1.060660171779821*dxvl[2]-2.121320343559642*wl[2])*fl[3]+fl[0]*(1.224744871391589*wl[2]-0.6123724356957944*dxvl[2])))*mufac; 
    incr2[6] = (nuVtSq[1]*((1.060660171779821*dxvl[2]-2.121320343559642*wl[2])*fl[7]+(1.224744871391589*wl[2]-0.6123724356957944*dxvl[2])*fl[4])+nuVtSq[0]*((1.060660171779821*dxvl[2]-2.121320343559642*wl[2])*fl[6]+fl[2]*(1.224744871391589*wl[2]-0.6123724356957944*dxvl[2])))*mufac; 
    incr2[7] = (nuVtSq[0]*((1.060660171779821*dxvl[2]-2.121320343559642*wl[2])*fl[7]+(1.224744871391589*wl[2]-0.6123724356957944*dxvl[2])*fl[4])+nuVtSq[1]*((1.060660171779821*dxvl[2]-2.121320343559642*wl[2])*fl[6]+fl[2]*(1.224744871391589*wl[2]-0.6123724356957944*dxvl[2])))*mufac; 

    outl[3] += incr2[3]*rdvSq4l; 
    outl[5] += incr2[5]*rdvSq4l; 
    outl[6] += incr2[6]*rdvSq4l; 
    outl[7] += incr2[7]*rdvSq4l; 

  }
  return 0.0; 
} 
double GkLBOconstNuBoundarySurf1x2vSer_Mu_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double mufac, const double nu, const double vMuMidMax, const double *nuU, const double *nuVtSq, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[NDIM]:    Cell-center coordinates. dxv[NDIM]: Cell spacing. idx[NDIM]: current grid index.
  // nu:         constant collisionality. 
  // vMuMidMax:  maximum midpoint value of v-u. 
  // nuU[2*3]: bulk velocity (in 2 directions) times nu. 
  // nuVtSq[3]: thermal speed squared times nu. 
  // fl/fr:      Distribution function in left/right cells 
  // outl/outr:  Incremented distribution function in left/right cells 
  double rdvSq4l = 4.0/(dxvl[2]*dxvl[2]); 
  double rdvSq4r = 4.0/(dxvr[2]*dxvr[2]); 

  double incr2[20]; 
  if (idxr[2] == 1) {

    incr2[3] = (nuVtSq[1]*(((-2.738612787525831*wl[2])-1.369306393762915*dxvl[2])*fr[15]+((-2.121320343559642*wl[2])-1.060660171779821*dxvl[2])*fr[5]+fr[1]*((-1.224744871391589*wl[2])-0.6123724356957944*dxvl[2]))+nuVtSq[2]*(((-2.121320343559642*wl[2])-1.060660171779821*dxvl[2])*fr[13]+((-1.224744871391589*wl[2])-0.6123724356957944*dxvl[2])*fr[7])+nuVtSq[0]*(((-2.738612787525831*wl[2])-1.369306393762915*dxvl[2])*fr[9]+((-2.121320343559642*wl[2])-1.060660171779821*dxvl[2])*fr[3]+fr[0]*((-1.224744871391589*wl[2])-0.6123724356957944*dxvl[2])))*mufac; 
    incr2[5] = (nuVtSq[2]*(((-2.449489742783178*wl[2])-1.224744871391589*dxvl[2])*fr[15]+((-1.897366596101028*wl[2])-0.9486832980505137*dxvl[2])*fr[5]+fr[1]*((-1.095445115010332*wl[2])-0.5477225575051661*dxvl[2]))+nuVtSq[0]*(((-2.738612787525831*wl[2])-1.369306393762915*dxvl[2])*fr[15]+((-2.121320343559642*wl[2])-1.060660171779821*dxvl[2])*fr[5]+fr[1]*((-1.224744871391589*wl[2])-0.6123724356957944*dxvl[2]))+nuVtSq[1]*(((-1.897366596101028*wl[2])-0.9486832980505138*dxvl[2])*fr[13]+((-2.738612787525831*wl[2])-1.369306393762915*dxvl[2])*fr[9]+((-1.095445115010332*wl[2])-0.5477225575051661*dxvl[2])*fr[7]+((-2.121320343559642*wl[2])-1.060660171779821*dxvl[2])*fr[3]+fr[0]*((-1.224744871391589*wl[2])-0.6123724356957944*dxvl[2])))*mufac; 
    incr2[6] = (nuVtSq[1]*(((-2.738612787525831*wl[2])-1.369306393762915*dxvl[2])*fr[19]+((-2.121320343559642*wl[2])-1.060660171779821*dxvl[2])*fr[10]+((-1.224744871391589*wl[2])-0.6123724356957944*dxvl[2])*fr[4])+nuVtSq[2]*(((-2.121320343559642*wl[2])-1.060660171779821*dxvl[2])*fr[17]+((-1.224744871391589*wl[2])-0.6123724356957944*dxvl[2])*fr[11])+nuVtSq[0]*(((-2.738612787525831*wl[2])-1.369306393762915*dxvl[2])*fr[16]+((-2.121320343559642*wl[2])-1.060660171779821*dxvl[2])*fr[6]+fr[2]*((-1.224744871391589*wl[2])-0.6123724356957944*dxvl[2])))*mufac; 
    incr2[9] = (nuVtSq[1]*(((-10.60660171779821*wl[2])-5.303300858899106*dxvl[2])*fr[15]+((-8.21583836257749*wl[2])-4.107919181288745*dxvl[2])*fr[5]+fr[1]*((-4.743416490252569*wl[2])-2.371708245126284*dxvl[2]))+nuVtSq[2]*(((-8.215838362577491*wl[2])-4.107919181288745*dxvl[2])*fr[13]+((-4.743416490252569*wl[2])-2.371708245126284*dxvl[2])*fr[7])+nuVtSq[0]*(((-10.60660171779821*wl[2])-5.303300858899105*dxvl[2])*fr[9]+((-8.21583836257749*wl[2])-4.107919181288745*dxvl[2])*fr[3]+fr[0]*((-4.743416490252569*wl[2])-2.371708245126284*dxvl[2])))*mufac; 
    incr2[10] = (nuVtSq[2]*(((-2.449489742783178*wl[2])-1.224744871391589*dxvl[2])*fr[19]+((-1.897366596101028*wl[2])-0.9486832980505137*dxvl[2])*fr[10]+((-1.095445115010332*wl[2])-0.5477225575051661*dxvl[2])*fr[4])+nuVtSq[0]*(((-2.738612787525831*wl[2])-1.369306393762915*dxvl[2])*fr[19]+((-2.121320343559642*wl[2])-1.060660171779821*dxvl[2])*fr[10]+((-1.224744871391589*wl[2])-0.6123724356957944*dxvl[2])*fr[4])+nuVtSq[1]*(((-1.897366596101028*wl[2])-0.9486832980505137*dxvl[2])*fr[17]+((-2.738612787525831*wl[2])-1.369306393762915*dxvl[2])*fr[16]+((-1.095445115010332*wl[2])-0.5477225575051661*dxvl[2])*fr[11]+((-2.121320343559642*wl[2])-1.060660171779821*dxvl[2])*fr[6]+fr[2]*((-1.224744871391589*wl[2])-0.6123724356957944*dxvl[2])))*mufac; 
    incr2[13] = (nuVtSq[1]*(((-2.449489742783178*wl[2])-1.224744871391589*dxvl[2])*fr[15]+((-1.897366596101028*wl[2])-0.9486832980505138*dxvl[2])*fr[5]+fr[1]*((-1.095445115010332*wl[2])-0.5477225575051661*dxvl[2]))+nuVtSq[2]*(((-1.355261854357877*wl[2])-0.6776309271789384*dxvl[2])*fr[13]+((-2.738612787525831*wl[2])-1.369306393762915*dxvl[2])*fr[9]+((-0.7824607964359517*wl[2])-0.3912303982179757*dxvl[2])*fr[7]+((-2.121320343559642*wl[2])-1.060660171779821*dxvl[2])*fr[3]+fr[0]*((-1.224744871391589*wl[2])-0.6123724356957944*dxvl[2]))+nuVtSq[0]*(((-2.121320343559642*wl[2])-1.060660171779821*dxvl[2])*fr[13]+((-1.224744871391589*wl[2])-0.6123724356957944*dxvl[2])*fr[7]))*mufac; 
    incr2[14] = (nuVtSq[1]*(((-2.121320343559642*wl[2])-1.060660171779821*dxvl[2])*fr[18]+((-1.224744871391589*wl[2])-0.6123724356957944*dxvl[2])*fr[12])+nuVtSq[0]*(((-2.121320343559642*wl[2])-1.060660171779821*dxvl[2])*fr[14]+((-1.224744871391589*wl[2])-0.6123724356957944*dxvl[2])*fr[8]))*mufac; 
    incr2[15] = (nuVtSq[2]*(((-9.48683298050514*wl[2])-4.743416490252569*dxvl[2])*fr[15]+((-7.348469228349535*wl[2])-3.674234614174767*dxvl[2])*fr[5]+fr[1]*((-4.242640687119286*wl[2])-2.121320343559642*dxvl[2]))+nuVtSq[0]*(((-10.60660171779821*wl[2])-5.303300858899105*dxvl[2])*fr[15]+((-8.215838362577493*wl[2])-4.107919181288746*dxvl[2])*fr[5]+fr[1]*((-4.743416490252569*wl[2])-2.371708245126284*dxvl[2]))+nuVtSq[1]*(((-7.348469228349534*wl[2])-3.674234614174766*dxvl[2])*fr[13]+((-10.60660171779821*wl[2])-5.303300858899106*dxvl[2])*fr[9]+((-4.242640687119286*wl[2])-2.121320343559642*dxvl[2])*fr[7]+((-8.215838362577493*wl[2])-4.107919181288746*dxvl[2])*fr[3]+fr[0]*((-4.743416490252569*wl[2])-2.371708245126284*dxvl[2])))*mufac; 
    incr2[16] = (nuVtSq[1]*(((-10.60660171779821*wl[2])-5.303300858899106*dxvl[2])*fr[19]+((-8.215838362577491*wl[2])-4.107919181288745*dxvl[2])*fr[10]+((-4.743416490252569*wl[2])-2.371708245126284*dxvl[2])*fr[4])+nuVtSq[2]*(((-8.215838362577491*wl[2])-4.107919181288745*dxvl[2])*fr[17]+((-4.743416490252569*wl[2])-2.371708245126284*dxvl[2])*fr[11])+nuVtSq[0]*(((-10.60660171779821*wl[2])-5.303300858899105*dxvl[2])*fr[16]+((-8.215838362577491*wl[2])-4.107919181288745*dxvl[2])*fr[6]+fr[2]*((-4.743416490252569*wl[2])-2.371708245126284*dxvl[2])))*mufac; 
    incr2[17] = (nuVtSq[1]*(((-2.449489742783178*wl[2])-1.224744871391589*dxvl[2])*fr[19]+((-1.897366596101028*wl[2])-0.9486832980505137*dxvl[2])*fr[10]+((-1.095445115010332*wl[2])-0.5477225575051661*dxvl[2])*fr[4])+nuVtSq[2]*(((-1.355261854357877*wl[2])-0.6776309271789384*dxvl[2])*fr[17]+((-2.738612787525831*wl[2])-1.369306393762915*dxvl[2])*fr[16]+((-0.7824607964359517*wl[2])-0.3912303982179757*dxvl[2])*fr[11]+((-2.121320343559642*wl[2])-1.060660171779821*dxvl[2])*fr[6]+fr[2]*((-1.224744871391589*wl[2])-0.6123724356957944*dxvl[2]))+nuVtSq[0]*(((-2.121320343559642*wl[2])-1.060660171779821*dxvl[2])*fr[17]+((-1.224744871391589*wl[2])-0.6123724356957944*dxvl[2])*fr[11]))*mufac; 
    incr2[18] = (nuVtSq[2]*(((-1.897366596101028*wl[2])-0.9486832980505137*dxvl[2])*fr[18]+((-1.095445115010332*wl[2])-0.5477225575051661*dxvl[2])*fr[12])+nuVtSq[0]*(((-2.121320343559642*wl[2])-1.060660171779821*dxvl[2])*fr[18]+((-1.224744871391589*wl[2])-0.6123724356957944*dxvl[2])*fr[12])+nuVtSq[1]*(((-2.121320343559642*wl[2])-1.060660171779821*dxvl[2])*fr[14]+((-1.224744871391589*wl[2])-0.6123724356957944*dxvl[2])*fr[8]))*mufac; 
    incr2[19] = (nuVtSq[2]*(((-9.48683298050514*wl[2])-4.743416490252569*dxvl[2])*fr[19]+((-7.348469228349534*wl[2])-3.674234614174766*dxvl[2])*fr[10]+((-4.242640687119286*wl[2])-2.121320343559642*dxvl[2])*fr[4])+nuVtSq[0]*(((-10.60660171779821*wl[2])-5.303300858899105*dxvl[2])*fr[19]+((-8.21583836257749*wl[2])-4.107919181288745*dxvl[2])*fr[10]+((-4.743416490252569*wl[2])-2.371708245126284*dxvl[2])*fr[4])+nuVtSq[1]*(((-7.348469228349534*wl[2])-3.674234614174766*dxvl[2])*fr[17]+((-10.60660171779821*wl[2])-5.303300858899106*dxvl[2])*fr[16]+((-4.242640687119286*wl[2])-2.121320343559642*dxvl[2])*fr[11]+((-8.21583836257749*wl[2])-4.107919181288745*dxvl[2])*fr[6]+fr[2]*((-4.743416490252569*wl[2])-2.371708245126284*dxvl[2])))*mufac; 

    outr[3] += incr2[3]*rdvSq4r; 
    outr[5] += incr2[5]*rdvSq4r; 
    outr[6] += incr2[6]*rdvSq4r; 
    outr[9] += incr2[9]*rdvSq4r; 
    outr[10] += incr2[10]*rdvSq4r; 
    outr[13] += incr2[13]*rdvSq4r; 
    outr[14] += incr2[14]*rdvSq4r; 
    outr[15] += incr2[15]*rdvSq4r; 
    outr[16] += incr2[16]*rdvSq4r; 
    outr[17] += incr2[17]*rdvSq4r; 
    outr[18] += incr2[18]*rdvSq4r; 
    outr[19] += incr2[19]*rdvSq4r; 

  } else {

    incr2[3] = (nuVtSq[1]*((2.738612787525831*wl[2]-1.369306393762915*dxvl[2])*fl[15]+(1.060660171779821*dxvl[2]-2.121320343559642*wl[2])*fl[5]+fl[1]*(1.224744871391589*wl[2]-0.6123724356957944*dxvl[2]))+nuVtSq[2]*((1.060660171779821*dxvl[2]-2.121320343559642*wl[2])*fl[13]+(1.224744871391589*wl[2]-0.6123724356957944*dxvl[2])*fl[7])+nuVtSq[0]*((2.738612787525831*wl[2]-1.369306393762915*dxvl[2])*fl[9]+(1.060660171779821*dxvl[2]-2.121320343559642*wl[2])*fl[3]+fl[0]*(1.224744871391589*wl[2]-0.6123724356957944*dxvl[2])))*mufac; 
    incr2[5] = (nuVtSq[0]*((2.738612787525831*wl[2]-1.369306393762915*dxvl[2])*fl[15]+(1.060660171779821*dxvl[2]-2.121320343559642*wl[2])*fl[5]+fl[1]*(1.224744871391589*wl[2]-0.6123724356957944*dxvl[2]))+nuVtSq[2]*((2.449489742783178*wl[2]-1.224744871391589*dxvl[2])*fl[15]+(0.9486832980505137*dxvl[2]-1.897366596101028*wl[2])*fl[5]+fl[1]*(1.095445115010332*wl[2]-0.5477225575051661*dxvl[2]))+nuVtSq[1]*((0.9486832980505138*dxvl[2]-1.897366596101028*wl[2])*fl[13]+(2.738612787525831*wl[2]-1.369306393762915*dxvl[2])*fl[9]+(1.095445115010332*wl[2]-0.5477225575051661*dxvl[2])*fl[7]+(1.060660171779821*dxvl[2]-2.121320343559642*wl[2])*fl[3]+fl[0]*(1.224744871391589*wl[2]-0.6123724356957944*dxvl[2])))*mufac; 
    incr2[6] = (nuVtSq[1]*((2.738612787525831*wl[2]-1.369306393762915*dxvl[2])*fl[19]+(1.060660171779821*dxvl[2]-2.121320343559642*wl[2])*fl[10]+(1.224744871391589*wl[2]-0.6123724356957944*dxvl[2])*fl[4])+nuVtSq[2]*((1.060660171779821*dxvl[2]-2.121320343559642*wl[2])*fl[17]+(1.224744871391589*wl[2]-0.6123724356957944*dxvl[2])*fl[11])+nuVtSq[0]*((2.738612787525831*wl[2]-1.369306393762915*dxvl[2])*fl[16]+(1.060660171779821*dxvl[2]-2.121320343559642*wl[2])*fl[6]+fl[2]*(1.224744871391589*wl[2]-0.6123724356957944*dxvl[2])))*mufac; 
    incr2[9] = (nuVtSq[1]*((5.303300858899106*dxvl[2]-10.60660171779821*wl[2])*fl[15]+(8.21583836257749*wl[2]-4.107919181288745*dxvl[2])*fl[5]+fl[1]*(2.371708245126284*dxvl[2]-4.743416490252569*wl[2]))+nuVtSq[2]*((8.215838362577491*wl[2]-4.107919181288745*dxvl[2])*fl[13]+(2.371708245126284*dxvl[2]-4.743416490252569*wl[2])*fl[7])+nuVtSq[0]*((5.303300858899105*dxvl[2]-10.60660171779821*wl[2])*fl[9]+(8.21583836257749*wl[2]-4.107919181288745*dxvl[2])*fl[3]+fl[0]*(2.371708245126284*dxvl[2]-4.743416490252569*wl[2])))*mufac; 
    incr2[10] = (nuVtSq[0]*((2.738612787525831*wl[2]-1.369306393762915*dxvl[2])*fl[19]+(1.060660171779821*dxvl[2]-2.121320343559642*wl[2])*fl[10]+(1.224744871391589*wl[2]-0.6123724356957944*dxvl[2])*fl[4])+nuVtSq[2]*((2.449489742783178*wl[2]-1.224744871391589*dxvl[2])*fl[19]+(0.9486832980505137*dxvl[2]-1.897366596101028*wl[2])*fl[10]+(1.095445115010332*wl[2]-0.5477225575051661*dxvl[2])*fl[4])+nuVtSq[1]*((0.9486832980505137*dxvl[2]-1.897366596101028*wl[2])*fl[17]+(2.738612787525831*wl[2]-1.369306393762915*dxvl[2])*fl[16]+(1.095445115010332*wl[2]-0.5477225575051661*dxvl[2])*fl[11]+(1.060660171779821*dxvl[2]-2.121320343559642*wl[2])*fl[6]+fl[2]*(1.224744871391589*wl[2]-0.6123724356957944*dxvl[2])))*mufac; 
    incr2[13] = (nuVtSq[1]*((2.449489742783178*wl[2]-1.224744871391589*dxvl[2])*fl[15]+(0.9486832980505138*dxvl[2]-1.897366596101028*wl[2])*fl[5]+fl[1]*(1.095445115010332*wl[2]-0.5477225575051661*dxvl[2]))+nuVtSq[2]*((0.6776309271789384*dxvl[2]-1.355261854357877*wl[2])*fl[13]+(2.738612787525831*wl[2]-1.369306393762915*dxvl[2])*fl[9]+(0.7824607964359517*wl[2]-0.3912303982179757*dxvl[2])*fl[7]+(1.060660171779821*dxvl[2]-2.121320343559642*wl[2])*fl[3]+fl[0]*(1.224744871391589*wl[2]-0.6123724356957944*dxvl[2]))+nuVtSq[0]*((1.060660171779821*dxvl[2]-2.121320343559642*wl[2])*fl[13]+(1.224744871391589*wl[2]-0.6123724356957944*dxvl[2])*fl[7]))*mufac; 
    incr2[14] = (nuVtSq[1]*((1.060660171779821*dxvl[2]-2.121320343559642*wl[2])*fl[18]+(1.224744871391589*wl[2]-0.6123724356957944*dxvl[2])*fl[12])+nuVtSq[0]*((1.060660171779821*dxvl[2]-2.121320343559642*wl[2])*fl[14]+(1.224744871391589*wl[2]-0.6123724356957944*dxvl[2])*fl[8]))*mufac; 
    incr2[15] = (nuVtSq[2]*((4.743416490252569*dxvl[2]-9.48683298050514*wl[2])*fl[15]+(7.348469228349535*wl[2]-3.674234614174767*dxvl[2])*fl[5]+fl[1]*(2.121320343559642*dxvl[2]-4.242640687119286*wl[2]))+nuVtSq[0]*((5.303300858899105*dxvl[2]-10.60660171779821*wl[2])*fl[15]+(8.215838362577493*wl[2]-4.107919181288746*dxvl[2])*fl[5]+fl[1]*(2.371708245126284*dxvl[2]-4.743416490252569*wl[2]))+nuVtSq[1]*((7.348469228349534*wl[2]-3.674234614174766*dxvl[2])*fl[13]+(5.303300858899106*dxvl[2]-10.60660171779821*wl[2])*fl[9]+(2.121320343559642*dxvl[2]-4.242640687119286*wl[2])*fl[7]+(8.215838362577493*wl[2]-4.107919181288746*dxvl[2])*fl[3]+fl[0]*(2.371708245126284*dxvl[2]-4.743416490252569*wl[2])))*mufac; 
    incr2[16] = (nuVtSq[1]*((5.303300858899106*dxvl[2]-10.60660171779821*wl[2])*fl[19]+(8.215838362577491*wl[2]-4.107919181288745*dxvl[2])*fl[10]+(2.371708245126284*dxvl[2]-4.743416490252569*wl[2])*fl[4])+nuVtSq[2]*((8.215838362577491*wl[2]-4.107919181288745*dxvl[2])*fl[17]+(2.371708245126284*dxvl[2]-4.743416490252569*wl[2])*fl[11])+nuVtSq[0]*((5.303300858899105*dxvl[2]-10.60660171779821*wl[2])*fl[16]+(8.215838362577491*wl[2]-4.107919181288745*dxvl[2])*fl[6]+fl[2]*(2.371708245126284*dxvl[2]-4.743416490252569*wl[2])))*mufac; 
    incr2[17] = (nuVtSq[1]*((2.449489742783178*wl[2]-1.224744871391589*dxvl[2])*fl[19]+(0.9486832980505137*dxvl[2]-1.897366596101028*wl[2])*fl[10]+(1.095445115010332*wl[2]-0.5477225575051661*dxvl[2])*fl[4])+nuVtSq[2]*((0.6776309271789384*dxvl[2]-1.355261854357877*wl[2])*fl[17]+(2.738612787525831*wl[2]-1.369306393762915*dxvl[2])*fl[16]+(0.7824607964359517*wl[2]-0.3912303982179757*dxvl[2])*fl[11]+(1.060660171779821*dxvl[2]-2.121320343559642*wl[2])*fl[6]+fl[2]*(1.224744871391589*wl[2]-0.6123724356957944*dxvl[2]))+nuVtSq[0]*((1.060660171779821*dxvl[2]-2.121320343559642*wl[2])*fl[17]+(1.224744871391589*wl[2]-0.6123724356957944*dxvl[2])*fl[11]))*mufac; 
    incr2[18] = (nuVtSq[2]*((0.9486832980505137*dxvl[2]-1.897366596101028*wl[2])*fl[18]+(1.095445115010332*wl[2]-0.5477225575051661*dxvl[2])*fl[12])+nuVtSq[0]*((1.060660171779821*dxvl[2]-2.121320343559642*wl[2])*fl[18]+(1.224744871391589*wl[2]-0.6123724356957944*dxvl[2])*fl[12])+nuVtSq[1]*((1.060660171779821*dxvl[2]-2.121320343559642*wl[2])*fl[14]+(1.224744871391589*wl[2]-0.6123724356957944*dxvl[2])*fl[8]))*mufac; 
    incr2[19] = (nuVtSq[2]*((4.743416490252569*dxvl[2]-9.48683298050514*wl[2])*fl[19]+(7.348469228349534*wl[2]-3.674234614174766*dxvl[2])*fl[10]+(2.121320343559642*dxvl[2]-4.242640687119286*wl[2])*fl[4])+nuVtSq[0]*((5.303300858899105*dxvl[2]-10.60660171779821*wl[2])*fl[19]+(8.21583836257749*wl[2]-4.107919181288745*dxvl[2])*fl[10]+(2.371708245126284*dxvl[2]-4.743416490252569*wl[2])*fl[4])+nuVtSq[1]*((7.348469228349534*wl[2]-3.674234614174766*dxvl[2])*fl[17]+(5.303300858899106*dxvl[2]-10.60660171779821*wl[2])*fl[16]+(2.121320343559642*dxvl[2]-4.242640687119286*wl[2])*fl[11]+(8.21583836257749*wl[2]-4.107919181288745*dxvl[2])*fl[6]+fl[2]*(2.371708245126284*dxvl[2]-4.743416490252569*wl[2])))*mufac; 

    outl[3] += incr2[3]*rdvSq4l; 
    outl[5] += incr2[5]*rdvSq4l; 
    outl[6] += incr2[6]*rdvSq4l; 
    outl[9] += incr2[9]*rdvSq4l; 
    outl[10] += incr2[10]*rdvSq4l; 
    outl[13] += incr2[13]*rdvSq4l; 
    outl[14] += incr2[14]*rdvSq4l; 
    outl[15] += incr2[15]*rdvSq4l; 
    outl[16] += incr2[16]*rdvSq4l; 
    outl[17] += incr2[17]*rdvSq4l; 
    outl[18] += incr2[18]*rdvSq4l; 
    outl[19] += incr2[19]*rdvSq4l; 

  }
  return 0.0; 
} 
