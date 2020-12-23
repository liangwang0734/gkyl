#include <ConstDiffusionModDecl.h> 
double ConstDiffusionVol3xTensorP2_diffDirs1(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[3]:   Cell-center coordinates.
  // dx[3]:  Cell spacing.
  // nu[3]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[1]; 
  rdxFnu[0] = 4.0*nu[0]/(dx[0]*dx[0]); 

  out[7] += 6.708203932499369*f[0]*rdxFnu[0]; 
  out[11] += 6.708203932499369*rdxFnu[0]*f[2]; 
  out[13] += 6.708203932499369*rdxFnu[0]*f[3]; 
  out[17] += 6.708203932499369*rdxFnu[0]*f[6]; 
  out[20] += 6.708203932499369*rdxFnu[0]*f[8]; 
  out[21] += 6.708203932499369*rdxFnu[0]*f[9]; 
  out[23] += 6.708203932499369*rdxFnu[0]*f[14]; 
  out[24] += 6.708203932499369*rdxFnu[0]*f[16]; 
  out[26] += 6.708203932499369*rdxFnu[0]*f[22]; 

  return (rdxFnu[0])*0.9;

} 
double ConstDiffusionVol3xTensorP2_diffDirs12(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[3]:   Cell-center coordinates.
  // dx[3]:  Cell spacing.
  // nu[3]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[2]; 
  rdxFnu[0] = 4.0*nu[0]/(dx[0]*dx[0]); 
  rdxFnu[1] = 4.0*nu[1]/(dx[1]*dx[1]); 

  out[7] += 6.708203932499369*f[0]*rdxFnu[0]; 
  out[8] += 6.708203932499369*f[0]*rdxFnu[1]; 
  out[11] += 6.708203932499369*rdxFnu[0]*f[2]; 
  out[12] += 6.708203932499369*f[1]*rdxFnu[1]; 
  out[13] += 6.708203932499369*rdxFnu[0]*f[3]; 
  out[14] += 6.708203932499369*rdxFnu[1]*f[3]; 
  out[17] += 6.708203932499369*rdxFnu[0]*f[6]; 
  out[18] += 6.708203932499369*rdxFnu[1]*f[5]; 
  out[20] += 6.708203932499369*rdxFnu[0]*f[8]+6.708203932499369*rdxFnu[1]*f[7]; 
  out[21] += 6.708203932499369*rdxFnu[0]*f[9]; 
  out[22] += 6.708203932499369*rdxFnu[1]*f[9]; 
  out[23] += 6.708203932499369*rdxFnu[0]*f[14]+6.708203932499369*rdxFnu[1]*f[13]; 
  out[24] += 6.708203932499369*rdxFnu[0]*f[16]; 
  out[25] += 6.708203932499369*rdxFnu[1]*f[15]; 
  out[26] += 6.708203932499369*rdxFnu[0]*f[22]+6.708203932499369*rdxFnu[1]*f[21]; 

  return (rdxFnu[0]+rdxFnu[1])*0.9;

} 
double ConstDiffusionVol3xTensorP2_diffDirs123(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[3]:   Cell-center coordinates.
  // dx[3]:  Cell spacing.
  // nu[3]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[3]; 
  rdxFnu[0] = 4.0*nu[0]/(dx[0]*dx[0]); 
  rdxFnu[1] = 4.0*nu[1]/(dx[1]*dx[1]); 
  rdxFnu[2] = 4.0*nu[2]/(dx[2]*dx[2]); 

  out[7] += 6.708203932499369*f[0]*rdxFnu[0]; 
  out[8] += 6.708203932499369*f[0]*rdxFnu[1]; 
  out[9] += 6.708203932499369*f[0]*rdxFnu[2]; 
  out[11] += 6.708203932499369*rdxFnu[0]*f[2]; 
  out[12] += 6.708203932499369*f[1]*rdxFnu[1]; 
  out[13] += 6.708203932499369*rdxFnu[0]*f[3]; 
  out[14] += 6.708203932499369*rdxFnu[1]*f[3]; 
  out[15] += 6.708203932499369*f[1]*rdxFnu[2]; 
  out[16] += 6.708203932499369*f[2]*rdxFnu[2]; 
  out[17] += 6.708203932499369*rdxFnu[0]*f[6]; 
  out[18] += 6.708203932499369*rdxFnu[1]*f[5]; 
  out[19] += 6.708203932499369*rdxFnu[2]*f[4]; 
  out[20] += 6.708203932499369*rdxFnu[0]*f[8]+6.708203932499369*rdxFnu[1]*f[7]; 
  out[21] += 6.708203932499369*rdxFnu[0]*f[9]+6.708203932499369*rdxFnu[2]*f[7]; 
  out[22] += 6.708203932499369*rdxFnu[1]*f[9]+6.708203932499369*rdxFnu[2]*f[8]; 
  out[23] += 6.708203932499369*rdxFnu[0]*f[14]+6.708203932499369*rdxFnu[1]*f[13]; 
  out[24] += 6.708203932499369*rdxFnu[0]*f[16]+6.708203932499369*rdxFnu[2]*f[11]; 
  out[25] += 6.708203932499369*rdxFnu[1]*f[15]+6.708203932499369*rdxFnu[2]*f[12]; 
  out[26] += 6.708203932499369*rdxFnu[0]*f[22]+6.708203932499369*rdxFnu[1]*f[21]+6.708203932499369*rdxFnu[2]*f[20]; 

  return (rdxFnu[0]+rdxFnu[1]+rdxFnu[2])*0.9;

} 
double ConstDiffusionVol3xTensorP2_diffDirs13(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[3]:   Cell-center coordinates.
  // dx[3]:  Cell spacing.
  // nu[3]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[2]; 
  rdxFnu[0] = 4.0*nu[0]/(dx[0]*dx[0]); 
  rdxFnu[1] = 4.0*nu[2]/(dx[2]*dx[2]); 

  out[7] += 6.708203932499369*f[0]*rdxFnu[0]; 
  out[9] += 6.708203932499369*f[0]*rdxFnu[1]; 
  out[11] += 6.708203932499369*rdxFnu[0]*f[2]; 
  out[13] += 6.708203932499369*rdxFnu[0]*f[3]; 
  out[15] += 6.708203932499369*f[1]*rdxFnu[1]; 
  out[16] += 6.708203932499369*rdxFnu[1]*f[2]; 
  out[17] += 6.708203932499369*rdxFnu[0]*f[6]; 
  out[19] += 6.708203932499369*rdxFnu[1]*f[4]; 
  out[20] += 6.708203932499369*rdxFnu[0]*f[8]; 
  out[21] += 6.708203932499369*rdxFnu[0]*f[9]+6.708203932499369*rdxFnu[1]*f[7]; 
  out[22] += 6.708203932499369*rdxFnu[1]*f[8]; 
  out[23] += 6.708203932499369*rdxFnu[0]*f[14]; 
  out[24] += 6.708203932499369*rdxFnu[0]*f[16]+6.708203932499369*rdxFnu[1]*f[11]; 
  out[25] += 6.708203932499369*rdxFnu[1]*f[12]; 
  out[26] += 6.708203932499369*rdxFnu[0]*f[22]+6.708203932499369*rdxFnu[1]*f[20]; 

  return (rdxFnu[0]+rdxFnu[1])*0.9;

} 
double ConstDiffusionVol3xTensorP2_diffDirs2(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[3]:   Cell-center coordinates.
  // dx[3]:  Cell spacing.
  // nu[3]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[1]; 
  rdxFnu[0] = 4.0*nu[1]/(dx[1]*dx[1]); 

  out[8] += 6.708203932499369*f[0]*rdxFnu[0]; 
  out[12] += 6.708203932499369*rdxFnu[0]*f[1]; 
  out[14] += 6.708203932499369*rdxFnu[0]*f[3]; 
  out[18] += 6.708203932499369*rdxFnu[0]*f[5]; 
  out[20] += 6.708203932499369*rdxFnu[0]*f[7]; 
  out[22] += 6.708203932499369*rdxFnu[0]*f[9]; 
  out[23] += 6.708203932499369*rdxFnu[0]*f[13]; 
  out[25] += 6.708203932499369*rdxFnu[0]*f[15]; 
  out[26] += 6.708203932499369*rdxFnu[0]*f[21]; 

  return (rdxFnu[0])*0.9;

} 
double ConstDiffusionVol3xTensorP2_diffDirs23(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[3]:   Cell-center coordinates.
  // dx[3]:  Cell spacing.
  // nu[3]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[2]; 
  rdxFnu[0] = 4.0*nu[1]/(dx[1]*dx[1]); 
  rdxFnu[1] = 4.0*nu[2]/(dx[2]*dx[2]); 

  out[8] += 6.708203932499369*f[0]*rdxFnu[0]; 
  out[9] += 6.708203932499369*f[0]*rdxFnu[1]; 
  out[12] += 6.708203932499369*rdxFnu[0]*f[1]; 
  out[14] += 6.708203932499369*rdxFnu[0]*f[3]; 
  out[15] += 6.708203932499369*f[1]*rdxFnu[1]; 
  out[16] += 6.708203932499369*rdxFnu[1]*f[2]; 
  out[18] += 6.708203932499369*rdxFnu[0]*f[5]; 
  out[19] += 6.708203932499369*rdxFnu[1]*f[4]; 
  out[20] += 6.708203932499369*rdxFnu[0]*f[7]; 
  out[21] += 6.708203932499369*rdxFnu[1]*f[7]; 
  out[22] += 6.708203932499369*rdxFnu[0]*f[9]+6.708203932499369*rdxFnu[1]*f[8]; 
  out[23] += 6.708203932499369*rdxFnu[0]*f[13]; 
  out[24] += 6.708203932499369*rdxFnu[1]*f[11]; 
  out[25] += 6.708203932499369*rdxFnu[0]*f[15]+6.708203932499369*rdxFnu[1]*f[12]; 
  out[26] += 6.708203932499369*rdxFnu[0]*f[21]+6.708203932499369*rdxFnu[1]*f[20]; 

  return (rdxFnu[0]+rdxFnu[1])*0.9;

} 
double ConstDiffusionVol3xTensorP2_diffDirs3(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[3]:   Cell-center coordinates.
  // dx[3]:  Cell spacing.
  // nu[3]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[1]; 
  rdxFnu[0] = 4.0*nu[2]/(dx[2]*dx[2]); 

  out[9] += 6.708203932499369*f[0]*rdxFnu[0]; 
  out[15] += 6.708203932499369*rdxFnu[0]*f[1]; 
  out[16] += 6.708203932499369*rdxFnu[0]*f[2]; 
  out[19] += 6.708203932499369*rdxFnu[0]*f[4]; 
  out[21] += 6.708203932499369*rdxFnu[0]*f[7]; 
  out[22] += 6.708203932499369*rdxFnu[0]*f[8]; 
  out[24] += 6.708203932499369*rdxFnu[0]*f[11]; 
  out[25] += 6.708203932499369*rdxFnu[0]*f[12]; 
  out[26] += 6.708203932499369*rdxFnu[0]*f[20]; 

  return (rdxFnu[0])*0.9;

} 
double ConstHyperDiffusion4Vol3xTensorP2_diffDirs1(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[3]:   Cell-center coordinates.
  // dx[3]:  Cell spacing.
  // nu[3]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[1]; 
  rdxFnu[0] = 16.0*nu[0]/(dx[0]*dx[0]*dx[0]*dx[0]); 


  return (rdxFnu[0])*0.9;

} 
double ConstHyperDiffusion4Vol3xTensorP2_diffDirs12(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[3]:   Cell-center coordinates.
  // dx[3]:  Cell spacing.
  // nu[3]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[2]; 
  rdxFnu[0] = 16.0*nu[0]/(dx[0]*dx[0]*dx[0]*dx[0]); 
  rdxFnu[1] = 16.0*nu[1]/(dx[1]*dx[1]*dx[1]*dx[1]); 


  return (rdxFnu[0]+rdxFnu[1])*0.9;

} 
double ConstHyperDiffusion4Vol3xTensorP2_diffDirs123(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[3]:   Cell-center coordinates.
  // dx[3]:  Cell spacing.
  // nu[3]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[3]; 
  rdxFnu[0] = 16.0*nu[0]/(dx[0]*dx[0]*dx[0]*dx[0]); 
  rdxFnu[1] = 16.0*nu[1]/(dx[1]*dx[1]*dx[1]*dx[1]); 
  rdxFnu[2] = 16.0*nu[2]/(dx[2]*dx[2]*dx[2]*dx[2]); 


  return (rdxFnu[0]+rdxFnu[1]+rdxFnu[2])*0.9;

} 
double ConstHyperDiffusion4Vol3xTensorP2_diffDirs13(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[3]:   Cell-center coordinates.
  // dx[3]:  Cell spacing.
  // nu[3]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[2]; 
  rdxFnu[0] = 16.0*nu[0]/(dx[0]*dx[0]*dx[0]*dx[0]); 
  rdxFnu[1] = 16.0*nu[2]/(dx[2]*dx[2]*dx[2]*dx[2]); 


  return (rdxFnu[0]+rdxFnu[1])*0.9;

} 
double ConstHyperDiffusion4Vol3xTensorP2_diffDirs2(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[3]:   Cell-center coordinates.
  // dx[3]:  Cell spacing.
  // nu[3]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[1]; 
  rdxFnu[0] = 16.0*nu[1]/(dx[1]*dx[1]*dx[1]*dx[1]); 


  return (rdxFnu[0])*0.9;

} 
double ConstHyperDiffusion4Vol3xTensorP2_diffDirs23(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[3]:   Cell-center coordinates.
  // dx[3]:  Cell spacing.
  // nu[3]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[2]; 
  rdxFnu[0] = 16.0*nu[1]/(dx[1]*dx[1]*dx[1]*dx[1]); 
  rdxFnu[1] = 16.0*nu[2]/(dx[2]*dx[2]*dx[2]*dx[2]); 


  return (rdxFnu[0]+rdxFnu[1])*0.9;

} 
double ConstHyperDiffusion4Vol3xTensorP2_diffDirs3(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[3]:   Cell-center coordinates.
  // dx[3]:  Cell spacing.
  // nu[3]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[1]; 
  rdxFnu[0] = 16.0*nu[2]/(dx[2]*dx[2]*dx[2]*dx[2]); 


  return (rdxFnu[0])*0.9;

} 
double ConstHyperDiffusion6Vol3xTensorP2_diffDirs1(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[3]:   Cell-center coordinates.
  // dx[3]:  Cell spacing.
  // nu[3]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[1]; 
  rdxFnu[0] = 64.0*nu[0]/(dx[0]*dx[0]*dx[0]*dx[0]*dx[0]*dx[0]); 


  return (rdxFnu[0])*0.9;

} 
double ConstHyperDiffusion6Vol3xTensorP2_diffDirs12(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[3]:   Cell-center coordinates.
  // dx[3]:  Cell spacing.
  // nu[3]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[2]; 
  rdxFnu[0] = 64.0*nu[0]/(dx[0]*dx[0]*dx[0]*dx[0]*dx[0]*dx[0]); 
  rdxFnu[1] = 64.0*nu[1]/(dx[1]*dx[1]*dx[1]*dx[1]*dx[1]*dx[1]); 


  return (rdxFnu[0]+rdxFnu[1])*0.9;

} 
double ConstHyperDiffusion6Vol3xTensorP2_diffDirs123(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[3]:   Cell-center coordinates.
  // dx[3]:  Cell spacing.
  // nu[3]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[3]; 
  rdxFnu[0] = 64.0*nu[0]/(dx[0]*dx[0]*dx[0]*dx[0]*dx[0]*dx[0]); 
  rdxFnu[1] = 64.0*nu[1]/(dx[1]*dx[1]*dx[1]*dx[1]*dx[1]*dx[1]); 
  rdxFnu[2] = 64.0*nu[2]/(dx[2]*dx[2]*dx[2]*dx[2]*dx[2]*dx[2]); 


  return (rdxFnu[0]+rdxFnu[1]+rdxFnu[2])*0.9;

} 
double ConstHyperDiffusion6Vol3xTensorP2_diffDirs13(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[3]:   Cell-center coordinates.
  // dx[3]:  Cell spacing.
  // nu[3]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[2]; 
  rdxFnu[0] = 64.0*nu[0]/(dx[0]*dx[0]*dx[0]*dx[0]*dx[0]*dx[0]); 
  rdxFnu[1] = 64.0*nu[2]/(dx[2]*dx[2]*dx[2]*dx[2]*dx[2]*dx[2]); 


  return (rdxFnu[0]+rdxFnu[1])*0.9;

} 
double ConstHyperDiffusion6Vol3xTensorP2_diffDirs2(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[3]:   Cell-center coordinates.
  // dx[3]:  Cell spacing.
  // nu[3]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[1]; 
  rdxFnu[0] = 64.0*nu[1]/(dx[1]*dx[1]*dx[1]*dx[1]*dx[1]*dx[1]); 


  return (rdxFnu[0])*0.9;

} 
double ConstHyperDiffusion6Vol3xTensorP2_diffDirs23(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[3]:   Cell-center coordinates.
  // dx[3]:  Cell spacing.
  // nu[3]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[2]; 
  rdxFnu[0] = 64.0*nu[1]/(dx[1]*dx[1]*dx[1]*dx[1]*dx[1]*dx[1]); 
  rdxFnu[1] = 64.0*nu[2]/(dx[2]*dx[2]*dx[2]*dx[2]*dx[2]*dx[2]); 


  return (rdxFnu[0]+rdxFnu[1])*0.9;

} 
double ConstHyperDiffusion6Vol3xTensorP2_diffDirs3(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[3]:   Cell-center coordinates.
  // dx[3]:  Cell spacing.
  // nu[3]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[1]; 
  rdxFnu[0] = 64.0*nu[2]/(dx[2]*dx[2]*dx[2]*dx[2]*dx[2]*dx[2]); 


  return (rdxFnu[0])*0.9;

} 
