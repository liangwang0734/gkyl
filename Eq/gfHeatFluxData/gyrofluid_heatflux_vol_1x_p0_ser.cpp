#include <gyrofluid_heatflux_mod_decl.h>

double gyrofluid_heatflux_vol_1x_p0_ser(const double q_, const double m_, const double kappaPar, const double kappaPerp, const double kperpSq, const double *w, const double *dx, const double *jac, const double *jacDbmag, const double *sMom, const double *phi, double *primMom, double *out) 
{ 
  // q_,m_:   species charge and mass.
  // kappa:   heat conductivity coefficients.
  // kperpSq: k_perp^2.
  // w:       cell-center.
  // dx:      cell length.
  // jac:     jacobian.
  // jacDbmag:jacobian divided by B (J/B).
  // sMom:    stepped moments (times Jacobian).
  // phi:     electrostatic potential.
  // primMom: primitive moments (upar, Tpar, Tperp).
  // out:     output increment.

  double rdx2 = 2.0/dx[0];
  double rdxSq4 = rdx2*rdx2;



  return 1.0*fmax((jac[0]*kappaPar*m_*rdxSq4)/sMom[0],(jac[0]*kappaPerp*m_*rdxSq4)/sMom[0]); 
}
