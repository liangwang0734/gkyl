// Gkyl ------------------------------------------------------------------------
//
// C++ back-end for five-moment source terms
//    _______     ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------

#include <FiveMomentSrcImpl.h>
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <Eigen/Eigen>

// Makes indexing cleaner
static const unsigned RHO = 0;
static const unsigned MX = 1;
static const unsigned MY = 2;
static const unsigned MZ = 3;
static const unsigned ER = 4;

#define sq(x) ((x) * (x))

void
gkylFiveMomentSrcRk3(MomentSrcData_t *sd, FluidData_t *fd, double dt, double **ff, double *em, double *staticEm, double *sigma, double *auxSrc)
{
  unsigned nFluids = sd->nFluids;

  std::vector<double> keOld(nFluids);
  for (unsigned n = 0; n < nFluids; ++n)
  {
    double *f = ff[n];
    if (!fd[n].evolve)
      continue;
    keOld[n] = 0.5 * (sq(f[MX]) + sq(f[MY]) + sq(f[MZ])) / f[RHO];
  }

  // update momenta and E field
  gkylMomentSrcRk3(sd, fd, dt, ff, em, staticEm, sigma, auxSrc);

  if (sd->hasPressure)
  {
    for (unsigned n = 0; n < nFluids; ++n)
    {
      if (!fd[n].evolve)
        continue;
      double *f = ff[n];
      double keNew = 0.5 * (sq(f[MX]) + sq(f[MY]) + sq(f[MZ])) / f[RHO];
      f[ER] += keNew - keOld[n];
    }
  } 
}

void
gkylFiveMomentSrcAxisym(MomentSrcData_t *sd, FluidData_t *fd, double dt, double **ff, double *em, double *staticEm, double *sigma, double *auxSrc, double *xc)
{
  // update momenta and E field
  gkylMomentSrcAxisym(sd, fd, dt, ff, em, staticEm, sigma, auxSrc, xc);
}

void
gkylFiveMomentSrcTimeCentered(MomentSrcData_t *sd, FluidData_t *fd, double dt, double **ff, double *em, double *staticEm, double *sigma, double *auxSrc)
{
  unsigned nFluids = sd->nFluids;

  std::vector<double> keOld(nFluids);
  for (unsigned n = 0; n < nFluids; ++n)
  {
    double *f = ff[n];
    if (!fd[n].evolve)
      continue;
    keOld[n] = 0.5 * (sq(f[MX]) + sq(f[MY]) + sq(f[MZ])) / f[RHO];
  }

  // update momenta and E field
  gkylMomentSrcTimeCentered(sd, fd, dt, ff, em, staticEm, sigma, auxSrc);

  if (sd->hasPressure)
  {
    for (unsigned n = 0; n < nFluids; ++n)
    {
      if (!fd[n].evolve)
        continue;
      double *f = ff[n];
      double keNew = 0.5 * (sq(f[MX]) + sq(f[MY]) + sq(f[MZ])) / f[RHO];
      f[ER] += keNew - keOld[n];
    }
  } 
}

void
gkylFiveMomentSrcTimeCenteredDirect2(MomentSrcData_t *sd, FluidData_t *fd, double dt, double **ff, double *em, double *staticEm, double *sigma, double *auxSrc)
{
  unsigned nFluids = sd->nFluids;

  std::vector<double> keOld(nFluids);
  for (unsigned n = 0; n < nFluids; ++n)
  {
    double *f = ff[n];
    if (!fd[n].evolve)
      continue;
    keOld[n] = 0.5 * (sq(f[MX]) + sq(f[MY]) + sq(f[MZ])) / f[RHO];
  }

  // update momenta and E field
  gkylMomentSrcTimeCenteredDirect2(sd, fd, dt, ff, em, staticEm, sigma, auxSrc);

  if (sd->hasPressure)
  {
    for (unsigned n = 0; n < nFluids; ++n)
    {
      if (!fd[n].evolve)
        continue;
      double *f = ff[n];
      double keNew = 0.5 * (sq(f[MX]) + sq(f[MY]) + sq(f[MZ])) / f[RHO];
      f[ER] += keNew - keOld[n];
    }
  } 
}

void
gkylFiveMomentSrcTimeCenteredDirect(MomentSrcData_t *sd, FluidData_t *fd, double dt, double **ff, double *em, double *staticEm, double *sigma, double *auxSrc)
{
  unsigned nFluids = sd->nFluids;

  std::vector<double> keOld(nFluids);
  for (unsigned n = 0; n < nFluids; ++n)
  {
    double *f = ff[n];
    if (!fd[n].evolve)
      continue;
    keOld[n] = 0.5 * (sq(f[MX]) + sq(f[MY]) + sq(f[MZ])) / f[RHO];
  }

  // update momenta and E field
  gkylMomentSrcTimeCenteredDirect(sd, fd, dt, ff, em, staticEm, sigma, auxSrc);

  if (sd->hasPressure)
  {
    for (unsigned n = 0; n < nFluids; ++n)
    {
      if (!fd[n].evolve)
        continue;
      double *f = ff[n];
      double keNew = 0.5 * (sq(f[MX]) + sq(f[MY]) + sq(f[MZ])) / f[RHO];
      f[ER] += keNew - keOld[n];
    }
  } 
}

void
gkylFiveMomentSrcExact(MomentSrcData_t *sd, FluidData_t *fd, double dt,
                       double **ff, double *em, double *staticEm, double *sigma,
                       double *auxSrc)
{
  unsigned nFluids = sd->nFluids;

  std::vector<double> keOld(nFluids);
  for (unsigned n = 0; n < nFluids; ++n)
  {
    double *f = ff[n];
    if (!fd[n].evolve)
      continue;
    keOld[n] = 0.5 * (sq(f[MX]) + sq(f[MY]) + sq(f[MZ])) / f[RHO];
  }

  // update momenta and E field
  gkylMomentSrcExact(sd, fd, dt, ff, em, staticEm, sigma, auxSrc);

  if (sd->hasPressure)
  {
    for (unsigned n = 0; n < nFluids; ++n)
    {
      if (!fd[n].evolve)
        continue;
      double *f = ff[n];
      double keNew = 0.5 * (sq(f[MX]) + sq(f[MY]) + sq(f[MZ])) / f[RHO];
      f[ER] += keNew - keOld[n];
    }
  } 
}
