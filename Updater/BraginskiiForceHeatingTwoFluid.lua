-- Gkyl ------------------------------------------------------------------------
--
-- Updater to apply Braginskii-like force and heating.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local UpdaterBase = require "Updater.Base"
local Lin = require "Lib.Linalg"
local Proto = require "Lib.Proto"

local BraginskiiForceHeatingTwoFluid = Proto(UpdaterBase)

function BraginskiiForceHeatingTwoFluid:init(tbl)
   -- setup base object
   BraginskiiForceHeatingTwoFluid.super.init(self, tbl)

   self._onGrid = assert(tbl.onGrid,
      "Updater.BraginskiiForceHeatingTwoFluid:" ..
      " Must provide grid object using 'onGrid'")

   -- Calculate tau_e with plasma parameters vs. using a preset value
   self.calcTauElc = tbl.calcTauElectron ~= nil and tbl.calcTauElectron or false
   self._tauElc = tbl.tauElectron
   assert(not (type(self._tauElc)=='number' and self.calcTauElc),
      "Updater.BraginskiiForceHeatingTwoFluid:" ..
      " Cannot specify 'tauElectron' and 'calcTauElectron' simultaneously")

   self._gasGamma = assert(tbl.gasGamma,
      "Updater.BraginskiiForceHeatingTwoFluid:" ..
      " Must provide 'gasGamma'.")

   self._mass = assert(tbl.mass,
      "Updater.BraginskiiForceHeatingTwoFluid:" ..
      " Must provide 'mass'.")

   self._charge = assert(tbl.charge,
      "Updater.BraginskiiForceHeatingTwoFluid:" ..
      " Must provide 'charge'.")

   self._epsilon0 = assert(tbl.epsilon0,
      "Updater.BraginskiiForceHeatingTwoFluid:" ..
      " Must provide 'epsilon0'.")

   self._logA = tbl.coulombLogarithm and tbl.coulombLogarithm or 10

   assert(self._gasGamma==5./3.,
      "Updater.BraginskiiForceHeatingTwoFluid:" ..
      " gasGamma must be 5/3.")

   assert(#self._mass==2 and #self._charge==2,
      "Updater.BraginskiiForceHeatingTwoFluid:" ..
      " Only two-fluid electron-ion plasmas are supprted.")

   assert(self._charge[1]<0 and self._charge[1]==-self._charge[2],
      "Updater.BraginskiiForceHeatingTwoFluid:" ..
      " Electron and ion charge magnitude must equal.")
end

local temperature = function (q, gasGamma, mass)
	 local pr = (gasGamma-1)*(q[5]-0.5*(q[2]*q[2]+q[3]*q[3]+q[4]*q[4])/q[1])
    return pr * mass / q[1]
end

function BraginskiiForceHeatingTwoFluid:_forwardEuler(
      self, tCurr, inFld, outFld)
   local grid = self._onGrid
   local dt = self._dt
   local gasGamma = self._gasGamma
   local elcMass = self._mass[1]
   local ionMass = self._mass[2]
   local charge = self._charge[2]
   local epsilon0 = self._epsilon0
   local logA = self._logA

   local dtSuggested = GKYL_MAX_DOUBLE
   local status = true

   -- Field objects.
   local elc = outFld[1]
   local ion = outFld[2]
   local emf = outFld[3]

   local elcBuf = inFld[1]
   local ionBuf = inFld[2]

   -- Field indexers.
   local elcIdxr = elc:genIndexer()
   local ionIdxr = ion:genIndexer()
   local emfIdxr = emf:genIndexer()

   local elcBufIdxr = elcBuf:genIndexer()
   local ionBufIdxr = ionBuf:genIndexer()

   -- Movable and re-usable field pointers to underlying C data.
   local elcPtr = elc:get(1)
   local ionPtr = ion:get(1)
   local emfPtr = emf:get(1)

   local elcBufPtr = elcBuf:get(1)
   local elcBufPtrM = elcBuf:get(1)
   local elcBufPtrP = elcBuf:get(1)

   local ionBufPtr = ionBuf:get(1)
   local ionBufPtrM = ionBuf:get(1)
   local ionBufPtrP = ionBuf:get(1)

   -- Indices.
   local ndim = grid:ndim()
   local idxm = Lin.IntVec(grid:ndim())
   local idxp = Lin.IntVec(grid:ndim())

   -- Compute temperature in internal and ghost cells.
   local localExtRange = elc:localExtRange()
   for idx in localExtRange:rowMajorIter() do
      elc:fill(elcIdxr(idx), elcPtr)
      ion:fill(ionIdxr(idx), ionPtr)
      elcBuf:fill(elcIdxr(idx), elcBufPtr)
      ionBuf:fill(ionIdxr(idx), ionBufPtr)

      elcBufPtr[5] = temperature(elcPtr, gasGamma, elcMass)
      ionBufPtr[5] = temperature(ionPtr, gasGamma, ionMass)
   end

   -- Comptue grad(T) in internal cells.
   local localRange = elc:localRange()
   for d = 1, ndim do
      for idx in localRange:rowMajorIter() do
         idx:copyInto(idxp)
         idx:copyInto(idxm)
         idxm[d] = idx[d]-1
         idxp[d] = idx[d]+1

         local __2dx = 0.5 / grid:dx(d)
        
         elcBuf:fill(elcBufIdxr(idx), elcBufPtr)
         elcBuf:fill(elcBufIdxr(idxm), elcBufPtrM)
         elcBuf:fill(elcBufIdxr(idxp), elcBufPtrP)
         elcBufPtr[d+1] = (elcBufPtrP[5] - elcBufPtrM[5]) * __2dx
        
         ionBuf:fill(ionBufIdxr(idx), ionBufPtr)
         ionBuf:fill(ionBufIdxr(idxm), ionBufPtrM)
         ionBuf:fill(ionBufIdxr(idxp), ionBufPtrP)
         ionBufPtr[d+1] = (ionBufPtrP[5] - ionBufPtrM[5]) * __2dx
      end
   end

   for idx in localRange:rowMajorIter() do
      elc:fill(elcIdxr(idx), elcPtr)
      ion:fill(ionIdxr(idx), ionPtr)
      emf:fill(emfIdxr(idx), emfPtr)

      -- Buffers that contain grad(T) and will store force and heating.
      elcBuf:fill(elcBufIdxr(idx), elcBufPtr)
      ionBuf:fill(ionBufIdxr(idx), ionBufPtr)

      local bx = emfPtr[4]
      local by = emfPtr[5]
      local bz = emfPtr[6]
      local bmag = math.sqrt(bx*bx + by*by + bz*bz)
      bx = bx / bmag
      by = by / bmag
      bz = bz / bmag

      local nElc = elcPtr[1] / elcMass
      local nIon = ionPtr[1] / ionMass
      local TElc = elcBufPtr[5]
      local TIon = ionBufPtr[5]

      local bDotGradTe = bx * elcBufPtr[2] + by * elcBufPtr[3] + bz * elcBufPtr[4]
      local bDotGradTi = bx * ionBufPtr[2] + by * ionBufPtr[3] + bz * ionBufPtr[4]

      -- Thermal force.
      -- FIXME: Neglecting b x grad(T) force.
      local RElcx = -0.71 * nElc * bx * bDotGradTe
      local RElcy = -0.71 * nElc * by * bDotGradTe
      local RElcz = -0.71 * nElc * bz * bDotGradTe

      local dVx = ionPtr[2]/ionPtr[1] - elcPtr[2]/elcPtr[1]
      local dVy = ionPtr[3]/ionPtr[1] - elcPtr[3]/elcPtr[1]
      local dVz = ionPtr[4]/ionPtr[1] - elcPtr[4]/elcPtr[1]
      local bDotDV = bx * dVx + by * dVy + bz * dVz
      local dVparx = bx * bDotDV
      local dVpary = by * bDotDV
      local dVparz = bz * bDotDV
      local dVperpx = dVx - dVparx
      local dVperpy = dVy - dVpary
      local dVperpz = dVz - dVparz

      -- Friction force.
      local tauElc = self._tauElc
      if self.calcTauElc then
         tauElc = 6*math.sqrt(2)*(math.pi^1.5)* 
                  epsilon0^2*math.sqrt(elcMass)*TElc^1.5/(logA*charge^4*nIon)
      end
      local coeff = elcPtr[1] / tauElc
      RElcx = RElcx + coeff * (0.5 * dVparx + dVperpx)
      RElcy = RElcy + coeff * (0.5 * dVpary + dVperpy)
      RElcz = RElcz + coeff * (0.5 * dVparz + dVperpz)

      -- Heating.
      local QIon = 3 * (elcMass/ionMass) * nElc * (TElc-TIon) / tauElc
      local QElc = -QIon + dVparx*RElcx + dVpary*RElcy + dVparz*RElcz

      -- Final updates
      local kElcOld = 0.5 * (elcPtr[2]^2 + elcPtr[3]^2 + elcPtr[4]^2) / elcPtr[1]
      elcPtr[2] = elcPtr[2] + RElcx
      elcPtr[3] = elcPtr[3] + RElcy
      elcPtr[4] = elcPtr[4] + RElcz
      local kElcNew = 0.5 * (elcPtr[2]^2 + elcPtr[3]^2 + elcPtr[4]^2) / elcPtr[1]
      elcPtr[5] = elcPtr[5] + kElcNew - kElcOld + QElc

      local kIonOld = 0.5 * (ionPtr[2]^2 + ionPtr[3]^2 + ionPtr[4]^2) / ionPtr[1]
      ionPtr[2] = ionPtr[2] - RElcx
      ionPtr[3] = ionPtr[3] - RElcy
      ionPtr[4] = ionPtr[4] - RElcz
      local kIonNew = 0.5 * (ionPtr[2]^2 + ionPtr[3]^2 + ionPtr[4]^2) / ionPtr[1]
      ionPtr[5] = ionPtr[5] + kIonNew - kIonOld + QIon
   end

   return status, dtSuggested
end

function BraginskiiForceHeatingTwoFluid:_advance(tCurr, inFld, outFld)
   return self:_forwardEuler(self, tCurr, inFld, outFld)
end

return BraginskiiForceHeatingTwoFluid
