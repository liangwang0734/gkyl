-- Gkyl ------------------------------------------------------------------------
--
-- Updater to apply Braginskii-like force and heating.
--
-- No heat conduction or viscosity.
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local UpdaterBase = require "Updater.Base"
local Lin = require "Lib.Linalg"
local Proto = require "Lib.Proto"

local BraginskiiForceAndHeating = Proto(UpdaterBase)

function BraginskiiForceAndHeating:init(tbl)
   -- setup base object
   BraginskiiForceAndHeating.super.init(self, tbl)

   local pfx = "Updater.BraginskiiForceAndHeating: "

   self._onGrid = assert(tbl.onGrid, pfx.." Must provide 'onGrid'.")

   -- Calculate tau_e with plasma parameters vs. using a preset value
   self._calcTauElectron = tbl.calcTauElectron ~= nil and tbl.calcTauElectron or false
   self._tauElectron = tbl.tauElectron
   assert(not (type(self._tauElectron)=='number' and self.calcTauElectron),
          pfx ..  "Cannot specify 'tauElectron' and 'calcTauElectron'" ..
          " simultaneously.")

   self._gasGamma = assert(tbl.gasGamma, pfx .. "Must provide 'gasGamma'.")

   self._mass = assert(tbl.mass, pfx .. "Must provide 'mass'.")

   self._charge = assert(tbl.charge, pfx .. "Must provide 'charge'.")

   self._epsilon0 = assert(tbl.epsilon0, pfx .. "Must provide 'epsilon0'.")

   self._logA = tbl.coulombLogarithm and tbl.coulombLogarithm or 10

   assert(self._gasGamma==5./3., pfx .. " gasGamma must be 5/3.")

   assert(#self._mass==2 and #self._charge==2,
          pfx .. "Lengths of mass and charge must be 2.")

   assert(self._charge[1]<0 and self._charge[1]==-self._charge[2],
          pfx .. "Electron and ion charge magnitude must equal.")
end

local temperature = function (q, gasGamma, mass)
	 local pr = (gasGamma-1)*(q[5]-0.5*(q[2]*q[2]+q[3]*q[3]+q[4]*q[4])/q[1])
    return pr * mass / q[1]
end

function BraginskiiForceAndHeating:_forwardEuler(
      self, tCurr, inFld, outFld)
   local grid = self._onGrid
   local dt = self._dt
   local gasGamma = self._gasGamma
   local me = self._mass[1]
   local mi = self._mass[2]
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
      elcBuf:fill(elcBufIdxr(idx), elcBufPtr)
      ionBuf:fill(ionBufIdxr(idx), ionBufPtr)

      elcBufPtr[5] = temperature(elcPtr, gasGamma, me)
      ionBufPtr[5] = temperature(ionPtr, gasGamma, mi)
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

      local ne = elcPtr[1] / me
      local ni = ionPtr[1] / mi
      local Te = elcBufPtr[5]
      local Ti = ionBufPtr[5]

      local bDotGradTe = bx*elcBufPtr[2] + by*elcBufPtr[3] + bz*elcBufPtr[4]
      local bDotGradTi = bx*ionBufPtr[2] + by*ionBufPtr[3] + bz*ionBufPtr[4]

      -- Thermal force.
      -- FIXME: Neglecting b x grad(T) force.
      local Rex = -0.71 * ne * bx * bDotGradTe
      local Rey = -0.71 * ne * by * bDotGradTe
      local Rez = -0.71 * ne * bz * bDotGradTe

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
      local tau_e = self._tauElectron
      if self._calcTauElectron then
         tau_e = 6*math.sqrt(2)*(math.pi^1.5)* 
                 epsilon0^2*math.sqrt(me)*Te^1.5/(logA*charge^4*ni)
      end
      local coeff = elcPtr[1] / tau_e
      Rex = Rex + coeff * (0.5 * dVparx + dVperpx)
      Rey = Rey + coeff * (0.5 * dVpary + dVperpy)
      Rez = Rez + coeff * (0.5 * dVparz + dVperpz)

      -- Heating.
      local Qi = 3 * (me/mi) * ne * (Te-Ti) / tau_e
      local Qe = -Qi + dVparx*Rex + dVpary*Rey + dVparz*Rez

      -- Final updates
      local keOld = 0.5 * (elcPtr[2]^2 + elcPtr[3]^2 + elcPtr[4]^2) / elcPtr[1]
      elcPtr[2] = elcPtr[2] + Rex
      elcPtr[3] = elcPtr[3] + Rey
      elcPtr[4] = elcPtr[4] + Rez
      local keNew = 0.5 * (elcPtr[2]^2 + elcPtr[3]^2 + elcPtr[4]^2) / elcPtr[1]
      elcPtr[5] = elcPtr[5] + keNew - keOld + Qe

      local kiOld = 0.5 * (ionPtr[2]^2 + ionPtr[3]^2 + ionPtr[4]^2) / ionPtr[1]
      ionPtr[2] = ionPtr[2] - Rex
      ionPtr[3] = ionPtr[3] - Rey
      ionPtr[4] = ionPtr[4] - Rez
      local kiNew = 0.5 * (ionPtr[2]^2 + ionPtr[3]^2 + ionPtr[4]^2) / ionPtr[1]
      ionPtr[5] = ionPtr[5] + kiNew - kiOld + Qi
   end

   return status, dtSuggested
end

function BraginskiiForceAndHeating:_advance(tCurr, inFld, outFld)
   return self:_forwardEuler(self, tCurr, inFld, outFld)
end

return BraginskiiForceAndHeating
