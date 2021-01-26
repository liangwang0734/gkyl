-- Gkyl ------------------------------------------------------------------------
--
-- Updater to update constant-coefficient diffusion equation using simple
-- central difference for non-DG data.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local UpdaterBase = require "Updater.Base"
local Lin = require "Lib.Linalg"
local Proto = require "Lib.Proto"

local ConstDiffusionSimple = Proto(UpdaterBase)

function ConstDiffusionSimple:init(tbl)
   ConstDiffusionSimple.super.init(self, tbl) -- setup base object

   self._onGrid = assert(tbl.onGrid,
      "Updater.ConstDiffusionSimple: Must provide grid object using 'onGrid'")
   self._nu = assert(tbl.nu,
      "Updater.ConstDiffusionSimple: Must specify diffusion coefficient using 'nu'")
   self._comps = assert(tbl.components,
      "Updater.ConstDiffusionSimple: Must specify components to work on using 'components'")
   self._cfl = tbl.cfl and tbl.cfl or 1.
end

function ConstDiffusionSimple:_forwardEuler(tCurr, dt, inFld, outFld)
   local qIn = inFld[1]
   local qOut = outFld[1]

   local grid = self._onGrid
   local dt = self._dt
   local nu = self._nu
   local comps = self._comps
   local cfl = self._cfl

   local ndim = grid:ndim()

   local qInIdxr = qIn:genIndexer()
   local qOutIdxr = qOut:genIndexer()
   local idxm = Lin.IntVec(grid:ndim())
   local idxp = Lin.IntVec(grid:ndim())
   -- movable and re-usable pointers
   local qInPtr = qIn:get(1)
   local qInPtrM = qIn:get(1)
   local qInPtrP = qIn:get(1)
   local qOutPtr = qOut:get(1)

   qOut:copy(qIn)

   local status = true
   local dtSuggested = GKYL_MAX_DOUBLE

   -- accumulate changes along different directions
   for d = 1, ndim do
      local dx = grid:dx(d)
      local coeff = nu * dt / (dx*dx)

      dtSuggested = math.min(dtSuggested, cfl * 0.5 * dx * dx / nu)
      status = status and coeff <= 0.5

      local localRange = qIn:localRange()   
      for idx in localRange:rowMajorIter() do
         idx:copyInto(idxp)
         idx:copyInto(idxm)
         idxm[d] = idx[d]-1
         idxp[d] = idx[d]+1

         grid:setIndex(idx)
        
         qIn:fill(qInIdxr(idx), qInPtr)
         qIn:fill(qInIdxr(idxm), qInPtrM)
         qIn:fill(qInIdxr(idxp), qInPtrP)
        
         qOut:fill(qOutIdxr(idx), qOutPtr)
        
         for icomp, c in ipairs(comps) do  -- FIXME faster to move to inner loop?
           qOutPtr[c] = qOutPtr[c] + coeff * (qInPtrM[c] - 2 * qInPtr[c] + qInPtrP[c])
         end
      end
   end

   return status, dtSuggested
end

function ConstDiffusionSimple:_advance(tCurr, inFld, outFld)
   return self:_forwardEuler(self, tCurr, inFld, outFld)
end

return ConstDiffusionSimple
