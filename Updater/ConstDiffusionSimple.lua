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

   local ndim = self._onGrid:ndim()
   self._nu = Lin.Vec(ndim)
   local nuIn = assert(tbl.nu,
      "Updater.ConstDiffusionSimple: Must specify diffusion coefficient using 'nu'")
   local nuInType = type(nuIn)
   if (nuInType == "number") then
      for d = 1, ndim do
         self._nu[d] = nuIn
      end
   elseif (nuInType == "table") then
      assert(#nuIn == ndim, "Updater.ConstDiffusionSimple: 'nu' table entry numbers does not match dimensionality.")
      for d = 1, ndim do
         self._nu[d] = nuIn[d]
      end
   end

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
   -- Movable and re-usable pointers to underlying C data.
   local qInPtr = qIn:get(1)
   local qInPtrM = qIn:get(1)
   local qInPtrP = qIn:get(1)
   local qOutPtr = qOut:get(1)

   qOut:copy(qIn)

   -- Suggest a proper dt and quit if the dt fed in is too big.
   local dx = {}
   local nu__dx2_sum = 0
   local nu_dt__dx2 = {}
   for d = 1, ndim do
      dx[d] = grid:dx(d)
      nu__dx2_sum = nu__dx2_sum + nu[d] / (dx[d]^2)
      nu_dt__dx2[d] = dt * nu[d] / (dx[d]^2)
   end

   local dtSuggested = cfl * 0.5 / nu__dx2_sum
   local status = dt <= dtSuggested

   if not status then
      return false, dtSuggested
   end

   -- Accumulate changes along different directions.
   for d = 1, ndim do
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
        
         for icomp, c in ipairs(comps) do  -- FIXME Faster to move to inner loop?
           qOutPtr[c] = qOutPtr[c] + nu_dt__dx2[d] * (qInPtrM[c] - 2 * qInPtr[c] + qInPtrP[c])
         end
      end
   end

   return status, dtSuggested
end

function ConstDiffusionSimple:_advance(tCurr, inFld, outFld)
   return self:_forwardEuler(self, tCurr, inFld, outFld)
end

return ConstDiffusionSimple
