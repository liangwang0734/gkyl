-- Gkyl ------------------------------------------------------------------------
--
-- Updater to compute moments of distribution function on a
-- rectangular (but potentially non-uniform) grid.
--
-- If collisions (LBO for now) are included, this updater also computes the
-- boundary corrections and, if using a piecewise polynomial basis, the star
-- moments.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Gkyl libraries.
local Lin          = require "Lib.Linalg"
local LinearDecomp = require "Lib.LinearDecomp"
local MomDecl      = require "Updater.momentCalcData.DistFuncMomentCalcModDecl"
local Proto        = require "Lib.Proto"
local UpdaterBase  = require "Updater.Base"
local lume         = require "Lib.lume"
local xsys         = require "xsys"

local cudaRunTime
if GKYL_HAVE_CUDA then
   cudaRunTime = require "Cuda.RunTime"
end

-- Moments updater object.
local DistFuncMomentCalc = Proto(UpdaterBase)

-- Valid moment names for Vlasov and GK equations.
local goodMomNames = {
   "M0", "M1i", "M2ij", "M2", "M3i", "FiveMoments", "FiveMomentsLBO",
}
local goodGkMomNames = {
   "GkM0", "GkM1", "GkM1proj", "GkM2par", "GkM2perp", "GkM2", "GkM3par", "GkM3perp",
   "GkThreeMoments", "GkThreeMomentsLBO"
}
local goodPartialMomNames = {
   -- Partial velocity moments, integrating over the region
   -- where one of the velocities if positive or negative.
   "M0Pvx", "M0Pvy", "M0Pvz","M0Nvx", "M0Nvy", "M0Nvz",
   "M1iPvx", "M1iPvy", "M1iPvz","M1iNvx", "M1iNvy", "M1iNvz",
   "M2Pvx", "M2Pvy", "M2Pvz","M2Nvx", "M2Nvy", "M2Nvz",
   "M3iPvx", "M3iPvy", "M3iPvz","M3iNvx", "M3iNvy", "M3iNvz"
}

function DistFuncMomentCalc:isMomentNameGood(nm)
   if lume.find(goodMomNames, nm) then return true end
   return false
end

function DistFuncMomentCalc:isGkMomentNameGood(nm)
   if lume.find(goodGkMomNames, nm) then return true end
   return false
end

function DistFuncMomentCalc:isPartialMomentNameGood(nm)
   if lume.find(goodPartialMomNames, nm) then return true end
   return false
end

function DistFuncMomentCalc:init(tbl)
   DistFuncMomentCalc.super.init(self, tbl)    -- Setup base object.

   self._onGrid = assert(
      tbl.onGrid, "Updater.DistFuncMomentCalc: Must provide grid object using 'onGrid'")

   local phaseBasis = assert(
      tbl.phaseBasis, "Updater.DistFuncMomentCalc: Must provide phase-space basis object using 'phaseBasis'")
   local confBasis = assert(
      tbl.confBasis, "Updater.DistFuncMomentCalc: Must provide configuration-space basis object using 'confBasis'")

   self._basisID   = phaseBasis:id()
   self._polyOrder = phaseBasis:polyOrder()

   -- Dimension of spaces.
   self._pDim = phaseBasis:ndim() 
   self._cDim = confBasis:ndim()
   self._vDim = self._pDim - self._cDim

   -- Number of basis functions.
   self._numBasisP = phaseBasis:numBasis()

   -- Ensure sanity.
   assert(self._polyOrder == confBasis:polyOrder(),
	  "Polynomial orders of phase-space and config-space basis must match")
   assert(self._basisID == confBasis:id(),
	  "Type of phase-space and config-space basis must match")

   local mom = assert(
      tbl.moment, "Updater.DistFuncMomentCalc: Must provide moment to compute using 'moment'.")

   if mom == "FiveMoments" or mom == "GkThreeMoments" or
      mom == "FiveMomentsLBO" or mom == "GkThreeMomentsLBO" then
      self._fiveMoments = true
      if mom == "FiveMomentsLBO" or mom == "GkThreeMomentsLBO" then
         self._fiveMomentsLBO = true
         if mom == "FiveMomentsLBO" then    -- Rename this variable to call the right kernel below.
            mom = "FiveMoments"
         elseif mom == "GkThreeMomentsLBO" then
            mom = "GkThreeMoments"
         end
      end
   end

   local calcOnDevice = false
   if GKYL_HAVE_CUDA then
      -- This allows us to force an updater to run on the host, even for a GPU simulation.
      self._calcOnHost = tbl.onHost
      if self._calcOnHost then
         self._advanceFunc = self._advance 
      end
   else
      self._calcOnHost = true
   end

   -- Cell index, center and length right of a cell-boundary (also used for current cell for p>1).
   self.idxP = Lin.IntVec(self._pDim)
   self.xcP  = Lin.Vec(self._pDim)
   self.dxP  = Lin.Vec(self._pDim)

   if self:isPartialMomentNameGood(mom) then
      self.isPartialMom = true
      local baseMom
      for _, nm in ipairs(goodMomNames) do 
         baseMom = nm
         if string.find(mom, baseMom) then break end
      end
      -- Extract the direction and whether to integrate over the positive or negative region 
      local velTrans       = {vx=1, vy=2, vz=3}
      local partialMomReg  = string.sub(string.gsub(mom, baseMom, ""),1,1)
      self.partialMomDir   = velTrans[string.sub(string.gsub(mom, baseMom, ""),2)]
      local partialMomDirP = self._cDim+self.partialMomDir
      mom = baseMom
      -- Compute the offsets used to shorten the velocity range. For now assume that the zero
      -- along any velocity dimension is located at a cell boundary and not inside of a cell.
      self.partialMomDirExts   = {0,0}
      local partialMomDirCells = self._onGrid:numCells(partialMomDirP)
      for d = 1,self._pDim do self.idxP[d]=1 end   -- Could be any cell in other directions.
      for idx=1,partialMomDirCells do
         self.idxP[partialMomDirP] = idx
         self._onGrid:setIndex(self.idxP)
         self._onGrid:cellCenter(self.xcP)
         if (partialMomReg == "P") and (self.xcP[partialMomDirP] > 0.0) then
            self.partialMomDirExts[1] = -(idx-1)
            break
         elseif (partialMomReg == "N") and (self.xcP[partialMomDirP] > 0.0) then
            self.partialMomDirExts[2] = -(partialMomDirCells-(idx-1))
            break
         end
      end
   end

   -- Function to compute specified moment.
   self._isGk = false
   if self:isMomentNameGood(mom) then
      self._kinSpecies = "Vm"
      self._momCalcFun = MomDecl.selectMomCalc(mom, self._basisID, self._cDim, self._vDim, self._polyOrder, calcOnDevice)
   elseif self:isGkMomentNameGood(mom) then
      self._kinSpecies = "Gk"
      self._momCalcFun = MomDecl.selectGkMomCalc(mom, self._basisID, self._cDim, self._vDim, self._polyOrder)
      self._isGk       = true
      assert(tbl.gkfacs, [[DistFuncMomentCalc: must provide a gkfacs table 
                        containing the species mass and the background magnetic field
                        to calculate a Gk moment]])
   else
      assert(false, "DistFuncMomentCalc: Moments must be one of M0, M1i, M2ij, M2, M3i, FiveMoments, or FiveMomentsLBO")
   end

   self.momfac = 1.0
   if tbl.momfac then self.momfac = tbl.momfac end
   if tbl.gkfacs then
      self.mass    = tbl.gkfacs[1]
      self.bmag    = assert(tbl.gkfacs[2], "DistFuncMomentCalc: must provide bmag in gkfacs")
      self.bmagItr = self.bmag:get(1)
   end

   if self._fiveMomentsLBO then
      -- If vDim>1, intFac=2*pi/m or 4*pi/m.
      self._intFac = Lin.Vec(self._vDim)
      for d = 1,self._vDim do
        self._intFac[d] = 1.0
      end
      if self._isGk and (self._vDim > 1) then -- A (vpar,mu) simulation has 3 physical velocity dimensions.
         self._intFac[1] = 2.0*math.pi/self.mass
         self._intFac[2] = 4.0*math.pi/self.mass
      end
      self._isFirst   = true
      self._perpRange = {}    -- Perp ranges in velocity directions.
      if self._polyOrder == 1 then
         self._StarM1iM2Calc = MomDecl.selectStarM1iM2Calc(self._kinSpecies, self._basisID, self._cDim, self._vDim)
         -- Cell index, center and length left of a cell-boundary.
         self.idxM = Lin.IntVec(self._pDim)
         self.xcM  = Lin.Vec(self._pDim)
         self.dxM  = Lin.Vec(self._pDim)
      end
   end

   self.applyPositivity = xsys.pickBool(tbl.positivity,false)   -- Positivity preserving option.

   self._StarM0Calc = {} 
   self._uCorrection = {}
   self._vtSqCorrection = {}
   local uDim = self._vDim
   if (self._isGk) then uDim=1 end
   for vDir = 1, uDim do
      self._StarM0Calc[vDir] = MomDecl.selectStarM0Calc(vDir, self._kinSpecies, self._basisID, self._cDim, self._vDim, self.applyPositivity)
      self._uCorrection[vDir] = MomDecl.selectBoundaryFintegral(vDir, self._kinSpecies, self._basisID, self._cDim, self._vDim, self._polyOrder)
      self._vtSqCorrection[vDir] = MomDecl.selectBoundaryVFintegral(vDir, self._kinSpecies, self._basisID, self._cDim, self._vDim, self._polyOrder)
   end
   for vDir = 1, self._vDim do
      self._vtSqCorrection[vDir] = MomDecl.selectBoundaryVFintegral(vDir, self._kinSpecies, self._basisID, self._cDim, self._vDim, self._polyOrder)
   end

   self.onGhosts = xsys.pickBool(tbl.onGhosts, true)

   -- Option to compute moments only once per timestep, based on tCurr input parameter.
   -- NOTE: this should not be used if the updater is used to compute several different quantities in the same timestep.
   self.oncePerTime = xsys.pickBool(tbl.oncePerTime, false)

end

function DistFuncMomentCalc:initDevice(tbl)

   if not self._calcOnHost then 
      mom = tbl.moment
   
      if mom == "FiveMoments" or mom == "GkThreeMoments" or
         mom == "FiveMomentsLBO" or mom == "GkThreeMomentsLBO" then
         if mom == "FiveMomentsLBO" or mom == "GkThreeMomentsLBO" then
            if mom == "FiveMomentsLBO" then
               mom = "FiveMoments"
            elseif mom == "GkThreeMomentsLBO" then
               mom = "GkThreeMoments"
            end
         end
      end
   
      local calcOnDevice = true
      -- Select device functions/kernels.
      if not self._isGk then
         self._momCalcFun = MomDecl.selectMomCalc(mom, self._basisID, self._cDim, self._vDim, self._polyOrder, calcOnDevice)
      else
         self._momCalcFun = MomDecl.selectGkMomCalc(mom, self._basisID, self._cDim, self._vDim, self._polyOrder)
      end

   
      local vDim = self._vDim
   
      local phaseRange    = self._onGrid:localRange()
      local numCellsLocal = phaseRange:volume()
   
      local deviceNumber  = cudaRunTime.GetDevice()
      self.deviceProps    = cudaRunTime.GetDeviceProperties(deviceNumber)
   
      self.distFuncMomThreads = math.min(GKYL_DEFAULT_NUM_THREADS, phaseRange:selectLast(vDim):volume())
      self.distFuncMomBlocks  = math.floor(numCellsLocal/self.distFuncMomThreads) --+1
   end
end

-- Advance method.
function DistFuncMomentCalc:_advance(tCurr, inFld, outFld)
   if self.oncePerTime and self.tCurr == tCurr then return end -- Do nothing, already computed on this step.

   local grid        = self._onGrid
   local distf, mom1 = inFld[1], outFld[1]

   local pDim, cDim, vDim = self._pDim, self._cDim, self._vDim

   local phaseRange = distf:localRange()
   if self.onGhosts then -- Extend range to config-space ghosts.
      for dir = 1, cDim do 
         phaseRange = phaseRange:extendDir(dir, distf:lowerGhost(), distf:upperGhost())
      end
   end

   local distfItr, mom1Itr = distf:get(1), mom1:get(1)
   
   -- Construct ranges for nested loops.
   local confRangeDecomp = LinearDecomp.LinearDecompRange {
      range = phaseRange:selectFirst(cDim), numSplit = grid:numSharedProcs() }
   local velRange = phaseRange:selectLast(vDim)
   local tId      = grid:subGridSharedId()    -- Local thread ID.

   if self.isPartialMom then
      -- For partial moments, we only wish to integrate over part
      -- of one of the velocity dimensions.
      velRange = velRange:extendDir(self.partialMomDir,self.partialMomDirExts[1],self.partialMomDirExts[2])
   end
   
   local phaseIndexer = distf:genIndexer()
   local confIndexer  = mom1:genIndexer()

   local mom2, mom3
   local mom2Itr, mom3Itr
   mom1:clear(0.0) -- Zero out moments.
   
   local cMomB, cEnergyB
   local m0Star, m1Star, m2Star
   local cMomBItr, cEnergyBItr
   local m0StarItr, m1StarItr, m2StarItr
   local distfItrP, distfItrM
   if self._fiveMoments then 
      mom2 = outFld[2]
      mom3 = outFld[3] 
   
      mom2Itr = mom2:get(1)
      mom3Itr = mom3:get(1) 
   
      mom2:clear(0.0)
      mom3:clear(0.0)
      if self._fiveMomentsLBO then 
         cMomB    = outFld[4]
         cEnergyB = outFld[5] 
   
         cMomBItr    = cMomB:get(1) 
         cEnergyBItr = cEnergyB:get(1) 
   
         cMomB:clear(0.0)
         cEnergyB:clear(0.0)
   
         -- Added for corrections and star moments.
         -- Distribution functions left and right of a cell-boundary.
         distfItrP = distf:get(1)
         distfItrM = distf:get(1)
   
         if self._polyOrder == 1 then
            m0Star = outFld[6]
            m1Star = outFld[7]
            m2Star = outFld[8]
   
            m0StarItr = m0Star:get(1) 
            m1StarItr = m1Star:get(1) 
            m2StarItr = m2Star:get(1) 
   
            m0Star:clear(0.0)
            m1Star:clear(0.0)
            m2Star:clear(0.0)
         end
      end
   end
   
   -- Separate the case with and without LBO collisions to reduce number of if statements.
   if (self._fiveMomentsLBO and (self._polyOrder==1)) then 
   
      -- Outer loop is threaded and over configuration space.
      for cIdx in confRangeDecomp:rowMajorIter(tId) do
   
         cIdx:copyInto(self.idxP)
   
         mom1:fill(confIndexer(cIdx), mom1Itr)
         mom2:fill(confIndexer(cIdx), mom2Itr)
         mom3:fill(confIndexer(cIdx), mom3Itr)
   
         if self._isGk then
            self.bmag:fill(confIndexer(cIdx), self.bmagItr)
         end
   
         -- Now loop over velocity space boundary surfaces to compute boundary corrections.
         cMomB:fill(confIndexer(cIdx), cMomBItr)
         cEnergyB:fill(confIndexer(cIdx), cEnergyBItr)
   
         -- Only when the contributions to m0Star from the first direction
         -- are collected, do we collect contributions to m1Star and m2Star.
         -- Also, since Gk velocities are organized as (vpar,mu) the velocity
         -- correction is only computed for the first velocity direction.
         local firstDir = true
   
         -- isLo=true current cell is the lower boundary cell.
         -- isLo=false current cell is the upper boundary cell.
         local isLo = true
   
         -- polyOrder=1 and >1 each use separate velocity grid loops to
         -- avoid evaluating (if polyOrder==1) at each velocity coordinate.
      
         -- To have energy conservation with piece-wise linear, we must use
         -- star moments in the second equation of the weak system solved
         -- in SelfPrimMoments.
         m0Star:fill(confIndexer(cIdx), m0StarItr)
         m1Star:fill(confIndexer(cIdx), m1StarItr)
         m2Star:fill(confIndexer(cIdx), m2StarItr)
   
         for vDir = 1, vDim do
            -- Lower/upper bounds in direction 'vDir': edge indices (including outer edges).
            local dirLoIdx, dirUpIdx = phaseRange:lower(cDim+vDir), phaseRange:upper(cDim+vDir)+1
   
            if self._isFirst then
               -- Restricted velocity range.
               -- Velocity integral in m0Star does not include last cell.
               self._perpRange[vDir] = phaseRange
               for cd = 1, cDim do
                  self._perpRange[vDir] = self._perpRange[vDir]:shorten(cd) -- shorten configuration range.
               end
               self._perpRange[vDir] = self._perpRange[vDir]:shorten(cDim+vDir) -- velocity range orthogonal to 'vDir'.
            end
            local perpRange = self._perpRange[vDir]
   
            -- Outer loop is over directions orthogonal to 'vDir' and
            -- inner loop is over 1D slice in 'vDir'.
            for vPerpIdx in perpRange:rowMajorIter() do
               vPerpIdx:copyInto(self.idxM); vPerpIdx:copyInto(self.idxP)
               for d = 1, cDim do self.idxM[d] = cIdx[d] end
               for d = 1, cDim do self.idxP[d] = cIdx[d] end
      
               for i = dirLoIdx, dirUpIdx do     -- This loop is over edges.
                  self.idxM[cDim+vDir], self.idxP[cDim+vDir] = i-1, i -- Cell left/right of edge 'i'.
      
                  grid:setIndex(self.idxM)
                  grid:getDx(self.dxM)
                  grid:cellCenter(self.xcM)
      
                  grid:setIndex(self.idxP)
                  grid:getDx(self.dxP)
                  grid:cellCenter(self.xcP)
      
                  distf:fill(phaseIndexer(self.idxM), distfItrM)
                  distf:fill(phaseIndexer(self.idxP), distfItrP)
      
                  if i>dirLoIdx and i<dirUpIdx then
                     if (self._isGk) then
                        if (firstDir) then
                           self._StarM0Calc[vDir](self._intFac[1], self.xcM:data(), self.xcP:data(), self.dxM:data(), self.dxP:data(), distfItrM:data(), distfItrP:data(), m0StarItr:data())
                        end
                     else
                        self._StarM0Calc[vDir](self.xcM:data(), self.xcP:data(), self.dxM:data(), self.dxP:data(), distfItrM:data(), distfItrP:data(), m0StarItr:data())
                     end
                  end
                  if firstDir and i<dirUpIdx then
                     if self._isGk then
                        self._momCalcFun(self.xcP:data(), self.dxP:data(), self.mass, self.bmagItr:data(), distfItrP:data(), 
                         		 mom1Itr:data(), mom2Itr:data(), mom3Itr:data())
                        self._StarM1iM2Calc(self.xcP:data(), self.dxP:data(), self._intFac[1], self.mass, self.bmagItr:data(), distfItrP:data(), m1StarItr:data(), m2StarItr:data())
                     else
                        self._momCalcFun(self.xcP:data(), self.dxP:data(), distfItrP:data(), mom1Itr:data(), mom2Itr:data(), mom3Itr:data())
                        self._StarM1iM2Calc(self.xcP:data(), self.dxP:data(), distfItrP:data(), m1StarItr:data(), m2StarItr:data())
                     end
                  end
   
                  if i==dirLoIdx or i==dirUpIdx-1 then
                     local vBound = 0.0
                     -- Careful: for vBound below we assume idxP was set after idxM above.
                     if isLo then
                        vBound = grid:cellLowerInDir(cDim + vDir)
                     else
                        vBound = grid:cellUpperInDir(cDim + vDir)
                     end
                     if (self._isGk) then
                        if (firstDir) then
                           self._uCorrection[vDir](isLo, self._intFac[1], vBound, self.dxP:data(), distfItrP:data(), cMomBItr:data())
                        end
                        self._vtSqCorrection[vDir](isLo, self._intFac[vDir], vBound, self.dxP:data(), distfItrP:data(), cEnergyBItr:data())
                     else
                        self._uCorrection[vDir](isLo, vBound, self.dxP:data(), distfItrP:data(), cMomBItr:data())
                        self._vtSqCorrection[vDir](isLo, vBound, self.dxP:data(), distfItrP:data(), cEnergyBItr:data())
                     end
      
                     isLo = not isLo
                  end    -- i==dirLoIdx or i==dirUpIdx-1.
      
               end    -- Loop over edges.
            end    -- Loop over directions perpendicular to vDir.
            firstDir = false
   
         end    -- vDir loop.
      end    -- Loop over configuration space.
   
   
   else    -- if self._fiveMomentsLBO and self._polyOrder=1. 
   
      -- Outer loop is threaded and over configuration space.
      for cIdx in confRangeDecomp:rowMajorIter(tId) do
   
         cIdx:copyInto(self.idxP)
   
         -- Inner loop is over velocity space: no threading to avoid race conditions.
         for vIdx in velRange:rowMajorIter() do
            --print("vIdx = ",vIdx[1])
            for d = 1, vDim do self.idxP[cDim+d] = vIdx[d] end
         
            grid:setIndex(self.idxP)
            grid:cellCenter(self.xcP)
            grid:getDx(self.dxP)
         
            distf:fill(phaseIndexer(self.idxP), distfItr)
            mom1:fill(confIndexer(cIdx), mom1Itr)
   
            if self._isGk then
               self.bmag:fill(confIndexer(cIdx), self.bmagItr)
               if self._fiveMoments then
                  mom2:fill(confIndexer(cIdx), mom2Itr)
                  mom3:fill(confIndexer(cIdx), mom3Itr)
                  self._momCalcFun(self.xcP:data(), self.dxP:data(), self.mass, self.bmagItr:data(), distfItr:data(), 
                   		mom1Itr:data(), mom2Itr:data(), mom3Itr:data())
               else
                  self._momCalcFun(self.xcP:data(), self.dxP:data(), self.mass, self.bmagItr:data(), distfItr:data(), mom1Itr:data())
               end
            elseif self._fiveMoments then 
               mom2:fill(confIndexer(cIdx), mom2Itr)
               mom3:fill(confIndexer(cIdx), mom3Itr)
               self._momCalcFun(self.xcP:data(), self.dxP:data(), distfItr:data(), mom1Itr:data(), mom2Itr:data(), mom3Itr:data())
            else
               self._momCalcFun(self.xcP:data(), self.dxP:data(), distfItr:data(), mom1Itr:data())
            end
         end
   
         if self._fiveMomentsLBO then  -- and polyOrder>1.
            -- Now loop over velocity space boundary surfaces to compute boundary corrections.
            cMomB:fill(confIndexer(cIdx), cMomBItr)
            cEnergyB:fill(confIndexer(cIdx), cEnergyBItr)
   
            -- Only when the contributions to m0Star from the first direction
            -- are collected, do we collect contributions to m1Star and m2Star.
            -- Also, since Gk velocities are organized as (vpar,mu) the velocity
            -- correction is only computed for the first velocity direction.
            local firstDir = true
   
            -- isLo=true current cell is the lower boundary cell.
            -- isLo=false current cell is the upper boundary cell.
            local isLo = true
   
            for vDir = 1, vDim do
               -- Lower/upper bounds in direction 'vDir': cell indices.
               local dirLoIdx, dirUpIdx = phaseRange:lower(cDim+vDir), phaseRange:upper(cDim+vDir)
   
               if self._isFirst then
                  self._perpRange[vDir] = phaseRange
                  for cd = 1, cDim do
                     self._perpRange[vDir] = self._perpRange[vDir]:shorten(cd) -- shorten configuration range.
                  end
                  self._perpRange[vDir] = self._perpRange[vDir]:shorten(cDim+vDir) -- velocity range orthogonal to 'vDir'.
               end
               local perpRange = self._perpRange[vDir]
   
               for vPerpIdx in perpRange:rowMajorIter() do
                  vPerpIdx:copyInto(self.idxP)
                  for d = 1, cDim do self.idxP[d] = cIdx[d] end
      
                  for _, i in ipairs({dirLoIdx, dirUpIdx}) do     -- This loop is over edges.
                     self.idxP[cDim+vDir] = i
      
                     grid:setIndex(self.idxP)
                     grid:getDx(self.dxP)
                     grid:cellCenter(self.xcP)
      
                     distf:fill(phaseIndexer(self.idxP), distfItrP)
      
                     local vBound = 0.0
                     if isLo then
                        vBound = grid:cellLowerInDir(cDim + vDir)
                     else
                        vBound = grid:cellUpperInDir(cDim + vDir)
                     end
      
                     if (self._isGk) then
                        if (firstDir) then
                          self._uCorrection[vDir](isLo, self._intFac[1], vBound, self.dxP:data(), distfItrP:data(), cMomBItr:data())
                        end
                        self._vtSqCorrection[vDir](isLo, self._intFac[vDir], vBound, self.dxP:data(), distfItrP:data(), cEnergyBItr:data())
                     else
                        self._uCorrection[vDir](isLo, vBound, self.dxP:data(), distfItrP:data(), cMomBItr:data())
                        self._vtSqCorrection[vDir](isLo, vBound, self.dxP:data(), distfItrP:data(), cEnergyBItr:data())
                     end
      
                     isLo = not isLo
                  end
               end    -- vPerpIdx loop.
               firstDir = false
            end    -- vDir loop.
      
         end    -- if self._fiveMomentsLBO.
      end    -- Loop over configuration space.
   
   end    -- if self._fiveMomentsLBO and polyOrder=1.
   if self.momfac ~= 1.0 then mom1:scale(self.momfac) end

   if self.oncePerTime then self.tCurr = tCurr end
end

function DistFuncMomentCalc:_advanceOnDevice(tCurr, inFld, outFld)
   -- MF: current device kernels may be limited to number of velocity-space cells
   -- that are multiples of warpSize(=32), or smaller than the warpSize.
   local distf, momOut = inFld[1], outFld[1]

   momOut:deviceClear(0.0)

   self._momCalcFun(self.deviceProps, self.distFuncMomBlocks, self.distFuncMomThreads, distf._onDevice, momOut._onDevice)
end

return DistFuncMomentCalc
