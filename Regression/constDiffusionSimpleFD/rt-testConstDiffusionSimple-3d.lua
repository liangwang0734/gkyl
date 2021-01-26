-- Gkyl ------------------------------------------------------------------------
--
-- Test the 3D ConstDiffusion equation object to solve a diffusion equation
--    u_t = D_x*u_{xx}+D_y*u_{yy}+D_z*u_{zz}.
-- with constant diffusion coefficient D.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- PREAMBLE.
local Grid       = require "Grid"
local DataStruct = require "DataStruct"
local Time       = require "Lib.Time"
local Basis      = require "Basis"
local Updater    = require "Updater"

-- USER INPUT.
local pOrder       = 0                 -- Polynomial order.
local basisName    = "Ser"             -- Type of polynomial basis.
local numCells     = {32, 32, 32}      -- Number of cells.
local lower        = {-1./2., -1./2., -1./2.}  -- Lower domain boundaries.
local upper        = { 1./2.,  1./2.,  1./2.}  -- Upper domain boundaries.

local periodicDirs = {1, 2, 3}         -- Periodic directions.

local diffCoeff    = {0.01, 0.01, 0.01}-- Diffusion coefficient.

local tStart       = 0.0               -- Start time.
local tEnd         = 1.5               -- End time.
local nFrames      = 2                 -- Number of frames to write out.
local cflNum       = 1./(2*pOrder+1)   -- CFL factor.

local initDt       = tEnd-tStart       -- initial time-step

-- Initial condition for simulation.
local function fInitial(xn, lower, upper)
   local x    = xn[1]
   local y    = xn[2]
   local z    = xn[3]

   -- A single sine wave
   local Pi   = math.pi
   local Lx      = {upper[1] - lower[1], upper[2] - lower[2], upper[3] - lower[3]}
   local kNum    = {2*Pi/Lx[1], 2*Pi/Lx[2], 2*Pi/Lx[3]}
   return math.sin(kNum[1]*x)*math.cos(kNum[2]*y)*math.sin(kNum[3]*z)
end

-- Analytic solution for verification.
local function fAnalytic(t, xn, lower, upper, Dcoeff)
   local x    = xn[1]
   local y    = xn[2]
   local z    = xn[3]

   -- Analytic answer to diffusing a single sine wave.
   local Pi   = math.pi
   local Lx      = {upper[1] - lower[1], upper[2] - lower[2], upper[3] - lower[3]}
   local kNum    = {2*Pi/Lx[1], 2*Pi/Lx[2], 2*Pi/Lx[3]}
   -- Single sine mode.
   return math.sin(kNum[1]*x)*math.cos(kNum[2]*y)*math.sin(kNum[3]*z)
         *math.exp(-(Dcoeff[1]*kNum[1]^2+Dcoeff[2]*kNum[2]^2+Dcoeff[3]*kNum[3]^2)*t)
end

-- SETUP.
local frame        = 1
local step         = 1
local tCurr        = tStart
local frameT       = tEnd/nFrames
local frameCount   = 1

local function createGrid(lo, up, nCells, pDirs)
   local gridOut = Grid.RectCart {
      lower        = lo,
      upper        = up,
      cells        = nCells,
      periodicDirs = pDirs,
   }
   return gridOut
end

local function createBasis(dim, pOrder, bKind)
   local basis
   if (bKind=="Ser") then
      basis = Basis.CartModalSerendipity { ndim = dim, polyOrder = pOrder }
   elseif (bKind=="Max") then
      basis = Basis.CartModalMaxOrder { ndim = dim, polyOrder = pOrder }
   elseif (bKind=="Tensor") then
      basis = Basis.CartModalTensor { ndim = dim, polyOrder = pOrder }
   else
      assert(false,"Invalid basis")
   end
   return basis
end

local function createField(grid, basis, numComp, vComp)
   numComp = numComp or basis:numBasis()
   vComp   = vComp or 1
   local fld = DataStruct.Field {
      onGrid        = grid,
      numComponents = numComp*vComp,
      ghost         = {1, 1},
      metaData = {
         polyOrder = basis:polyOrder(),
         basisType = basis:id()
      }
   }
   fld:clear(0.0)
   return fld
end

local grid  = createGrid(lower, upper, numCells, periodicDirs)
local basis = createBasis(grid:ndim(), pOrder, basisName)

-- States bfore and after a time-integration step.
local q    = createField(grid, basis)
local qNew = createField(grid, basis)
-- Analytic solution.
local qA   = createField(grid, basis)

local constDiffusionSlvr = Updater.ConstDiffusionSimple {
   onGrid      = grid,
   basis       = basis,
   cfl         = cflNum,
   coefficient = diffCoeff,
   components  = {1},
}

local project = Updater.ProjectOnBasis {
   onGrid   = grid,
   basis    = basis,
   evaluate = function (t, xn) return 1.0 end -- Set a dummy function for now.
}

-- RUNNING
-- Projecting initial condition.
project:setFunc(function(t,xn) return fInitial(xn, lower, upper) end)
project:advance(0.0, {}, {q})
q:write("qSimulation_0.bp", tStart, 0)  -- Write out initial condition.
q:sync()

-- Projecting. the initial analytic solution.
project:setFunc(function(t,xn) return fAnalytic(t, xn, lower, upper, diffCoeff) end)
project:advance(tCurr, {}, {qA})
qA:write("qAnalytic_0.bp", tCurr)

-- Function to take a single forward-euler time-step.
local function forwardEuler(tCurr, dt, fIn, fOut)
   constDiffusionSlvr:setDtAndCflRate(dt)
   local status, dtSuggested = constDiffusionSlvr:advance(tCurr, {fIn}, {fOut})
   fOut:sync()
   return status, dtSuggested
end

-- Main loop.
local tmSimStart = Time.clock()
local myDt = initDt
while true do
   if tCurr+myDt > tEnd then myDt = tEnd-tCurr end

   local status, dtSuggested = forwardEuler(tCurr, myDt, q, qNew)

   if status then
      tCurr = tCurr + myDt
      myDt = dtSuggested
      step = step + 1
      q:copy(qNew)

      if tCurr >= frameCount*frameT or math.abs(tCurr-frameCount*frameT) < 1e-10 then
         q:write(string.format("qSimulation_%d.bp",frameCount), tCurr, frameCount)

         project:advance(tCurr, {}, {qA})
         qA:write(string.format("qAnalytic_%d.bp",frameCount), tCurr, frameCount)

         frameCount = frameCount+1
      end
      if (tCurr >= tEnd) then
         break
      end
   else
      print(string.format(" ** Time step %g too large! Will retake with dt %g\n", myDt, dtSuggested))
      myDt = dtSuggested
   end
end
local tmSimEnd = Time.clock()
print(string.format(" Number of steps taken: %i", step))
print(string.format(" Total time for simulation %g\n", tmSimEnd - tmSimStart))
