local Grid       = require "Grid"
local DataStruct = require "DataStruct"
local Updater    = require "Updater"
local Basis      = require "Basis"

local gasGamma    = 5./3.
local epsilon0    = 1.
local elcMass     = 1.
local ionMass     = 100.
local ionCharge   = 1.
local elcCharge   = -ionCharge

local calcTauElectron = false
local tauElectron = 1.

if calcTauElectron then
   tauElectron = nil
end

local lower        = { 1./2.,  1./2.,  1./2.}
local upper        = {-1./2., -1./2., -1./2.}
local cells        = {32, 32, 32}
local cells        = {4, 4, 4}

local periodicDirs = {1, 2, 3}

local function initElc (t, xn)
   local x, y, z = xn[1], xn[2], xn[3]
   return 1, 0, 0, 0, 1
end

local function initIon (t, xn)
   local x, y, z = xn[1], xn[2], xn[3]
   return 1, 0, 0, 0, 1
end

local function initEmf (t, xn)
   local x, y, z = xn[1], xn[2], xn[3]
   return 0, 0, 0, 1, 1, 1, 0, 0
end

local grid = Grid.RectCart {
   lower        = lower,
   upper        = upper,
   cells        = cells,
   periodicDirs = periodicDirs,
}
local basis = Basis.CartModalSerendipity { ndim=grid:ndim(), polyOrder=0 }

local function createField(grid, basis, numComponents)
   local fld = DataStruct.Field {
      onGrid        = grid,
      numComponents = numComponents,
      ghost         = {2, 2},
      metaData = {
         polyOrder = basis:polyOrder(),
         basisType = basis:id()
      }
   }
   fld:clear(0.0)
   return fld
end

local braginskiiHeatConduction = Updater.BraginskiiHeatConduction {
   onGrid      = grid,
   nFluids     = 2,
   tau         = {1,1},
   gasGamma = gasGamma,
   mass = {elcMass, ionMass},
   charge = {elcCharge, ionCharge},
   epsilon0 = epsilon0,
}

local project = Updater.ProjectOnBasis {
   onGrid   = grid,
   basis    = basis,
   evaluate = function (t, xn) return 1.0 end -- Set a dummy function for now.
}

local elc = createField(grid, basis, 5)
local ion = createField(grid, basis, 5)
local emf = createField(grid, basis, 8)
local elcBuf = createField(grid, basis, 5)
local ionBuf = createField(grid, basis, 5)
local emfBuf = createField(grid, basis, 8)

local tStart = 0
project:setFunc(initElc)
project:advance(0.0, {}, {elc})
elc:write("elc_0.bp", tStart, 0)
elc:sync()

project:setFunc(initIon)
project:advance(0.0, {}, {ion})
ion:write("ion_0.bp", tStart, 0)
ion:sync()

project:setFunc(initEmf)
project:advance(0.0, {}, {emf})
emf:write("emf_0.bp", tStart, 0)
emf:sync()

local tCurr = tStart
local dt = 1
braginskiiHeatConduction:setDtAndCflRate(dt)
local status, dtSuggested = braginskiiHeatConduction:advance(
   tCurr, {elcBuf,ionBuf,emfBuf}, {elc,ion,emf})
elc:sync()
ion:sync()
emf:sync()

elc:write("elc_1.bp", tStart+dt, 1)
ion:write("ion_1.bp", tStart+dt, 1)
emf:write("emf_1.bp", tStart+dt, 1)
