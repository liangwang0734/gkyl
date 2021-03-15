-- Gkyl ------------------------------------------------------------------------
--
-- PlasmaOnCartGrid support code:  Braginskii friction and thermal force and
-- heating, i.e., the R and Q terms.
--
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local SourceBase = require "App.Sources.SourceBase"
local Proto = require "Lib.Proto"
local Updater = require "Updater"

local BraginskiiForceAndHeatingSource = Proto(SourceBase)

function BraginskiiForceAndHeatingSource:init(tbl)
   self.tbl = tbl
end

function BraginskiiForceAndHeatingSource:fullInit(appTbl)
   local tbl = self.tbl -- previously stored table

   self.cfl = nil
   self.grid = nil
   self.slvr = nil
end

function BraginskiiForceAndHeatingSource:setName(nm)
   self.name = nm
end

function BraginskiiForceAndHeatingSource:setCfl(cfl)
   self.cfl = cfl
end
function BraginskiiForceAndHeatingSource:setConfGrid(grid)
   self.grid = grid
end

function BraginskiiForceAndHeatingSource:createSolver(species, field)
   local tbl = self.tbl

   local mass, charge = {}, {}
   for i, nm in ipairs(tbl.species) do
      mass[i] = species[nm]:getMass()
      charge[i] = species[nm]:getCharge()
   end

   self.slvr = Updater.BraginskiiForceAndHeating {
      onGrid           = self.grid,
      mass             = mass,
      charge           = charge,
      gasGamma         = tbl.gasGamma,
      epsilon0         = tbl.epsilon0,
      tauElectron      = tbl.tauElectron,
      calcTauElectron  = tbl.calcTauElectron,
      coulombLogarithm = tbl.coulombLogarithm,
   }
end

function BraginskiiForceAndHeatingSource:updateSource(
      tCurr, dt, speciesVar, fieldVar, speciesVarBuf, fieldVarBuf)
   local tbl = self.tbl

   local outVars = {}
   for i, nm in ipairs(tbl.species) do
      outVars[i] = speciesVar[nm]
   end
   outVars[#tbl.species+1] = fieldVar

   local inVars = {}
   for i, nm in ipairs(tbl.species) do
      inVars[i] = speciesVarBuf[nm]
   end
   inVars[#tbl.species+1] = fieldVarBuf

   self.slvr:setDtAndCflRate(dt, nil)
   return self.slvr:advance(tCurr, inVars, outVars)
end

function BraginskiiForceAndHeatingSource:write(tm, frame)
end

function BraginskiiForceAndHeatingSource:totalTime()
   return self.slvr.totalTime
end

return BraginskiiForceAndHeatingSource
