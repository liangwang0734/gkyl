-- Gkyl ------------------------------------------------------------------------
--
-- PlasmaOnCartGrid support code: Viscosity as velocity diffusion and optional
-- viscous heating.
--
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local SourceBase = require "App.Sources.SourceBase"
local Proto = require "Lib.Proto"
local Updater = require "Updater"

local BraginskiiViscosityDiffusionSource = Proto(SourceBase)

function BraginskiiViscosityDiffusionSource:init(tbl)
   self.tbl = tbl
end

function BraginskiiViscosityDiffusionSource:fullInit(appTbl)
   local tbl = self.tbl -- previously stored table

   self.cfl = nil
   self.grid = nil
   self.slvr = nil
end

function BraginskiiViscosityDiffusionSource:setName(nm)
   self.name = nm
end

function BraginskiiViscosityDiffusionSource:setCfl(cfl)
   self.cfl = cfl
end
function BraginskiiViscosityDiffusionSource:setConfGrid(grid)
   self.grid = grid
end

function BraginskiiViscosityDiffusionSource:createSolver(species, field)
   local tbl = self.tbl

   local mass, charge = {}, {}
   for i, nm in ipairs(tbl.species) do
      mass[i] = species[nm]:getMass()
      charge[i] = species[nm]:getCharge()
   end

   self.slvr = Updater.BraginskiiViscosityDiffusion {
      onGrid           = self.grid,
      numFluids        = #tbl.species,
      mass             = mass,
      charge           = charge,
      gasGamma         = tbl.gasGamma,
      epsilon0         = tbl.epsilon0,
      tau              = tbl.tau,
      calcTau          = tbl.calcTau,
      coulombLogarithm = tbl.coulombLogarithm,
      hasHeating       = tbl.hasHeating,
      coordinate       = tbl.coordinate,
   }
end

function BraginskiiViscosityDiffusionSource:updateSource(
      tCurr, dt, speciesVar, fieldVar, speciesBuf, fieldBuf, allSpecies)
   local tbl = self.tbl

   local outVars = {}
   for i, nm in ipairs(tbl.species) do
      outVars[i] = speciesVar[nm]

      for d=1,self.grid:ndim() do
         allSpecies[nm]:applyBc(tCurr, speciesVar[nm], d)
      end
   end
   outVars[#tbl.species+1] = fieldVar

   local inVars = {}
   for i, nm in ipairs(tbl.species) do
      inVars[i] = speciesBuf[nm]
   end
   inVars[#tbl.species+1] = fieldBuf

   self.slvr:setDtAndCflRate(dt, nil)
   -- return true, GKYL_MAX_DOUBLE
   return self.slvr:advance(tCurr, inVars, outVars)
end

function BraginskiiViscosityDiffusionSource:write(tm, frame)
end

function BraginskiiViscosityDiffusionSource:totalTime()
   return self.slvr.totalTime
end

return BraginskiiViscosityDiffusionSource
