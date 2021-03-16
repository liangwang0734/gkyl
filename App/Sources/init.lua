-- Gkyl ------------------------------------------------------------------------
--
-- For accessing source objects
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local SourceBase = require "App.Sources.SourceBase"
local CollisionlessEmSource = require "App.Sources.CollisionlessEmSource"
local TenMomentRelaxSource = require "App.Sources.TenMomentRelaxSource"
local BraginskiiForceAndHeatingSource = require "App.Sources.BraginskiiForceAndHeatingSource"
local BraginskiiHeatConductionSource = require "App.Sources.BraginskiiHeatConductionSource"
local BraginskiiViscosityDiffusionSource = require "App.Sources.BraginskiiViscosityDiffusionSource"

return {
   SourceBase = SourceBase,
   CollisionlessEmSource = CollisionlessEmSource,
   TenMomentRelaxSource = TenMomentRelaxSource,
   BraginskiiForceAndHeatingSource = BraginskiiForceAndHeatingSource,
   BraginskiiHeatConductionSource = BraginskiiHeatConductionSource,
   BraginskiiViscosityDiffusionSource = BraginskiiViscosityDiffusionSource,
}
