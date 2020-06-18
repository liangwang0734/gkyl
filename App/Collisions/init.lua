-- Gkyl ------------------------------------------------------------------------
--
-- For accessing collisions objects
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local CollisionsBase    = require "App.Collisions.CollisionsBase"
local FluidDiffusion    = require "App.Collisions.FluidDiffusion"
local GkBGKCollisions   = require "App.Collisions.GkBGKCollisions"
local GkLBOCollisions   = require "App.Collisions.GkLBOCollisions"
local GkChargeExchange  = require "App.Collisions.GkChargeExchange"
local GkIonization      = require "App.Collisions.GkIonization"
local VmBGKCollisions   = require "App.Collisions.VmBGKCollisions"
local VmLBOCollisions   = require "App.Collisions.VmLBOCollisions"
local VmChargeExchange  = require "App.Collisions.VmChargeExchange"
local VmIonization      = require "App.Collisions.VmIonization"
local VoronovIonization = require "App.Collisions.VoronovIonization"

return {
  CollisionsBase    = CollisionsBase,
  FluidDiffusion    = FluidDiffusion,
  GkBGKCollisions   = GkBGKCollisions,
  GkLBOCollisions   = GkLBOCollisions,
  GkChargeExchange  = GkChargeExchange,
  GkIonization      = GkIonization,
  VmBGKCollisions   = VmBGKCollisions,
  VmLBOCollisions   = VmLBOCollisions,
  VmChargeExchange  = VmChargeExchange,
  VmIonization      = VmIonization,
  VoronovIonization = VoronovIonization,
}
