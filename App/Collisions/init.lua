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
local GkIonization      = require "App.Collisions.GkIonization"
local VmBGKCollisions   = require "App.Collisions.VmBGKCollisions"
local VmLBOCollisions   = require "App.Collisions.VmLBOCollisions"
local VmIonization      = require "App.Collisions.VmIonization"
local VoronovIonization = require "App.Collisions.VoronovIonization"

return {
  CollisionsBase    = CollisionsBase,
  FluidDiffusion    = FluidDiffusion,
  GkBGKCollisions   = GkBGKCollisions,
  GkLBOCollisions   = GkLBOCollisions,
  GkIonization      = GkIonization,
  VmBGKCollisions   = VmBGKCollisions,
  VmLBOCollisions   = VmLBOCollisions,
  VmIonization      = VmIonization,
  VoronovIonization = VoronovIonization,
}
