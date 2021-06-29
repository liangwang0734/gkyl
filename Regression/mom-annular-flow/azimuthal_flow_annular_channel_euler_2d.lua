-- Neutral fluid (5-moment) modeling of an azimuthal flow in the cross section
-- of an annular channel.

local Moments = require("App.PlasmaOnCartGrid").Moments()
local Euler = require "Eq.Euler"
local BoundaryCondition = require "Updater.BoundaryCondition"
local Logger = require "Lib.Logger"

--------------------------------------------------------------------------------

local gasGamma = 5/3
local rho0 = 1.
local p0 = 1.
local vtheta__cs = 0.1
local r_inn = 0.5
local r_out = 1.5
local Nx, Ny = 100, 100
local decompCuts = {1, 1}

local cs0 = math.sqrt(gasGamma*p0/rho0)
local vtheta0 = cs0 * vtheta__cs

local r_max = r_out * 1.1
local xlo, xup = -r_max, r_max
local ylo, yup = -r_max, r_max

local tEnd = 1
local nFrame = 1

local init = function (t, xn)
   local x, y = xn[1], xn[2]
   local r = math.sqrt(x*x+y*y)
   local rho = rho0
   local vx, vy, vz = 0, 0, 0
   if r < r_out and r > r_inn then
      local a = (r - r_inn) / (r_out-r_inn)
      -- Make the azimuthal flow velocity to vanish at inner and outer walls.
      local vtheta = vtheta0 * 0.5 * (math.cos(2 * math.pi * (a + 0.5)) + 1)
      vx = -vtheta * y / r
      vy = vtheta * x / r
   end
   local p = p0
   local er = p / (gasGamma - 1.) + 0.5 * rho * (vx^2 + vy^2 + vz^2)
   return rho, rho*vx, rho*vy, rho*vz, er
end

-- Define the annular channel between two radial walls.
local embedded_boundary_mask_func = function (t, xn)
   local x, y = xn[1], xn[2]
   local r = math.sqrt(x*x + y*y)
   if (r>r_inn) and (r<r_out) then return 1 else return -1 end
end

local embedded_bc = {
   function(dir, tm, idxSkin, qSkin, qGhost, xcGhost, xcSkin)
      local rhoS = qSkin[1]
      local vxS = qSkin[2] / rhoS
      local vyS = qSkin[3] / rhoS
      local vzS = qSkin[4] / rhoS
      local hS = qSkin[5] - 0.5 * rhoS * (vxS*vxS + vyS*vyS + vzS*vzS)

      -- "No-slip" flow boundary condition.
      -- TODO: Try to set vthetaGhost = - vthetaSkin.
      local rhoG = rhoS
      local vxG = 0
      local vyG = 0
      local vzG = 0
      local hG = hS
      
      qGhost[1] = rhoG
      qGhost[2] = rhoG * vxG
      qGhost[3] = rhoG * vyG
      qGhost[4] = rhoG * vzG
      qGhost[5] = hG + 0.5 * rhoG * (vxG*vxG + vyG*vyG + vzG*vzG)
   end
}

--------------------------------------------------------------------------------

local logger = Logger { logToFile = True }
local log = function(...) logger(string.format(...).."\n") end

log("%30s = %g", "r_inn", r_inn)
log("%30s = %g", "r_out", r_out)
log("%30s = %g", "vtheta0/cs0", vtheta0/cs0)
log("%30s = %g", "Nx", Nx)
log("%30s = %g", "Ny", Ny)
log("%30s = %g", "tEnd", tEnd)
log("%30s = %g", "nFrame", nFrame)

--------------------------------------------------------------------------------
app = Moments.App {
   logToFile = true,

   tEnd = tEnd,
   nFrame = nFrame,
   lower = {xlo, ylo},
   upper = {xup, yup},
   cells = {Nx, Ny},
   cflFrac = cfl,
   timeStepper = "fvDimSplit",
   decompCuts = decompCuts,

   fluid = Moments.Species {
      charge = 0.0, mass = 1.0,
      equation = Euler { gasGamma = gasGamma },
      equationInv = Euler { gasGamma = gasGamma, numericalFlux="lax" },
      forceInv = false,
      init = init,
      bcx = { Moments.Species.bcCopy, Moments.Species.bcCopy },
      bcy = { Moments.Species.bcCopy, Moments.Species.bcCopy },
      hasSsBnd = true,
      inOutFunc = embedded_boundary_mask_func,
      ssBc = { embedded_bc },
   },   
}

app:run()
