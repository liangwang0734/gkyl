-------------------------------------------------------------------------------
-- Two-fluid-five-moment-Maxwell modeling of an azimuthal flow in the cross  --
-- section of an annular channel.                                            --
-------------------------------------------------------------------------------

-- PREAMBLE - ------------------------------------------------------------------

local Moments = require("App.PlasmaOnCartGrid").Moments()
local Euler = require "Eq.Euler"
local BoundaryCondition = require "Updater.BoundaryCondition"
local Logger = require "Lib.Logger"

-- PARAMETERS ------------------------------------------------------------------

-- Basic normalizations.
local gasGamma = 5/3
local vA0 = 1.
local mu0 = 1.
local rho0 = 1.
local l0 = 1.
local n0 = 1.
local mi = 1.

-- Problem-specific mhd parameters.
local Bz0_B0 = 1
local T0_mhd = 0.002
local MA0 = 0.5  -- v_theta0/vA0; set to nonzero to set a Er x Bz drift.

-- Non-mhd parameters.
local lightSpeed__vA0 = 10
local mi__me = 25.
local Ti0__Te0 = 1.
local di0__l0 = 0.1

-- Derived mhd parameters.
local B0 = vA0 * math.sqrt(mu0 * rho0)
local Bz0 = B0
local pmag0 = B0*B0 / 2 / mu0
local p0 = n0 * T0_mhd
local tA0 = l0 / vA0
local v_theta0 = vA0 * MA0
 
-- Derived non-mhd parameters
local lightSpeed = vA0 * lightSpeed__vA0
local epsilon0 = 1. / lightSpeed^2 / mu0
local me = mi / mi__me
local rhoe0 = rho0 / (1 + mi__me)
local pe0 = p0 / (1 + Ti0__Te0)
local rhoi0 = rho0 - rhoe0
local pi0 = p0 - pe0
local p0 = pe0 + pi0
local rho0 = rhoe0 + rhoi0
local vS0 = math.sqrt(gasGamma*p0/rho0)
local Er0 = -v_theta0 * Bz0

local di0 = l0 * di0__l0
local wpi0 = lightSpeed__vA0 / di0
local qi = wpi0 * math.sqrt(epsilon0 * mi / n0)
local qe = -qi
local de0 = di0 / math.sqrt(mi/me)

-- Domain, spatial and time resolutions.
local r_inn = 0.45 * l0
local r_out = 1.45 * l0
local r_max = r_out * 1.1
local Nx, Ny = 100, 100
local decompCuts = {1, 1}
local cfl = 0.5
local tEnd = tA0*50
local nFrame = 50

local xlo, xup = -r_max, r_max
local ylo, yup = -r_max, r_max
local dx = (xup-xlo) / Nx

local logger = Logger { logToFile = True }
local log = function(...) logger(string.format(...).."\n") end

log("%30s = %g", "r_inn", r_inn)
log("%30s = %g", "r_out", r_out)
log("%30s = %g", "Nx", Nx)
log("%30s = %g", "Ny", Ny)
log("%30s = %g", "r_inn/di0", r_inn/di0)
log("%30s = %g", "r_out/di0", r_out/di0)
log("%30s = %g", "(r_out-r_inn)/di0", (r_out-r_inn)/di0)
log("%30s = %g = 1/%g", "dx/di0", dx/di0, di0/dx)
log("%30s = %g", "r_inn/de0", r_inn/de0)
log("%30s = %g", "r_out/de0", r_out/de0)
log("%30s = %g", "(r_out-r_inn)/de0", (r_out-r_inn)/de0)
log("%30s = %g = 1/%g", "dx/de0", dx/de0, de0/dx)
log("%30s = %g", "tEnd", tEnd)
log("%30s = %g", "tEnd/tA0", tEnd/tA0)
log("%30s = %g", "nFrame", nFrame)

log("%30s = %g", "v_theta0", v_theta0)
log("%30s = %g", "v_theta0/vA0", v_theta0/vA0)
log("%30s = %g", "v_theta0/vS0", v_theta0/vS0)
log("%30s = %g", "Bz0", Bz0)
log("%30s = %g", "Er0", Er0)

-- INITIAL CONDITIONS ----------------------------------------------------------

local init_v = function (x, y)
   local r = math.sqrt(x*x+y*y)
   local vx, vy, vz = 0, 0, 0
   if r < r_out and r > r_inn then
      -- Make the azimuthal flow velocity to vanish at inner and outer walls.
      local a = (r - r_inn) / (r_out-r_inn)
      local v_theta = v_theta0 * 0.5 * (math.cos(2 * math.pi * (a + 0.5)) + 1)
      vx = -v_theta * y / r
      vy = v_theta * x / r
   end
   return vx, vy, vz
end

local create_init_fluid = function (rho0, p0)
   return function (t, xn)
      local x, y = xn[1], xn[2]
      local rho = rho0
      local vx, vy, vz = init_v(x, y)
      local p = p0
      local er = p / (gasGamma - 1.) + 0.5 * rho * (vx^2 + vy^2 + vz^2)
      return rho, rho*vx, rho*vy, rho*vz, er
   end
end

local init_elc = create_init_fluid(rhoe0, pe0)
local init_ion = create_init_fluid(rhoe0, pe0)

local init_field = function (t, xn)
   local x, y = xn[1], xn[2]
   local vx, vy, vz = init_v(x, y)
   local Bx, By, Bz = 0, 0, Bz0
   local Ex = -vy * Bz + vz * By
   local Ey = -vz * Bx + vx * Bz
   local Ez = -vx * By + vy * Bx
   local phiE = 0
   local phiB = 0
   return Ex, Ey, Ez, Bx, By, Bz, phiE, phiB
end

-- EMBEDDED BOUNDARY CONDITIONS ------------------------------------------------

-- Define the annular channel between two radial walls.
local embedded_boundary_mask_func = function (t, xn)
   local x, y = xn[1], xn[2]
   local r = math.sqrt(x*x + y*y)
   if (r>r_inn) and (r<r_out) then return 1 else return -1 end
end

local embedded_bc_fluid = {
   function(dir, tm, idxSkin, qSkin, qGhost, xcGhost, xcSkin)
      local rhoS = qSkin[1]
      local vxS = qSkin[2] / rhoS
      local vyS = qSkin[3] / rhoS
      local vzS = qSkin[4] / rhoS
      local hS = qSkin[5] - 0.5 * rhoS * (vxS*vxS + vyS*vyS + vzS*vzS)

      -- "No-slip" flow boundary condition.
      -- TODO: Try to set v_thetaGhost = - v_thetaSkin.
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

local embedded_bc_elc = embedded_bc_fluid
local embedded_bc_ion = embedded_bc_fluid

local embedded_bc_field = Moments.Field.bcReflect

-- CREATE AND RUNNING THE APP --------------------------------------------------

local app = Moments.App {
   logToFile = true,

   tEnd = tEnd,
   nFrame = nFrame,
   lower = {xlo, ylo},
   upper = {xup, yup},
   cells = {Nx, Ny},
   cflFrac = cfl,
   timeStepper = "fvDimSplit",
   decompCuts = decompCuts,

   elc = Moments.Species {
      charge = qe, mass = me,
      equation = Euler { gasGamma = gasGamma },
      equationInv = Euler { gasGamma = gasGamma, numericalFlux="lax" },
      forceInv = false,
      init = init_elc,
      bcx = { {}, {} }, bcy = { {}, {} }, -- Dummy BCs for outer boundaries. 
      hasSsBnd = true,
      inOutFunc = embedded_boundary_mask_func,
      ssBc = { embedded_bc_elc },
   },

   ion = Moments.Species {
      charge = qi, mass = mi,
      equation = Euler { gasGamma = gasGamma },
      equationInv = Euler { gasGamma = gasGamma, numericalFlux="lax" },
      forceInv = false,
      init = init_ion,
      bcx = { {}, {} }, bcy = { {}, {} }, -- Dummy BCs for outer boundaries. 
      hasSsBnd = true,
      inOutFunc = embedded_boundary_mask_func,
      ssBc = { embedded_bc_ion },
   },

   field = Moments.Field {
      epsilon0 = epsilon0, mu0 = mu0,
      init = init_field,
      bcx = { {}, {} }, bcy = { {}, {} }, -- Dummy BCs for outer boundaries. 
      hasSsBnd = true,
      inOutFunc = embedded_boundary_mask_func,
      ssBc = { embedded_bc_field },
   },

   emSource = Moments.CollisionlessEmSource {
      species = {"elc", "ion"},
      timeStepper = "time-centered",
   },
}

app:run()
