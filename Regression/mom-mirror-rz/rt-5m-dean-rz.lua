-- Gkyl ------------------------------------------------------------------------
-- 5-moment modeling of dean flow in the 2d axisymmetric (r, z) plane
-- r -> x, phi -> y, z -> z; One cell along phi (y)

local Moments = require("App.PlasmaOnCartGrid").Moments()
local Euler = require "Eq.Euler"
local BoundaryCondition = require "Updater.BoundaryCondition"
local Logger = require "Lib.Logger"

local logger = Logger {
   logToFile = True
}

local log = function(...)
   logger(string.format(...))
   logger("\n")
end

-- basic normalizations
gasGamma = 1.4
vA0 = 1.
mu0 = 1.
rho0 = 1.
l0 = 1.

n0 = 1.
mi = 1.

-- problem-specific mhd parameters
T0_mhd = 0.002
Bz0_B0 = 1

-- non-mhd parameters
lightSpeed__vA0 = 10
mi__me = 25.
Ti0__Te0 = 1.
di0__l0 = 0.2

-- derived mhd parameters
B0 = vA0 * math.sqrt(mu0 * rho0)
Bz0 = B0
pmag0 = B0*B0 / 2 / mu0
p0 = n0 * T0_mhd
tA0 = l0 / vA0

-- derived non-mhd parameters
lightSpeed = vA0 * lightSpeed__vA0
epsilon0 = 1. / lightSpeed^2 / mu0
me = mi / mi__me
rhoe0 = rho0 / (1 + mi__me)
pe0 = p0 / (1 + Ti0__Te0)
rhoi0 = rho0 - rhoe0
pi0 = p0 - pe0

di0 = l0 * di0__l0
wpi0 = lightSpeed__vA0 / di0
qi = wpi0 * math.sqrt(epsilon0 * mi / n0)
qe = -qi

-- domain size
r_inn = 0.45 * l0
r_out = 1.45 * l0
Lz = 5 * l0
Nr, Nz = 60, 300

de0 = di0 / math.sqrt(mi/me)
dz = Lz / Nz
dr = (r_out - r_inn) / Nr

log("%30s = %g", "Lz/de0", Lz/de0)
log("%30s = %g", "(r_out-r_inn)/de0", (r_out-r_inn)/de0)
log("%30s = %g = 1/%g", "dz/de0", dz/de0, de0/dz)
log("%30s = %g = 1/%g", "dr/de0", dr/de0, de0/dr)

momentApp = Moments.App {
   logToFile = true,

   tEnd = 10. * tA0,
   nFrame = 1,
   lower = {r_inn, -1, 0},
   upper = {r_out, 1, Lz},
   cells = {Nr, 1, Nz},
   timeStepper = "fvDimSplit",

   periodicDirs = {2, 3},  -- periodic along phi and z
   decompCuts = {1, 1, 4},

   elc = Moments.Species {
      charge = qe, mass = me,
      equation = Euler { gasGamma = gasGamma },
      equationInv = Euler { gasGamma = gasGamma, numericalFlux = "lax" },
      forceInv = false,
      init = function (t, xn)
         local r, phi, z = xn[1], xn[2], xn[3]
         local rho = rhoe0
         local vr = 0.
         local vphi = 0.
         local vz = 0.
         local p = pe0
         local er = p / (gasGamma - 1.) + 0.5 * rho * (vr^2 + vphi^2 + vz^2)
         return rho, rho*vr, rho*vphi, rho*vz, er
      end,
      bcx = { Euler.bcWall, Euler.bcWall },
   },

   ion = Moments.Species {
      charge = qi, mass = mi,
      equation = Euler { gasGamma = gasGamma },
      equationInv = Euler { gasGamma = gasGamma, numericalFlux = "lax" },
      forceInv = false,
      init = function (t, xn)
         local r, phi, z = xn[1], xn[2], xn[3]
         local rho = rhoi0
         local vr = 0.
         local vphi = 0.
         local vz = 0.
         local p = pi0
         local er = p / (gasGamma - 1.) + 0.5 * rho * (vr^2 + vphi^2 + vz^2)
         return rho, rho*vr, rho*vphi, rho*vz, er
      end,
      bcx = { Euler.bcWall, Euler.bcWall },  -- radial
   },

   field = Moments.Field {
      epsilon0 = epsilon0, mu0 = epsilon0,
      init = function (t, xn)
         local r, phi, z = xn[1], xn[2], xn[3]
         local Er = 0
         local Ephi = 0
         local Ez = 0
         local Br = 0
         local Bphi = 0
         local Bz = Bz0
         return Er, Ephi, Ez, Br, Bphi, Bz
      end,
      bcx = { Moments.Field.bcReflect, Moments.Field.bcReflect },  -- radial
   },

   axisymSource = Moments.CollisionlessEmSource {
      species = {"elc", "ion"},
      timeStepper = "axisymmetric",
      epsilon0 = epsilon0, mu0 = epsilon0,
      gasGamma = gasGamma,
   },   

   emSource = Moments.CollisionlessEmSource {
      species = {"elc", "ion"},
      timeStepper = "direct",
   },   
}
-- run application
momentApp:run()
