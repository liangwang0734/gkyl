-- Two-fluid-five-moment-Maxwell modeling of an azimuthal flow in the cross   --
-- section of an annular channel.                                             --

-- PREAMBLE - ------------------------------------------------------------------

local Moments = require("App.PlasmaOnCartGrid").Moments()
local Euler = require "Eq.Euler"
local BoundaryCondition = require "Updater.BoundaryCondition"
local Logger = require "Lib.Logger"

-- PARAMETERS ------------------------------------------------------------------

local gasGamma = 5/3
local lightSpeed = 1.
local mu0 = 1.
local me = 1.
local qe = -1.
local n0 = 1.
local wpe0_wce0 = 10.

local epsilon0 = 1. / lightSpeed^2 / mu0
local wpe0 = math.sqrt(n0 * qe * qe / me / epsilon0)
local de0 = lightSpeed / wpe0
local wce0 = wpe0 / wpe0_wce0
local B0 = math.abs(wce0 * me / qe)

local Bz_inn = -B0
local Bz_out = B0
local betae_inn = 10

local pmag_inn = Bz_inn*Bz_inn / 2 / mu0
local pe_inn = betae_inn * pmag_inn
 
local r_inn = 0
local r_out = 50 * de0
local r_max = r_out * 1.1
local Nx, Ny = 201, 201
local decompCuts = {2, 2}
local cfl = 0.95
local tEnd = 10 * 2 * r_max / lightSpeed
local nFrame = 10

local xlo, xup = -r_max, r_max
local ylo, yup = -r_max, r_max
local dx = (xup-xlo) / Nx
local suggestedDt = dx / lightSpeed * cfl

local q = qe
local m = me
local n = n0
local eq_r0 = 10 * de0
local eq_hw = 1 * de0
local eq_p0 = pe_inn
local eq_nr = Nx * 5
local eq_dr = (r_out - r_inn) / eq_nr
local eq_pressure = {eq_p0}

local logger = Logger { logToFile = True }
local log = function(...) logger(string.format(...).."\n") end

log("%30s = %g", "lightSpeed", lightSpeed)
log("%30s = %g", "mu0", mu0)
log("%30s = %g", "wpe0", wpe0)
log("%30s = %g", "de0", de0)
log("%30s = %g", "r_inn", r_inn)
log("%30s = %g", "r_out", r_out)
log("%30s = %g", "r_max", r_max)
log("%30s = %g", "Nx", Nx)
log("%30s = %g", "Ny", Ny)
log("%30s = %g", "dx", dx)
log("%30s = %g", "r_inn/de0", r_inn/de0)
log("%30s = %g", "r_out/de0", r_out/de0)
log("%30s = %g", "(r_out-r_inn)/de0", (r_out-r_inn)/de0)
log("%30s = %g = 1/%g", "dx/de0", dx/de0, de0/dx)
log("%30s = %g", "tEnd", tEnd)
log("%30s = %g", "tEnd*wpe0", tEnd*wpe0)
log("%30s = %g", "tEnd*wce0", tEnd*wce0)
log("%30s = %g", "nFrame", nFrame)

-- INITIAL CONDITIONS ----------------------------------------------------------

-- The equilibrium following Hakim2007pop assuming stationary ions, and no theta
-- or z dependence. Magnetic field Bz reversal and electron current Jtheta are
-- balanced in the Ampere's law. The pressure gradient balances the centrifugal
-- and Lorentz forces along the radial direction.

local calc_eq_Bz = function (r)
   return Bz_inn + 0.5 * (Bz_out-Bz_inn) * (math.tanh((r-eq_r0)/eq_hw) + 1)
end

local calc_eq_vtheta = function (r)
   tanh_value = math.tanh((r-eq_r0)/eq_hw)
   return (0.5 * (Bz_out-Bz_inn) * ( 1 - tanh_value*tanh_value ) / eq_hw) * (-1./mu0/q)
end

-- Precompute the pressure at different radii by integrating dp/dr.
for i=2,eq_nr do
   local r = r_inn + eq_dr * i
   local Bz = calc_eq_Bz(r)
   local ut = calc_eq_vtheta(r)
   local dpdr = n*q*ut*Bz - m*n*ut*ut/r
   eq_pressure[i] = eq_pressure[i-1] + dpdr * eq_dr
end
print(eq_pressure[1], eq_pressure[eq_nr])

local calc_eq_p = function (r)
   local idx = math.floor((r-r_inn) / eq_dr)
   idx = math.max(idx, 1)
   idx = math.min(idx, eq_nr)
   return eq_pressure[idx]
end

log("-- Equilibrium parameters --")
log("%30s = %g", "eq mu0", mu0)
log("%30s = %g", "eq q", q)
log("%30s = %g", "eq n", n)
log("%30s = %g", "eq m", m)
log("%30s = %g", "eq r0", eq_r0)
log("%30s = %g", "eq hw", eq_hw)
log("%30s = %g", "eq r_inn", r_inn)
log("%30s = %g", "eq r_out", r_out)
log("%30s = %g", "eq nr", eq_nr)
log("%30s = %g", "eq dr", eq_dr)
log("%30s = %g", "eq vtheta(r0)", calc_eq_vtheta(eq_r0))
log("%30s = %g", "eq vtheta(r0)/lightSpeed", calc_eq_vtheta(eq_r0)/lightSpeed)

-- CREATE AND RUNNING THE APP --------------------------------------------------

local app = Moments.App {
   logToFile = true,

   tEnd = tEnd,
   nFrame = nFrame,
   lower = {xlo, ylo},
   upper = {xup, yup},
   cells = {Nx, Ny},
   cflFrac = cfl,
   suggestedDt = suggestedDt,
   timeStepper = "fvDimSplit",
   decompCuts = decompCuts,

   elc = Moments.Species {
      charge = qe, mass = me,
      equation = Euler { gasGamma = gasGamma },
      equationInv = Euler { gasGamma = gasGamma, numericalFlux="lax" },
      forceInv = false,
      init = function (t, xn)
         local x, y = xn[1], xn[2]
         local r = math.sqrt(x*x+y*y)
         local rho = n0 * me
         local vtheta = calc_eq_vtheta(r)
         local vx = -vtheta * y / r
         local vy = vtheta * x / r
         if r<eq_dr then
            vx, vy = 0, 0
         end
         local vz = 0.
         local p = calc_eq_p(r)
         local er = p / (gasGamma - 1.) + 0.5 * rho * (vx^2 + vy^2 + vz^2)
         return rho, rho*vx, rho*vy, rho*vz, er
      end,
      bcx = { Moments.Species.bcCopy, Moments.Species.bcCopy },
      bcy = { Moments.Species.bcCopy, Moments.Species.bcCopy },
   },

   field = Moments.Field {
      epsilon0 = epsilon0, mu0 = mu0,
      init = function (t, xn)
         local x, y = xn[1], xn[2]
         local r = math.sqrt(x*x+y*y)
         local Bx, By, Bz = 0, 0, calc_eq_Bz(r)
         local Ex, Ey, Ez = 0, 0, 0
         local phiE = 0
         local phiB = 0
         return Ex, Ey, Ez, Bx, By, Bz, phiE, phiB
      end,
      bcx = { Moments.Field.bcCopy, Moments.Field.bcCopy },
      bcy = { Moments.Field.bcCopy, Moments.Field.bcCopy },
   },

   emSource = Moments.CollisionlessEmSource {
      species = {"elc"},
      timeStepper = "time-centered",
   },
}

app:run()
