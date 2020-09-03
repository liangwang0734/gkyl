-- Gkyl ------------------------------------------------------------------------
--
-- Species objects for use in multi-moment linear solver
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Proto = require "Lib.Proto"
local complex = require "sci.complex"
local Lin = require "Lib.Linalg"

-- clear out contents of matrix
local function matrixClear(m, val)
   for i = 1, m:numRows() do
      for j = 1, m:numCols() do
	 m[i][j] = val
      end
   end
end

-- Base species type
local Species = Proto()

function Species:init(tbl)
   self.mass, self.charge = tbl.mass, tbl.charge
end

-- Isothermal species type
local IsoThermal = Proto(Species)

function IsoThermal:init(tbl)
   IsoThermal.super.init(self, tbl)
   
   self.density = tbl.density -- number density
   
   assert(#tbl.velocity == 3, "Must provide all 3 components of velocity")
   self.velocity = tbl.velocity -- velocity vector
   
   self.temperature = tbl.temperature -- species temperature (0.0 is okay)
end

-- Number of equations in system
function IsoThermal:numEquations()
   return 4
end

-- Euler (5M) species type
local Euler = Proto(Species)

function Euler:init(tbl)
   Euler.super.init(self, tbl)

   self.density = tbl.density -- number density
   
   assert(#tbl.velocity == 3, "Must provide all 3 components of velocity")
   self.velocity = tbl.velocity -- velocity vector
   self.pressure = tbl.pressure  -- pressure
   
   self.gasGamma = tbl.gasGamma and tbl.gasGamma or 5.0/3.0
end

-- Number of equations in system
function Euler:numEquations()
   return 5
end

-- Construct matrix with contributions to moment-equation part of
-- dispersion matrix: this is a numEquations X numEquations sized
-- block
function Euler:calcMomDispMat(k, E, B)
   local D = Lin.ComplexMat(5, 5)
   matrixClear(D, 0.0)

   local gamma = self.gasGamma
   local n, p = self.density, self.pressure
   local u = self.velocity
   local m, qbym = self.mass, self.charge/self.mass

   D[1][1] = k[3]*u[3]+k[2]*u[2]+k[1]*u[1] + 0*1i 
   D[1][2] = k[1]*n + 0*1i 
   D[1][3] = k[2]*n + 0*1i 
   D[1][4] = k[3]*n + 0*1i 
   D[1][5] = 0 + 0*1i 
   D[2][1] = 0 + 0*1i 
   D[2][2] = k[3]*u[3]+k[2]*u[2]+k[1]*u[1] + 0*1i 
   D[2][3] = 0 + B[3]*qbym*1i 
   D[2][4] = 0 + -B[2]*qbym*1i 
   D[2][5] = k[1]/(m*n) + 0*1i 
   D[3][1] = 0 + 0*1i 
   D[3][2] = 0 + -B[3]*qbym*1i 
   D[3][3] = k[3]*u[3]+k[2]*u[2]+k[1]*u[1] + 0*1i 
   D[3][4] = 0 + B[1]*qbym*1i 
   D[3][5] = k[2]/(m*n) + 0*1i 
   D[4][1] = 0 + 0*1i 
   D[4][2] = 0 + B[2]*qbym*1i 
   D[4][3] = 0 + -B[1]*qbym*1i 
   D[4][4] = k[3]*u[3]+k[2]*u[2]+k[1]*u[1] + 0*1i 
   D[4][5] = k[3]/(m*n) + 0*1i 
   D[5][1] = 0 + 0*1i 
   D[5][2] = k[1]*p*gamma + 0*1i 
   D[5][3] = k[2]*p*gamma + 0*1i 
   D[5][4] = k[3]*p*gamma + 0*1i 
   D[5][5] = k[3]*u[3]+k[2]*u[2]+k[1]*u[1] + 0*1i 
   
   return D
end

-- Construct matrix with contributions to field part of dispersion
-- matrix: this is a numEquations X 6 sized block (there are 6 maxwell
-- equations)
function Euler:calcFieldDispMat(k, E, B)
   local D = Lin.ComplexMat(5, 6)
   matrixClear(D, 0.0)

   local u = self.velocity
   local qbym = self.charge/self.mass

   D[2][1] = qbym*1i 
   D[2][5] = -u[3]*qbym*1i 
   D[2][6] = u[2]*qbym*1i 
   D[3][2] = qbym*1i 
   D[3][4] = u[3]*qbym*1i 
   D[3][6] = -u[1]*qbym*1i 
   D[4][3] = qbym*1i 
   D[4][4] = -u[2]*qbym*1i 
   D[4][5] = u[1]*qbym*1i
 
   return D
end

-- Ten-moment species type
local TenMoment = Proto(Species)

function TenMoment:init(tbl)
   TenMoment.super.init(self, tbl)

   self.density = tbl.density -- number density
   
   assert(#tbl.velocity == 3, "Must provide all 3 components of velocity")
   self.velocity = tbl.velocity -- velocity vector

   assert(#tbl.pressure == 6, "All six components of pressure (Pxx, Pxy, Pxz, Pyy, Pyz, Pzz) must be provided")
   self.pressureTensor = tbl.pressueTensor -- pressure
end

-- Number of equations in system
function TenMoment:numEquations()
   return 10
end

-- Field equations
local Maxwell = Proto()

function Maxwell:init(tbl)
   self.epsilon0, self.mu0 = tbl.epsilon0, tbl.mu0

   assert(#tbl.electricField == 3, "Must all 3 components of equilibrium electric field")
   self.electricFld = tbl.electricField

   assert(#tbl.magneticField == 3, "Must all 3 components of equilibrium magnetic field")
   self.magneticFld = tbl.magneticField
end

-- Number of equations in system
function Maxwell:numEquations()
   return 6
end

-- Calculate contribution to field part of dispersion matrix
function Maxwell:calcFieldDispMat(k)
   local D = Lin.ComplexMat(6, 6)
   matrixClear(D, 0.0)
   
   local c = 1/math.sqrt(self.epsilon0*self.mu0)

   D[1][5] = k[3]*c^2
   D[1][6] = -k[2]*c^2
   D[2][4] = -k[3]*c^2 
   D[2][6] = k[1]*c^2 
   D[3][4] = k[2]*c^2 
   D[3][5] = -k[1]*c^2
   D[4][2] = -k[3]
   D[4][3] = k[2]
   D[5][1] = k[3]
   D[5][3] = -k[1] 
   D[6][1] = -k[2]
   D[6][2] = k[1]

   return D
end

-- Calculate contribution to moment part of dispersion matrix
function Maxwell:calcMomDispMat(k, q, n, u)
   -- 3 components of current, contributing to n and u[1..3]
   local D = Lin.ComplexMat(3, 4)
   matrixClear(D, 0.0)
   
   local epsilon0 = self.epsilon0

   D[1][1] = -(u[1]*q)/epsilon0*1i 
   D[1][2] = -(n*q)/epsilon0*1i 
   D[1][3] = 0*1i 
   D[1][4] = 0*1i 
   D[2][1] = -(u[2]*q)/epsilon0*1i 
   D[2][2] = 0*1i 
   D[2][3] = -(n*q)/epsilon0*1i 
   D[2][4] = 0*1i 
   D[3][1] = -(u[3]*q)/epsilon0*1i 
   D[3][2] = 0*1i 
   D[3][3] = 0*1i 
   D[3][4] = -(n*q)/epsilon0*1i

   return D
end

return {
   IsoThermal = IsoThermal,
   Euler = Euler,
   TenMoment = TenMoment,
   Maxwell = Maxwell
}