-- Gkyl ------------------------------------------------------------------------
--
-- Test for updater to interpolate a field from one grid to another.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Unit       = require "Unit"
local Grid       = require "Grid"
local DataStruct = require "DataStruct"
local Basis      = require "Basis"
local Updater    = require "Updater"

local assert_equal = Unit.assert_equal
local assert_close = Unit.assert_close
local stats        = Unit.stats

local function createGrid(lo,up,nCells)
   local gridOut = Grid.RectCart {
      lower = lo,
      upper = up,
      cells = nCells,
   }
   return gridOut
end

local function createBasis(dim, pOrder, bKind)
   local basis
   if (bKind=="Ser") then
      basis = Basis.CartModalSerendipity { ndim = dim, polyOrder = pOrder }
   elseif (bKind=="Max") then
      basis = Basis.CartModalMaxOrder { ndim = dim, polyOrder = pOrder }
   elseif (bKind=="Tensor") then
      basis = Basis.CartModalTensor { ndim = dim, polyOrder = pOrder }
   else
      assert(false,"Invalid basis")
   end
   return basis
end

local function createField(grid, basis, vComp)
   vComp = vComp or 1
   local fld = DataStruct.Field {
      onGrid        = grid,
      numComponents = basis:numBasis()*vComp,
      ghost         = {0, 0},
   }
   return fld
end

local function createProject(grid,basis)
   local projUp = Updater.ProjectOnBasis {
      onGrid   = grid,
      basis    = basis,
      evaluate = function(t, xn) return 1 end
   }
   return projUp
end

local function intDGmom(gridIn,basisIn,momIn,addBasisIn)
   addBasisIn = false or addBasisIn
   local intMom
   if addBasisIn then
      intMom = Updater.IntegratedDGMoment {
         onGrid    = gridIn,
         basis     = basisIn,
         confBasis = addBasisIn,
         moment    = momIn,
      }
   else
      intMom = Updater.IntegratedDGMoment {
         onGrid = gridIn,
         basis  = basisIn,
         moment = momIn,
      }
   end
   return intMom
end

local function interpDG(gridIn,basisIn,gridOut,basisOut)
   local interpUp = Updater.CartFieldInterpolate {
      onGrid    = gridOut,
      onBasis   = basisOut,
      fromGrid  = gridIn,
      fromBasis = basisIn,
   }
   return interpUp
end

-- Function to project and transfer between 1D grids.
local function sampleFunc1D(x, aIn, c0In, c1In)
   return (x[1]^2)/2.0-aIn*(x[1]^4)/12.0+c0In*x[1]+c1In
end

-- Function to project and transfer between 2D grids.
local function sampleFunc2D(xIn, aIn, bIn, cIn, dIn)
   local x = xIn[1]
   local y = xIn[2]
   return ((x^2)/2.0+cIn[1]*x+cIn[2])*((y^2)/2.0+dIn[1]*y+dIn[2])
end

local function test_1x(pOrder, basisName, writeOut)
   local lower = {0.0}
   local upper = {1.0}

   local maxPower = 3
   local powBases = {2, 3}
   -- List all the possible grids we will consider.
   local possibleNs = {} 
   local nI = 0
   for cI = 1, (maxPower+1)^(#powBases) do
      local pow1 = math.floor(math.ceil(cI/(maxPower+1))-1)
      local pow2 = math.floor(((cI-1+(maxPower+1)) % (maxPower+1)))
      if not ((pow1==0) and (pow2==0)) then
         nI = nI+1
         possibleNs[nI] = math.floor(powBases[1]^pow1)*math.floor(powBases[2]^pow2)
      end
   end
   -- Generate the permutations P(n,2) where n=#possibleNs.
   local pairsN = {}
   for i = 1,#possibleNs do
      for j = 1,#possibleNs do
         pairsN[(i-1)*#possibleNs+j] = {{possibleNs[i]}, {possibleNs[j]}}
      end
   end

   for _, N in ipairs(pairsN) do
      for i = 1,2 do   -- Test interpolation in both directions.
         if i==1 then
            fromIdx = 1
            toIdx   = 2
         else
            fromIdx = 2
            toIdx   = 1
         end
         -- Grids and basis.
         local grid  = {createGrid(lower, upper, N[fromIdx]), createGrid(lower, upper, N[toIdx])}
         local basis = createBasis(grid[1]:ndim(), pOrder, basisName)
         -- Fields.
         local fld = {createField(grid[1],basis), createField(grid[2],basis)}
      
         -- Project a sample function on the field to interpolate from.
         local project = createProject(grid[1],basis)
         project:setFunc(function (t, xn)
               a  = 2.0
               c0 = (a/12.0-0.5)
               c1 = 0.0
               return sampleFunc1D(xn,a,c0,c1) 
            end)
         project:advance(0.0, {}, {fld[1]})
      
         -- Create interpolation operator and perform interpolation.
         local interpolation = interpDG(grid[1],basis,grid[2],basis)
         interpolation:advance(0.0,{fld[1]},{fld[2]})
      
         -- Compute the first 3 integrated moments of the fields.
         local intF     = {DataStruct.DynVector {numComponents = 1,}, DataStruct.DynVector {numComponents = 1,}}
         local intFcalc = {intDGmom(grid[1],basis,"one"), intDGmom(grid[2],basis,"one")}
         intFcalc[1]:advance(0.0,{fld[1]},{intF[1]})
         intFcalc[2]:advance(0.0,{fld[2]},{intF[2]})
         local jnk
         local fMom = {}
         jnk, fMom[1] = intF[1]:lastData()
         jnk, fMom[2] = intF[2]:lastData()
      
         local intxF     = {DataStruct.DynVector {numComponents = 1,}, DataStruct.DynVector {numComponents = 1,}}
         local intxFcalc = {intDGmom(grid[1],basis,"x1"), intDGmom(grid[2],basis,"x1")}
         intxFcalc[1]:advance(0.0,{fld[1]},{intxF[1]})
         intxFcalc[2]:advance(0.0,{fld[2]},{intxF[2]})
         local jnk
         local xfMom = {}
         jnk, xfMom[1] = intxF[1]:lastData()
         jnk, xfMom[2] = intxF[2]:lastData()
      
         local intxSqF     = {DataStruct.DynVector {numComponents = 1,}, DataStruct.DynVector {numComponents = 1,}}
         local intxSqFcalc = {intDGmom(grid[1],basis,"xSq"), intDGmom(grid[2],basis,"xSq")}
         intxSqFcalc[1]:advance(0.0,{fld[1]},{intxSqF[1]})
         intxSqFcalc[2]:advance(0.0,{fld[2]},{intxSqF[2]})
         local jnk
         local xSqfMom = {}
         jnk, xSqfMom[1] = intxSqF[1]:lastData()
         jnk, xSqfMom[2] = intxSqF[2]:lastData()

         -- Write data for plotting
         if writeOut then
            fld[1]:write(string.format("inFld_%sN%d-N%d_p%d.bp",basisName,N[fromIdx][1],N[toIdx][1],pOrder),0.0)
            fld[2]:write(string.format("outFld_%sN%d-N%d_p%d.bp",basisName,N[fromIdx][1],N[toIdx][1],pOrder),0.0)
         end
      
         assert_equal(fMom[1][1], fMom[2][1], string.format("Checking 'one' moment in N%d-N%d",N[fromIdx][1],N[toIdx][1]))
         assert_equal(xfMom[1][1], xfMom[2][1], string.format("Checking 'x_1' moment in N%d-N%d",N[fromIdx][1],N[toIdx][1]))
         if (pOrder >= 2) then
            assert_equal(xSqfMom[1][1], xSqfMom[2][1], string.format("Checking 'x^2' moment in N%d-N%d",N[fromIdx][1],N[toIdx][1]))
         else
            assert_close(xSqfMom[1][1], xSqfMom[2][1], 0.01*xSqfMom[1][1], string.format("Checking 'x^2' moment in N%d-N%d",N[fromIdx][1],N[toIdx][1]))
         end
      end
   end

end

local function test_2x(pOrder, basisName, writeOut, testN, off)
   off   = off or 0
   testN = testN or 500

   local lower = {0.0, 0.0}
   local upper = {1.0, 1.0}

   local maxPower = 2
   local powBases = {2, 3}
   -- List all the possible grids we will consider.
   local possibleNs = {}
   local nI = 0
   for cI = 1, (maxPower+1)^(#powBases) do
      local pow1 = math.floor(math.ceil(cI/(maxPower+1))-1)
      local pow2 = math.floor(((cI-1+(maxPower+1)) % (maxPower+1)))
      if not ((pow1==0) and (pow2==0)) then
         nI = nI+1
         possibleNs[nI] = math.floor(powBases[1]^pow1)*math.floor(powBases[2]^pow2)
      end
   end
   -- Generate the permutations P(n,2) where n=#possibleNs.
   local pairsN = {}
   for i = 1,#possibleNs do for j = 1,#possibleNs do
      pairsN[(i-1)*#possibleNs+j] = {possibleNs[i], possibleNs[j]}
   end end
   -- Generate the permutations of grid-pair combinations.
   local gridPairs = {}
   for i = 1,#pairsN do for j = 1,#pairsN do
      gridPairs[(i-1)*#pairsN+j] = {pairsN[i], pairsN[j]}
   end end

   testNpairs = testN
   offPairs   = off
   if (#gridPairs)*2 < 1000 then
      testNpairs = #gridPairs
      offPairs   = 0
   elseif ((off+1)*testN > #gridPairs) then
      testNpairs = #gridPairs - off*testN
   end

   print(" ")
   print(string.format(" Number of interpolations possible: %d",(#gridPairs)*2))
   print(string.format(" Number of possible grid combinations: %d",#gridPairs))
   print(string.format(" Testing grid combinations %d to %d.",offPairs*testN+1,offPairs*testN+testNpairs))
   print(" ")

   for k = 1, testNpairs do
      N = gridPairs[offPairs*testN+k]
      for i = 1,2 do   -- Test interpolation in both directions.
         if i==1 then
            fromIdx = 1
            toIdx   = 2
         else
            fromIdx = 2
            toIdx   = 1
         end
         -- Grids and basis.
         local grid  = {createGrid(lower, upper, N[fromIdx]), createGrid(lower, upper, N[toIdx])}
         local basis = createBasis(grid[1]:ndim(), pOrder, basisName)
         -- Fields.
         local fld = {createField(grid[1],basis), createField(grid[2],basis)}

         -- Project a sample function on the field to interpolate from.
         local project = createProject(grid[1],basis)
         project:setFunc(function (t, xn)
               local a = 2.0
               local b = 5.0
               local c = {a/12.0 - 0.5, 0.0}
               local d = {0.0, b/12.0 - 0.5}
               return sampleFunc2D(xn,a,b,c,d)
            end)
         project:advance(0.0, {}, {fld[1]})

         -- Create interpolation operator and perform interpolation.
         local interpolation = interpDG(grid[1],basis,grid[2],basis)
         interpolation:advance(0.0,{fld[1]},{fld[2]})

         -- Compute the first 3 integrated moments of the fields.
         local intF     = {DataStruct.DynVector {numComponents = 1,}, DataStruct.DynVector {numComponents = 1,}}
         local intFcalc1 = intDGmom(grid[1],basis,"one")
         local intFcalc2 = intDGmom(grid[2],basis,"one")
         intFcalc1:advance(0.0,{fld[1]},{intF[1]})
         intFcalc2:advance(0.0,{fld[2]},{intF[2]})
         local jnk
         local fMom = {}
         jnk, fMom[1] = intF[1]:lastData()
         jnk, fMom[2] = intF[2]:lastData()

         local intxF     = {DataStruct.DynVector {numComponents = grid[1]:ndim(),}, DataStruct.DynVector {numComponents = grid[2]:ndim(),}}
         local intxFcalc1 = intDGmom(grid[1],basis,"xi")
         local intxFcalc2 = intDGmom(grid[2],basis,"xi")
         intxFcalc1:advance(0.0,{fld[1]},{intxF[1]})
         intxFcalc2:advance(0.0,{fld[2]},{intxF[2]})
         local jnk
         local xfMom = {}
         jnk, xfMom[1] = intxF[1]:lastData()
         jnk, xfMom[2] = intxF[2]:lastData()

         local intxSqF     = {DataStruct.DynVector {numComponents = 1,}, DataStruct.DynVector {numComponents = 1,}}
         local intxSqFcalc1 = intDGmom(grid[1],basis,"xSq")
         local intxSqFcalc2 = intDGmom(grid[2],basis,"xSq")
         intxSqFcalc1:advance(0.0,{fld[1]},{intxSqF[1]})
         intxSqFcalc2:advance(0.0,{fld[2]},{intxSqF[2]})
         local jnk
         local xSqfMom = {}
         jnk, xSqfMom[1] = intxSqF[1]:lastData()
         jnk, xSqfMom[2] = intxSqF[2]:lastData()

         -- Write data for plotting
         if writeOut then
            fld[1]:write(string.format("inFld_%sN%dx%d-N%dx%d_p%d.bp",basisName,N[fromIdx][1],N[fromIdx][2],N[toIdx][1],N[toIdx][2],pOrder),0.0)
            fld[2]:write(string.format("outFld_%sN%dx%d-N%dx%d_p%d.bp",basisName,N[fromIdx][1],N[fromIdx][2],N[toIdx][1],N[toIdx][2],pOrder),0.0)
         end

         assert_equal(fMom[1][1], fMom[2][1], string.format("Checking 'one' moment in N%dx%d-N%dx%d",N[fromIdx][1],N[fromIdx][2],N[toIdx][1],N[toIdx][2]))
         assert_equal(xfMom[1][1], xfMom[2][1], string.format("Checking 'x_1' moment in N%dx%d-N%dx%d",N[fromIdx][1],N[fromIdx][2],N[toIdx][1],N[toIdx][2]))
         assert_equal(xfMom[1][2], xfMom[2][2], string.format("Checking 'x_2' moment in N%dx%d-N%dx%d",N[fromIdx][1],N[fromIdx][2],N[toIdx][1],N[toIdx][2]))
         if (pOrder >= 2) then
            assert_equal(xSqfMom[1][1], xSqfMom[2][1], string.format("Checking 'x^2' moment in N%dx%d-N%dx%d",N[fromIdx][1],N[fromIdx][2],N[toIdx][1],N[toIdx][2]))
         else
            assert_close(xSqfMom[1][1], xSqfMom[2][1], 0.02*xSqfMom[1][1], string.format("Checking 'x^2' moment in N%dx%d-N%dx%d",N[fromIdx][1],N[fromIdx][2],N[toIdx][1],N[toIdx][2]))
         end
      end
   end
end

local polyOrder    = 1
local basisType    = "Ser"
local writeFiles   = false

-- These are for 2D and higher dim.
local numInterpMax = 200   -- Maximum number of different interpolations to test (keep <201).
local interpOffset = 0     -- Offset to decide which interpolations to test (start with 0).

-- Due to the high number of possible interpolations, need to only perform one test at a time.
--test_1x( polyOrder, basisType, writeFiles)
test_2x( polyOrder, basisType, writeFiles, numInterpMax,  interpOffset )


if stats.fail > 0 then
   print(string.format("\nPASSED %d tests", stats.pass))
   print(string.format("**** FAILED %d tests", stats.fail))
else
   print(string.format("PASSED ALL %d tests!", stats.pass))
end
