-- Gkyl ------------------------------------------------------------------------
--
-- Updater to compute diagnostics from a field: diagnostics can be
-- computed using functions that return a single scalar.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Gkyl libraries
local Mpi = require "Comm.Mpi"
local Proto = require "Proto"
local UpdaterBase = require "Updater.Base"
local ffi = require "ffi"

-- Updater to computer diagnostic
local CalcDiagnostic = Proto(UpdaterBase)

-- Constructor: read data from tbl
function CalcDiagnostic:init(tbl)
   CalcDiagnostic.super.init(self, tbl) -- setup base object

   self._onGrid = assert(tbl.onGrid, "Updater.CalcDiagnostic: Must provide grid object using 'onGrid'")
   -- diagnostic function
   self._diagnostic = assert(tbl.diagnostic, "Updater.CalcDiagnostic: Must specify diagnostic using 'diagnostic'")
end

-- advance method computes diagnostic across MPI ranks
function CalcDiagnostic:_advance(tCurr, dt, inFld, outFld)
   local grid = self._onGrid
   local qIn = assert(inFld[1], "CalcDiagnostic.advance: Must specify an output field")
   local dynOut = assert(outFld[1], "CalcDiagnostic.advance: Must specify an output field")

   local qItr = qIn:get(1)
   local v = 0.0
   local indexer = qIn:genIndexer()
   -- loop, computing diagnostics in each cell
   for idx in qIn:localRangeIter() do
      grid:setIndex(idx)
      local vol = grid:cellVolume()
      qIn:fill(indexer(idx), qItr)
      v = v + vol*self._diagnostic(tCurr+dt, qItr)
   end

   -- all-reduce across global communicator. NOTE: If DynVector
   -- storage changes to use shared memory on shared nodes, the below
   -- code will need modification
   local valLocal, val = ffi.new("double[1]"), ffi.new("double[1]")
   valLocal[0] = v
   Mpi.Allreduce(valLocal, val, 1, Mpi.DOUBLE, Mpi.SUM, self:getComm())
   dynOut:appendData(tCurr+dt, { val[0] })

   return true, GKYL_MAX_DOUBLE
end

return CalcDiagnostic
