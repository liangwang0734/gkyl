-- Gkyl ------------------------------------------------------------------------
--
-- CartField I/O using ADIOS format
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- system libraries
local ffi = require "ffi"
local xsys = require "xsys"
local new, copy, fill, sizeof, typeof, metatype = xsys.from(ffi,
     "new, copy, fill, sizeof, typeof, metatype")

-- Gkyl libraries
local Mpi = require "Comm.Mpi"
local Adios = require "Io.Adios"
local Alloc = require "Lib.Alloc"
local Proto = require "Lib.Proto"

-- Code from Lua wiki to convert table to comma-seperated-values
-- string.
-- Used to escape "'s by toCSV
local function escapeCSV (s)
  if string.find(s, '[,"]') then
    s = '"' .. string.gsub(s, '"', '""') .. '"'
  end
  return s
end
-- Convert from table to CSV string
local function toCSV (tt)
  local s = ""
  -- ChM 23.02.2014: changed pairs to ipairs assumption is that
  -- fromCSV and toCSV maintain data as ordered array
  for _,p in ipairs(tt) do  
    s = s .. "," .. escapeCSV(p)
  end
  return string.sub(s, 2)      -- remove first comma
end

-- AdiosCartFieldIo ------------------------------------------------------------
--
-- CartField I/O using ADIOS
--------------------------------------------------------------------------------

-- CartField I/O using ADIOS
local AdiosCartFieldIo = Proto()

-- constructor to make a new uniform grid
function AdiosCartFieldIo:init(tbl)
   -- by default, write out doubles
   local elct = tbl.elemType and tbl.elemType or typeof("double")
   self._elemType = elct -- element type stored in field
   self._method = tbl.method and tbl.method or "MPI"

   -- set ADIOS data-types
   self._elctIoType = Adios.double
   if ffi.istype(new(elct), new("double")) then
      self._elctIoType = Adios.double
   elseif ffi.istype(new(elct), new("float")) then
      self._elctIoType = Adios.real
      isNumberType = true
   elseif ffi.istype(new(elct), new("int")) then
      self._elctIoType = Adios.integer
   else
      self._elctIoType = Adios.byte
   end

   -- create memory allocator
   self._allocator = Alloc.Alloc_meta_ctor(elct)
   -- allocate memory buffer for use in ADIOS I/O
   self._outBuff = self._allocator(1) -- this will be resized on an actual write()
end

-- Writes field to file.
-- fName: file name
-- tmStamp: time-stamp
-- frNum: frame number
function AdiosCartFieldIo:write(field, fName, tmStamp, frNum)
   local comm =  Mpi.getComm(field:grid():commSet().nodeComm)
   -- (the extra getComm() is needed as Lua has no concept of
   -- pointers and hence we don't know before hand if nodeComm is a
   -- pointer or an object)

   -- no need to do anything if communicator is not valid
   if not Mpi.Is_comm_valid(comm) then return end

   local ndim = field:ndim()
   local localRange, globalRange = field:localRange(), field:globalRange()
   
   -- for use in ADIOS output
   local _adLocalSz, _adGlobalSz, _adOffset = {}, {}, {}
   for d = 1, ndim do
      _adLocalSz[d] = localRange:shape(d)
      _adGlobalSz[d] = globalRange:shape(d)
      _adOffset[d] = localRange:lower(d)-1
   end
   _adLocalSz[ndim+1] = field:numComponents()
   _adGlobalSz[ndim+1] = field:numComponents()
   _adOffset[ndim+1] = 0

   -- convert tables to comma-seperated-string. For some strange
   -- reasons, this is what ADIOS expects
   adLocalSz = toCSV(_adLocalSz)
   adGlobalSz = toCSV(_adGlobalSz)
   adOffset = toCSV(_adOffset)

   if not frNum then frNum = 5000 end  -- default frame-number
   if not tmStamp then tmStamp = 0.0 end -- default time-stamp

   -- resize buffer (only done if needed. Alloc handles this automatically)
   self._outBuff:expand(field:size())

   local rank = Mpi.Comm_rank(comm)

   -- setup ADIOS for IO
   Adios.init_noxml(comm)
   --Adios.set_max_buffer_size(16) -- 16 MB chunks	 

   -- setup group and set I/O method
   local grpId = Adios.declare_group("CartField", "", Adios.flag_no)
   Adios.select_method(grpId, self._method, "", "")

   -- field attributes
   local cells = new("int[?]", ndim)
   for d = 1, ndim do cells[d-1] = globalRange:shape(d) end
   Adios.define_attribute_byvalue(grpId, "numCells", "", Adios.integer, ndim, cells)

   local lower = new("double[?]", ndim)
   for d = 1, ndim do lower[d-1] = field:grid():lower(d) end
   Adios.define_attribute_byvalue(grpId, "lowerBounds", "", Adios.double, ndim, lower)

   local upper = new("double[?]", ndim)
   for d = 1, ndim do upper[d-1] = field:grid():upper(d) end
   Adios.define_attribute_byvalue(grpId, "upperBounds", "", Adios.double, ndim, upper)

   -- define data to write
   Adios.define_var(
      grpId, "frame", "", Adios.integer, "", "", "")
   Adios.define_var(
      grpId, "time", "", Adios.double, "", "", "")
   Adios.define_var(
      grpId, "CartGridField", "", self._elctIoType, adLocalSz, adGlobalSz, adOffset)

   -- copy field into output buffer (this copy is needed as
   -- field also contains ghost-cell data, and, in addition,
   -- ADIOS expects data to be laid out in row-major order)
   field:_copy_from_field_region(field:localRange(), self._outBuff)

   local fullNm = GKYL_OUT_PREFIX .. "_" .. fName -- concatenate prefix
   -- open file to write out group
   local fd = Adios.open("CartField", fullNm, "w", comm)

   local tmStampBuff = new("double[1]"); tmStampBuff[0] = tmStamp
   Adios.write(fd, "time", tmStampBuff)

   local frNumBuff = new("int[1]"); frNumBuff[0] = frNum
   Adios.write(fd, "frame", frNumBuff)

   Adios.write(fd, "CartGridField", self._outBuff:data())
   Adios.close(fd)
   
   Adios.finalize(rank)      
end

-- Read field from file.
-- fName: file name
-- returns time-stamp, frame number
function AdiosCartFieldIo:read(field, fName)
   local comm =  Mpi.getComm(field:grid():commSet().nodeComm)
   -- (the extra getComm() is needed as Lua has no concept of
   -- pointers and hence we don't know before hand if nodeComm is a
   -- pointer or an object)

   -- no need to do anything if communicator is not valid
   if not Mpi.Is_comm_valid(comm) then return end

   local ndim = field:ndim()
   local localRange, globalRange = field:localRange(), field:globalRange()
   
   -- for use in ADIOS output
   local _adLocalSz, _adGlobalSz, _adOffset = {}, {}, {}
   for d = 1, ndim do
      _adLocalSz[d] = localRange:shape(d)
      _adGlobalSz[d] = globalRange:shape(d)
      _adOffset[d] = localRange:lower(d)-1
   end
   _adLocalSz[ndim+1] = field:numComponents()
   _adGlobalSz[ndim+1] = field:numComponents()
   _adOffset[ndim+1] = 0

   -- convert tables to comma-seperated-string. For some strange
   -- reasons, this is what ADIOS expects
   adLocalSz = toCSV(_adLocalSz)
   adGlobalSz = toCSV(_adGlobalSz)
   adOffset = toCSV(_adOffset)

   -- resize buffer (only done if needed. Alloc handles this automatically)
   self._outBuff:expand(field:size())

   local rank = Mpi.Comm_rank(comm)

   -- setup ADIOS for IO
   Adios.init_noxml(comm)
   --Adios.set_max_buffer_size(16) -- 16 MB chunks	 

   -- setup group and set I/O method
   local grpId = Adios.declare_group("CartField", "", Adios.flag_no)
   Adios.select_method(grpId, self._method, "", "")

   -- field attributes
   local cells = new("int[?]", ndim)
   for d = 1, ndim do cells[d-1] = globalRange:shape(d) end
   Adios.define_attribute_byvalue(grpId, "numCells", "", Adios.integer, ndim, cells)

   local lower = new("double[?]", ndim)
   for d = 1, ndim do lower[d-1] = field:grid():lower(d) end
   Adios.define_attribute_byvalue(grpId, "lowerBounds", "", Adios.double, ndim, lower)

   local upper = new("double[?]", ndim)
   for d = 1, ndim do upper[d-1] = field:grid():upper(d) end
   Adios.define_attribute_byvalue(grpId, "upperBounds", "", Adios.double, ndim, upper)

   -- define data to write
   Adios.define_var(
      grpId, "frame", "", Adios.integer, "", "", "")
   Adios.define_var(
      grpId, "time", "", Adios.double, "", "", "")
   Adios.define_var(
      grpId, "CartGridField", "", self._elctIoType, adLocalSz, adGlobalSz, adOffset)

   local fullNm = GKYL_OUT_PREFIX .. "_" .. fName -- concatenate prefix
   -- open file to write out group
   local fd = Adios.open("CartField", fullNm, "r", comm)

   local tmStampBuff = new("double[1]")
   local frNumBuff = new("int[1]")

   Adios.read(fd, "time", tmStampBuff, sizeof("double"))
   Adios.read(fd, "frame", frNumBuff, sizeof("int"))

   local fieldLocalSz = localRange:volume()*field:numComponents()*sizeof("double")
   Adios.read(fd, "CartGridField", self._outBuff:data(), fieldLocalSz)
   Adios.close(fd)

   -- copy output buffer into field (this copy is needed as field also
   -- contains ghost-cell data, and, in addition, ADIOS expects data
   -- to be laid out in row-major order)
   field:_copy_to_field_region(field:localRange(), self._outBuff)
   
   Adios.finalize(rank)

   return tmStampBuff[0], frNumBuff[0]
end

return AdiosCartFieldIo
