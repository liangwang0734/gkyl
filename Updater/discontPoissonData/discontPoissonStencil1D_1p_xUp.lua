local Lin = require "Lib.Linalg"
local function stencilFn(dx, D, N, val)
  local _M = {}

  _M[1] = Lin.Mat(2,2)
  _M[1][1][1] = (41472.0*N+3072.0*D)/(18432.0*dx[1]^2*N+4096.0*dx[1]^2*D)
  _M[1][1][2] = (39906.45060638692*N+1773.62002695053*D)/(18432.0*dx[1]^2*N+4096.0*dx[1]^2*D)
  _M[1][2][1] = -(1.0*(50548.17076809011*N+19509.82029645583*D))/(18432.0*dx[1]^2*N+4096.0*dx[1]^2*D)
  _M[1][2][2] = ((-44544.0*N)-19456.0*D)/(18432.0*dx[1]^2*N+4096.0*dx[1]^2*D)
  _M[2] = Lin.Mat(2,2)
  _M[2][1][1] = -(1.0*(41472.0*N+39936.0*D))/(18432.0*dx[1]^2*N+4096.0*dx[1]^2*D)
  _M[2][1][2] = (39906.45060638692*N-54982.22083546644*D)/(18432.0*dx[1]^2*N+4096.0*dx[1]^2*D)
  _M[2][2][1] = (50548.17076809011*N-30151.54045815902*D)/(18432.0*dx[1]^2*N+4096.0*dx[1]^2*D)
  _M[2][2][2] = ((-216576.0*N)-134144.0*D)/(18432.0*dx[1]^2*N+4096.0*dx[1]^2*D)
  _M[3] = Lin.Vec(2)
  _M[3][1] = 36.0/(12.72792206135786*dx[1]^2*N+2.828427124746191*dx[1]^2*D)*val
  _M[3][2] = 48.49742261192856/(12.72792206135786*dx[1]^2*N+2.828427124746191*dx[1]^2*D)*val
  return(_M)
end

return(stencilFn)