local Lin = require "Lib.Linalg"
local function stencilFn(dx, a, b, val)
  local _M = {}

  _M[1] = Lin.Vec(2)
  _M[1][1] = -33.9411254969543/(12.0*dx[1]^2*b-2.0*dx[1]^2*a)*val
  _M[1][2] = 97.97958971132718/(24.0*dx[1]^2*b-4.0*dx[1]^2*a)*val
  _M[2] = Lin.Mat(2,2)
  _M[2][1][1] = -(1.0*(27.0*b-21.0*a))/(12.0*dx[1]^2*b-2.0*dx[1]^2*a)
  _M[2][1][2] = -(1.0*(25.98076211353316*b+46.76537180435967*a))/(12.0*dx[1]^2*b-2.0*dx[1]^2*a)
  _M[2][2][1] = -(1.0*(77.94228634059947*b+34.64101615137754*a))/(24.0*dx[1]^2*b-4.0*dx[1]^2*a)
  _M[2][2][2] = -(1.0*(315.0*b-200.0*a))/(24.0*dx[1]^2*b-4.0*dx[1]^2*a)
  _M[3] = Lin.Mat(2,2)
  _M[3][1][1] = (27.0*b+3.0*a)/(12.0*dx[1]^2*b-2.0*dx[1]^2*a)
  _M[3][1][2] = -(1.0*(25.98076211353316*b+5.196152422706631*a))/(12.0*dx[1]^2*b-2.0*dx[1]^2*a)
  _M[3][2][1] = (77.94228634059947*b-34.64101615137754*a)/(24.0*dx[1]^2*b-4.0*dx[1]^2*a)
  _M[3][2][2] = -(1.0*(75.0*b-40.0*a))/(24.0*dx[1]^2*b-4.0*dx[1]^2*a)
  return(_M)
end

return(stencilFn)