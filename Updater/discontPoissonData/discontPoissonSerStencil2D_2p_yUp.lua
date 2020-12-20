local Lin = require "Lib.Linalg"
local function stencilFn(dx, D, N, val)
  local _M = {}

  _M[1] = Lin.Mat(8,8)
  _M[1][1][1] = (135.0*N+30.0*D)/(36.0*dx[2]^2*N+4.0*dx[2]^2*D)
  _M[1][1][2] = 0.0
  _M[1][1][3] = (342.9460598986376*N+81.40638795573723*D)/(72.0*dx[2]^2*N+8.0*dx[2]^2*D)
  _M[1][1][4] = 0.0
  _M[1][1][5] = 0.0
  _M[1][1][6] = (965.9813662799093*N+261.6199533674754*D)/(360.0*dx[2]^2*N+40.0*dx[2]^2*D)
  _M[1][1][7] = 0.0
  _M[1][1][8] = 0.0
  _M[1][2][1] = 0.0
  _M[1][2][2] = (135.0*N+30.0*D)/(36.0*dx[2]^2*N+4.0*dx[2]^2*D)
  _M[1][2][3] = 0.0
  _M[1][2][4] = (342.9460598986376*N+81.40638795573723*D)/(72.0*dx[2]^2*N+8.0*dx[2]^2*D)
  _M[1][2][5] = 0.0
  _M[1][2][6] = 0.0
  _M[1][2][7] = 0.0
  _M[1][2][8] = (965.9813662799089*N+261.6199533674754*D)/(360.0*dx[2]^2*N+40.0*dx[2]^2*D)
  _M[1][3][1] = -(1.0*(36.37306695894642*N-1.732050807568877*D))/(9.0*dx[2]^2*N+dx[2]^2*D)
  _M[1][3][2] = 0.0
  _M[1][3][3] = -(1.0*(687.0*N-57.0*D))/(144.0*dx[2]^2*N+16.0*dx[2]^2*D)
  _M[1][3][4] = 0.0
  _M[1][3][5] = 0.0
  _M[1][3][6] = -(1.0*(1591.796155291248*N-298.2197176579711*D))/(720.0*dx[2]^2*N+80.0*dx[2]^2*D)
  _M[1][3][7] = 0.0
  _M[1][3][8] = 0.0
  _M[1][4][1] = 0.0
  _M[1][4][2] = -(1.0*(36.37306695894642*N-1.732050807568877*D))/(9.0*dx[2]^2*N+dx[2]^2*D)
  _M[1][4][3] = 0.0
  _M[1][4][4] = -(1.0*(687.0*N-57.0*D))/(144.0*dx[2]^2*N+16.0*dx[2]^2*D)
  _M[1][4][5] = 0.0
  _M[1][4][6] = 0.0
  _M[1][4][7] = 0.0
  _M[1][4][8] = -(1.0*(1591.796155291248*N-298.2197176579711*D))/(720.0*dx[2]^2*N+80.0*dx[2]^2*D)
  _M[1][5][1] = 0.0
  _M[1][5][2] = 0.0
  _M[1][5][3] = 0.0
  _M[1][5][4] = 0.0
  _M[1][5][5] = (135.0*N+30.0*D)/(36.0*dx[2]^2*N+4.0*dx[2]^2*D)
  _M[1][5][6] = 0.0
  _M[1][5][7] = (1714.730299493189*N+407.0319397786862*D)/(360.0*dx[2]^2*N+40.0*dx[2]^2*D)
  _M[1][5][8] = 0.0
  _M[1][6][1] = (80.49844718999243*N+20.12461179749811*D)/(18.0*dx[2]^2*N+2.0*dx[2]^2*D)
  _M[1][6][2] = 0.0
  _M[1][6][3] = (755.2317525104463*N+213.014084041408*D)/(144.0*dx[2]^2*N+16.0*dx[2]^2*D)
  _M[1][6][4] = 0.0
  _M[1][6][5] = 0.0
  _M[1][6][6] = (333.0*N+129.0*D)/(144.0*dx[2]^2*N+16.0*dx[2]^2*D)
  _M[1][6][7] = 0.0
  _M[1][6][8] = 0.0
  _M[1][7][1] = 0.0
  _M[1][7][2] = 0.0
  _M[1][7][3] = 0.0
  _M[1][7][4] = 0.0
  _M[1][7][5] = (34.64101615137755*D-727.4613391789285*N)/(180.0*dx[2]^2*N+20.0*dx[2]^2*D)
  _M[1][7][6] = 0.0
  _M[1][7][7] = -(1.0*(687.0*N-57.0*D))/(144.0*dx[2]^2*N+16.0*dx[2]^2*D)
  _M[1][7][8] = 0.0
  _M[1][8][1] = 0.0
  _M[1][8][2] = (160.9968943799848*N+40.24922359499621*D)/(36.0*dx[2]^2*N+4.0*dx[2]^2*D)
  _M[1][8][3] = 0.0
  _M[1][8][4] = (755.2317525104464*N+213.0140840414079*D)/(144.0*dx[2]^2*N+16.0*dx[2]^2*D)
  _M[1][8][5] = 0.0
  _M[1][8][6] = 0.0
  _M[1][8][7] = 0.0
  _M[1][8][8] = (333.0*N+129.0*D)/(144.0*dx[2]^2*N+16.0*dx[2]^2*D)
  _M[2] = Lin.Mat(8,8)
  _M[2][1][1] = -(1.0*(135.0*N+102.0*D))/(36.0*dx[2]^2*N+4.0*dx[2]^2*D)
  _M[2][1][2] = 0.0
  _M[2][1][3] = (342.9460598986376*N-122.9756073373903*D)/(72.0*dx[2]^2*N+8.0*dx[2]^2*D)
  _M[2][1][4] = 0.0
  _M[2][1][5] = 0.0
  _M[2][1][6] = -(1.0*(965.9813662799093*N+1670.342779192343*D))/(360.0*dx[2]^2*N+40.0*dx[2]^2*D)
  _M[2][1][7] = 0.0
  _M[2][1][8] = 0.0
  _M[2][2][1] = 0.0
  _M[2][2][2] = -(1.0*(135.0*N+102.0*D))/(36.0*dx[2]^2*N+4.0*dx[2]^2*D)
  _M[2][2][3] = 0.0
  _M[2][2][4] = (342.9460598986376*N-122.9756073373903*D)/(72.0*dx[2]^2*N+8.0*dx[2]^2*D)
  _M[2][2][5] = 0.0
  _M[2][2][6] = 0.0
  _M[2][2][7] = 0.0
  _M[2][2][8] = -(1.0*(965.9813662799089*N+1670.342779192343*D))/(360.0*dx[2]^2*N+40.0*dx[2]^2*D)
  _M[2][3][1] = (36.37306695894642*N-29.44486372867091*D)/(9.0*dx[2]^2*N+dx[2]^2*D)
  _M[2][3][2] = 0.0
  _M[2][3][3] = -(1.0*(2097.0*N+729.0*D))/(144.0*dx[2]^2*N+16.0*dx[2]^2*D)
  _M[2][3][4] = 0.0
  _M[2][3][5] = 0.0
  _M[2][3][6] = -(1.0*(848.1833528194243*N+4907.069899644797*D))/(720.0*dx[2]^2*N+80.0*dx[2]^2*D)
  _M[2][3][7] = 0.0
  _M[2][3][8] = 0.0
  _M[2][4][1] = 0.0
  _M[2][4][2] = (36.37306695894642*N-29.44486372867091*D)/(9.0*dx[2]^2*N+dx[2]^2*D)
  _M[2][4][3] = 0.0
  _M[2][4][4] = -(1.0*(2097.0*N+729.0*D))/(144.0*dx[2]^2*N+16.0*dx[2]^2*D)
  _M[2][4][5] = 0.0
  _M[2][4][6] = 0.0
  _M[2][4][7] = 0.0
  _M[2][4][8] = -(1.0*(848.1833528194243*N+4907.069899644797*D))/(720.0*dx[2]^2*N+80.0*dx[2]^2*D)
  _M[2][5][1] = 0.0
  _M[2][5][2] = 0.0
  _M[2][5][3] = 0.0
  _M[2][5][4] = 0.0
  _M[2][5][5] = -(1.0*(135.0*N+102.0*D))/(36.0*dx[2]^2*N+4.0*dx[2]^2*D)
  _M[2][5][6] = 0.0
  _M[2][5][7] = (1714.730299493189*N-614.8780366869514*D)/(360.0*dx[2]^2*N+40.0*dx[2]^2*D)
  _M[2][5][8] = 0.0
  _M[2][6][1] = -(1.0*(160.9968943799849*N+147.5804865149861*D))/(36.0*dx[2]^2*N+4.0*dx[2]^2*D)
  _M[2][6][2] = 0.0
  _M[2][6][3] = -(1.0*(398.9172846593639*D-731.9938524332015*N))/(144.0*dx[2]^2*N+16.0*dx[2]^2*D)
  _M[2][6][4] = 0.0
  _M[2][6][5] = 0.0
  _M[2][6][6] = -(1.0*(6003.0*N+1599.0*D))/(144.0*dx[2]^2*N+16.0*dx[2]^2*D)
  _M[2][6][7] = 0.0
  _M[2][6][8] = 0.0
  _M[2][7][1] = 0.0
  _M[2][7][2] = 0.0
  _M[2][7][3] = 0.0
  _M[2][7][4] = 0.0
  _M[2][7][5] = (181.8653347947321*N-147.2243186433546*D)/(45.0*dx[2]^2*N+5.0*dx[2]^2*D)
  _M[2][7][6] = 0.0
  _M[2][7][7] = -(1.0*(2097.0*N+729.0*D))/(144.0*dx[2]^2*N+16.0*dx[2]^2*D)
  _M[2][7][8] = 0.0
  _M[2][8][1] = 0.0
  _M[2][8][2] = -(1.0*(160.9968943799848*N+147.5804865149861*D))/(36.0*dx[2]^2*N+4.0*dx[2]^2*D)
  _M[2][8][3] = 0.0
  _M[2][8][4] = (731.9938524332018*N-398.9172846593639*D)/(144.0*dx[2]^2*N+16.0*dx[2]^2*D)
  _M[2][8][5] = 0.0
  _M[2][8][6] = 0.0
  _M[2][8][7] = 0.0
  _M[2][8][8] = -(1.0*(6003.0*N+1599.0*D))/(144.0*dx[2]^2*N+16.0*dx[2]^2*D)
  _M[3] = Lin.Vec(8)
  _M[3][1] = 36.0/(9.0*dx[2]^2*N+dx[2]^2*D)*val
  _M[3][2] = 0.0*val
  _M[3][3] = 55.42562584220407/(9.0*dx[2]^2*N+dx[2]^2*D)*val
  _M[3][4] = 0.0*val
  _M[3][5] = 0.0*val
  _M[3][6] = 53.66563145999496/(9.0*dx[2]^2*N+dx[2]^2*D)*val
  _M[3][7] = 0.0*val
  _M[3][8] = 0.0*val
  return(_M)
end

return(stencilFn)