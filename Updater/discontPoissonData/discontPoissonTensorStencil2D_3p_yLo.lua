local Lin = require "Lib.Linalg"
local function stencilFn(dx, D, N, val)
  local _M = {}

  _M[1] = Lin.Vec(16)
  _M[1][1] = -60.0/(15.0*dx[2]^2*N-1.0*dx[2]^2*D)*val
  _M[1][2] = 0.0*val
  _M[1][3] = 96.99484522385713/(15.0*dx[2]^2*N-1.0*dx[2]^2*D)*val
  _M[1][4] = 0.0*val
  _M[1][5] = 0.0*val
  _M[1][6] = -107.3312629199899/(15.0*dx[2]^2*N-1.0*dx[2]^2*D)*val
  _M[1][7] = 0.0*val
  _M[1][8] = 0.0*val
  _M[1][9] = 0.0*val
  _M[1][10] = 95.24704719832526/(15.0*dx[2]^2*N-1.0*dx[2]^2*D)*val
  _M[1][11] = 0.0*val
  _M[1][12] = 0.0*val
  _M[1][13] = 0.0*val
  _M[1][14] = 0.0*val
  _M[1][15] = 0.0*val
  _M[1][16] = 0.0*val
  _M[2] = Lin.Mat(16,16)
  _M[2][1][1] = -(1.0*(2625.0*N-855.0*D))/(480.0*dx[2]^2*N-32.0*dx[2]^2*D)
  _M[2][1][2] = 0.0
  _M[2][1][3] = -(1.0*(3715.248982235241*N+1837.705906830579*D))/(480.0*dx[2]^2*N-32.0*dx[2]^2*D)
  _M[2][1][4] = 0.0
  _M[2][1][5] = 0.0
  _M[2][1][6] = -(1.0*(2985.15074996222*N-1120.270056727395*D))/(480.0*dx[2]^2*N-32.0*dx[2]^2*D)
  _M[2][1][7] = 0.0
  _M[2][1][8] = 0.0
  _M[2][1][9] = 0.0
  _M[2][1][10] = -(1.0*(9167.528292838808*N+17295.27632042923*D))/(3360.0*dx[2]^2*N-224.0*dx[2]^2*D)
  _M[2][1][11] = 0.0
  _M[2][1][12] = 0.0
  _M[2][1][13] = 0.0
  _M[2][1][14] = 0.0
  _M[2][1][15] = 0.0
  _M[2][1][16] = 0.0
  _M[2][2][1] = 0.0
  _M[2][2][2] = -(1.0*(2625.0*N-855.0*D))/(480.0*dx[2]^2*N-32.0*dx[2]^2*D)
  _M[2][2][3] = 0.0
  _M[2][2][4] = -(1.0*(11145.74694670572*N+5513.117720491736*D))/(1440.0*dx[2]^2*N-96.0*dx[2]^2*D)
  _M[2][2][5] = 0.0
  _M[2][2][6] = 0.0
  _M[2][2][7] = 0.0
  _M[2][2][8] = -(1.0*(2985.150749962219*N-1120.270056727395*D))/(480.0*dx[2]^2*N-32.0*dx[2]^2*D)
  _M[2][2][9] = 0.0
  _M[2][2][10] = 0.0
  _M[2][2][11] = 0.0
  _M[2][2][12] = 0.0
  _M[2][2][13] = -(1.0*(27502.58487851641*N+51885.82896128768*D))/(10080.0*dx[2]^2*N-672.0*dx[2]^2*D)
  _M[2][2][14] = 0.0
  _M[2][2][15] = 0.0
  _M[2][2][16] = 0.0
  _M[2][3][1] = -(1.0*(4200.223208354527*N+819.2600319800789*D))/(480.0*dx[2]^2*N-32.0*dx[2]^2*D)
  _M[2][3][2] = 0.0
  _M[2][3][3] = -(1.0*(11307.0*N-4125.0*D))/(480.0*dx[2]^2*N-32.0*dx[2]^2*D)
  _M[2][3][4] = 0.0
  _M[2][3][5] = 0.0
  _M[2][3][6] = -(1.0*(4736.658632411671*N+1173.513953900847*D))/(480.0*dx[2]^2*N-32.0*dx[2]^2*D)
  _M[2][3][7] = 0.0
  _M[2][3][8] = 0.0
  _M[2][3][9] = 0.0
  _M[2][3][10] = -(1.0*(51292.76975364071*N-32366.7321334731*D))/(3360.0*dx[2]^2*N-224.0*dx[2]^2*D)
  _M[2][3][11] = 0.0
  _M[2][3][12] = 0.0
  _M[2][3][13] = 0.0
  _M[2][3][14] = 0.0
  _M[2][3][15] = 0.0
  _M[2][3][16] = 0.0
  _M[2][4][1] = 0.0
  _M[2][4][2] = -(1.0*(4200.223208354527*N+819.2600319800789*D))/(480.0*dx[2]^2*N-32.0*dx[2]^2*D)
  _M[2][4][3] = 0.0
  _M[2][4][4] = -(1.0*(11307.0*N-4125.0*D))/(480.0*dx[2]^2*N-32.0*dx[2]^2*D)
  _M[2][4][5] = 0.0
  _M[2][4][6] = 0.0
  _M[2][4][7] = 0.0
  _M[2][4][8] = -(1.0*(4736.658632411671*N+1173.513953900847*D))/(480.0*dx[2]^2*N-32.0*dx[2]^2*D)
  _M[2][4][9] = 0.0
  _M[2][4][10] = 0.0
  _M[2][4][11] = 0.0
  _M[2][4][12] = 0.0
  _M[2][4][13] = -(1.0*(7327.538536234388*N-4623.818876210443*D))/(480.0*dx[2]^2*N-32.0*dx[2]^2*D)
  _M[2][4][14] = 0.0
  _M[2][4][15] = 0.0
  _M[2][4][16] = 0.0
  _M[2][5][1] = 0.0
  _M[2][5][2] = 0.0
  _M[2][5][3] = 0.0
  _M[2][5][4] = 0.0
  _M[2][5][5] = -(1.0*(2625.0*N-855.0*D))/(480.0*dx[2]^2*N-32.0*dx[2]^2*D)
  _M[2][5][6] = 0.0
  _M[2][5][7] = -(1.0*(18576.24491117621*N+9188.529534152894*D))/(2400.0*dx[2]^2*N-160.0*dx[2]^2*D)
  _M[2][5][8] = 0.0
  _M[2][5][9] = 0.0
  _M[2][5][10] = 0.0
  _M[2][5][11] = -(1.0*(14925.7537498111*N-5601.350283636974*D))/(2400.0*dx[2]^2*N-160.0*dx[2]^2*D)
  _M[2][5][12] = 0.0
  _M[2][5][13] = 0.0
  _M[2][5][14] = 0.0
  _M[2][5][15] = -(1.0*(45837.64146419404*N+86476.38160214614*D))/(16800.0*dx[2]^2*N-1120.0*dx[2]^2*D)
  _M[2][5][16] = 0.0
  _M[2][6][1] = -(1.0*(771.4434522374276*N-1267.850543242381*D))/(480.0*dx[2]^2*N-32.0*dx[2]^2*D)
  _M[2][6][2] = 0.0
  _M[2][6][3] = (801.7075526649351*N-3783.904729244646*D)/(480.0*dx[2]^2*N-32.0*dx[2]^2*D)
  _M[2][6][4] = 0.0
  _M[2][6][5] = 0.0
  _M[2][6][6] = -(1.0*(17355.0*N-2805.0*D))/(480.0*dx[2]^2*N-32.0*dx[2]^2*D)
  _M[2][6][7] = 0.0
  _M[2][6][8] = 0.0
  _M[2][6][9] = 0.0
  _M[2][6][10] = (82582.55769228755*N-37537.52622376708*D)/(3360.0*dx[2]^2*N-224.0*dx[2]^2*D)
  _M[2][6][11] = 0.0
  _M[2][6][12] = 0.0
  _M[2][6][13] = 0.0
  _M[2][6][14] = 0.0
  _M[2][6][15] = 0.0
  _M[2][6][16] = 0.0
  _M[2][7][1] = 0.0
  _M[2][7][2] = 0.0
  _M[2][7][3] = 0.0
  _M[2][7][4] = 0.0
  _M[2][7][5] = -(1.0*(21001.11604177265*N+4096.300159900395*D))/(2400.0*dx[2]^2*N-160.0*dx[2]^2*D)
  _M[2][7][6] = 0.0
  _M[2][7][7] = -(1.0*(11307.0*N-4125.0*D))/(480.0*dx[2]^2*N-32.0*dx[2]^2*D)
  _M[2][7][8] = 0.0
  _M[2][7][9] = 0.0
  _M[2][7][10] = 0.0
  _M[2][7][11] = -(1.0*(4736.658632411671*N+1173.513953900847*D))/(480.0*dx[2]^2*N-32.0*dx[2]^2*D)
  _M[2][7][12] = 0.0
  _M[2][7][13] = 0.0
  _M[2][7][14] = 0.0
  _M[2][7][15] = -(1.0*(36637.69268117194*N-23119.09438105221*D))/(2400.0*dx[2]^2*N-160.0*dx[2]^2*D)
  _M[2][7][16] = 0.0
  _M[2][8][1] = 0.0
  _M[2][8][2] = -(1.0*(771.4434522374274*N-1267.85054324238*D))/(480.0*dx[2]^2*N-32.0*dx[2]^2*D)
  _M[2][8][3] = 0.0
  _M[2][8][4] = (801.7075526649353*N-3783.904729244646*D)/(480.0*dx[2]^2*N-32.0*dx[2]^2*D)
  _M[2][8][5] = 0.0
  _M[2][8][6] = 0.0
  _M[2][8][7] = 0.0
  _M[2][8][8] = -(1.0*(17355.0*N-2805.0*D))/(480.0*dx[2]^2*N-32.0*dx[2]^2*D)
  _M[2][8][9] = 0.0
  _M[2][8][10] = 0.0
  _M[2][8][11] = 0.0
  _M[2][8][12] = 0.0
  _M[2][8][13] = (82582.55769228755*N-37537.52622376707*D)/(3360.0*dx[2]^2*N-224.0*dx[2]^2*D)
  _M[2][8][14] = 0.0
  _M[2][8][15] = 0.0
  _M[2][8][16] = 0.0
  _M[2][9][1] = 0.0
  _M[2][9][2] = 0.0
  _M[2][9][3] = 0.0
  _M[2][9][4] = 0.0
  _M[2][9][5] = 0.0
  _M[2][9][6] = 0.0
  _M[2][9][7] = 0.0
  _M[2][9][8] = 0.0
  _M[2][9][9] = -(1.0*(2625.0*N-855.0*D))/(480.0*dx[2]^2*N-32.0*dx[2]^2*D)
  _M[2][9][10] = 0.0
  _M[2][9][11] = 0.0
  _M[2][9][12] = -(1.0*(26006.74287564669*N+12863.94134781405*D))/(3360.0*dx[2]^2*N-224.0*dx[2]^2*D)
  _M[2][9][13] = 0.0
  _M[2][9][14] = -(1.0*(20896.05524973553*N-7841.890397091763*D))/(3360.0*dx[2]^2*N-224.0*dx[2]^2*D)
  _M[2][9][15] = 0.0
  _M[2][9][16] = -(1.0*(9167.528292838808*N+17295.27632042923*D))/(3360.0*dx[2]^2*N-224.0*dx[2]^2*D)
  _M[2][10][1] = -(1.0*(3770.195618267042*N+828.1201603632169*D))/(480.0*dx[2]^2*N-32.0*dx[2]^2*D)
  _M[2][10][2] = 0.0
  _M[2][10][3] = -(1.0*(10489.51576575392*N-4009.753733086362*D))/(480.0*dx[2]^2*N-32.0*dx[2]^2*D)
  _M[2][10][4] = 0.0
  _M[2][10][5] = 0.0
  _M[2][10][6] = -(1.0*(3922.360896195046*N+1200.964195969222*D))/(480.0*dx[2]^2*N-32.0*dx[2]^2*D)
  _M[2][10][7] = 0.0
  _M[2][10][8] = 0.0
  _M[2][10][9] = 0.0
  _M[2][10][10] = -(1.0*(49833.0*N-7383.0*D))/(480.0*dx[2]^2*N-32.0*dx[2]^2*D)
  _M[2][10][11] = 0.0
  _M[2][10][12] = 0.0
  _M[2][10][13] = 0.0
  _M[2][10][14] = 0.0
  _M[2][10][15] = 0.0
  _M[2][10][16] = 0.0
  _M[2][11][1] = 0.0
  _M[2][11][2] = 0.0
  _M[2][11][3] = 0.0
  _M[2][11][4] = 0.0
  _M[2][11][5] = -(1.0*(771.4434522374276*N-1267.850543242381*D))/(480.0*dx[2]^2*N-32.0*dx[2]^2*D)
  _M[2][11][6] = 0.0
  _M[2][11][7] = (801.7075526649353*N-3783.904729244646*D)/(480.0*dx[2]^2*N-32.0*dx[2]^2*D)
  _M[2][11][8] = 0.0
  _M[2][11][9] = 0.0
  _M[2][11][10] = 0.0
  _M[2][11][11] = -(1.0*(17355.0*N-2805.0*D))/(480.0*dx[2]^2*N-32.0*dx[2]^2*D)
  _M[2][11][12] = 0.0
  _M[2][11][13] = 0.0
  _M[2][11][14] = 0.0
  _M[2][11][15] = (82582.55769228755*N-37537.52622376706*D)/(3360.0*dx[2]^2*N-224.0*dx[2]^2*D)
  _M[2][11][16] = 0.0
  _M[2][12][1] = 0.0
  _M[2][12][2] = 0.0
  _M[2][12][3] = 0.0
  _M[2][12][4] = 0.0
  _M[2][12][5] = 0.0
  _M[2][12][6] = 0.0
  _M[2][12][7] = 0.0
  _M[2][12][8] = 0.0
  _M[2][12][9] = -(1.0*(29401.56245848169*N+5734.820223860553*D))/(3360.0*dx[2]^2*N-224.0*dx[2]^2*D)
  _M[2][12][10] = 0.0
  _M[2][12][11] = 0.0
  _M[2][12][12] = -(1.0*(11307.0*N-4125.0*D))/(480.0*dx[2]^2*N-32.0*dx[2]^2*D)
  _M[2][12][13] = 0.0
  _M[2][12][14] = -(1.0*(33156.61042688169*N+8214.59767730593*D))/(3360.0*dx[2]^2*N-224.0*dx[2]^2*D)
  _M[2][12][15] = 0.0
  _M[2][12][16] = -(1.0*(51292.76975364071*N-32366.7321334731*D))/(3360.0*dx[2]^2*N-224.0*dx[2]^2*D)
  _M[2][13][1] = 0.0
  _M[2][13][2] = -(1.0*(11310.58685480112*N+2484.36048108965*D))/(1440.0*dx[2]^2*N-96.0*dx[2]^2*D)
  _M[2][13][3] = 0.0
  _M[2][13][4] = -(1.0*(31468.54729726175*N-12029.26119925908*D))/(1440.0*dx[2]^2*N-96.0*dx[2]^2*D)
  _M[2][13][5] = 0.0
  _M[2][13][6] = 0.0
  _M[2][13][7] = 0.0
  _M[2][13][8] = -(1.0*(11767.08268858514*N+3602.892587907666*D))/(1440.0*dx[2]^2*N-96.0*dx[2]^2*D)
  _M[2][13][9] = 0.0
  _M[2][13][10] = 0.0
  _M[2][13][11] = 0.0
  _M[2][13][12] = 0.0
  _M[2][13][13] = -(1.0*(49833.0*N-7383.0*D))/(480.0*dx[2]^2*N-32.0*dx[2]^2*D)
  _M[2][13][14] = 0.0
  _M[2][13][15] = 0.0
  _M[2][13][16] = 0.0
  _M[2][14][1] = 0.0
  _M[2][14][2] = 0.0
  _M[2][14][3] = 0.0
  _M[2][14][4] = 0.0
  _M[2][14][5] = 0.0
  _M[2][14][6] = 0.0
  _M[2][14][7] = 0.0
  _M[2][14][8] = 0.0
  _M[2][14][9] = -(1.0*(14287.33267618557*N-23480.92065912238*D))/(8889.724405177025*dx[2]^2*N-592.6482936784683*dx[2]^2*D)
  _M[2][14][10] = 0.0
  _M[2][14][11] = 0.0
  _M[2][14][12] = (5611.952868654547*N-26487.33310471252*D)/(3360.0*dx[2]^2*N-224.0*dx[2]^2*D)
  _M[2][14][13] = 0.0
  _M[2][14][14] = -(1.0*(17355.0*N-2805.0*D))/(480.0*dx[2]^2*N-32.0*dx[2]^2*D)
  _M[2][14][15] = 0.0
  _M[2][14][16] = (82582.55769228755*N-37537.52622376706*D)/(3360.0*dx[2]^2*N-224.0*dx[2]^2*D)
  _M[2][15][1] = 0.0
  _M[2][15][2] = 0.0
  _M[2][15][3] = 0.0
  _M[2][15][4] = 0.0
  _M[2][15][5] = -(1.0*(18850.97809133521*N+4140.600801816085*D))/(2400.0*dx[2]^2*N-160.0*dx[2]^2*D)
  _M[2][15][6] = 0.0
  _M[2][15][7] = -(1.0*(52447.5788287696*N-20048.7686654318*D))/(2400.0*dx[2]^2*N-160.0*dx[2]^2*D)
  _M[2][15][8] = 0.0
  _M[2][15][9] = 0.0
  _M[2][15][10] = 0.0
  _M[2][15][11] = -(1.0*(3922.360896195045*N+1200.964195969222*D))/(480.0*dx[2]^2*N-32.0*dx[2]^2*D)
  _M[2][15][12] = 0.0
  _M[2][15][13] = 0.0
  _M[2][15][14] = 0.0
  _M[2][15][15] = -(1.0*(49833.0*N-7383.0*D))/(480.0*dx[2]^2*N-32.0*dx[2]^2*D)
  _M[2][15][16] = 0.0
  _M[2][16][1] = 0.0
  _M[2][16][2] = 0.0
  _M[2][16][3] = 0.0
  _M[2][16][4] = 0.0
  _M[2][16][5] = 0.0
  _M[2][16][6] = 0.0
  _M[2][16][7] = 0.0
  _M[2][16][8] = 0.0
  _M[2][16][9] = -(1.0*(3770.195618267042*N+828.1201603632169*D))/(480.0*dx[2]^2*N-32.0*dx[2]^2*D)
  _M[2][16][10] = 0.0
  _M[2][16][11] = 0.0
  _M[2][16][12] = -(1.0*(10489.51576575392*N-4009.75373308636*D))/(480.0*dx[2]^2*N-32.0*dx[2]^2*D)
  _M[2][16][13] = 0.0
  _M[2][16][14] = -(1.0*(3922.360896195045*N+1200.964195969222*D))/(480.0*dx[2]^2*N-32.0*dx[2]^2*D)
  _M[2][16][15] = 0.0
  _M[2][16][16] = -(1.0*(49833.0*N-7383.0*D))/(480.0*dx[2]^2*N-32.0*dx[2]^2*D)
  _M[3] = Lin.Mat(16,16)
  _M[3][1][1] = (2625.0*N+105.0*D)/(480.0*dx[2]^2*N-32.0*dx[2]^2*D)
  _M[3][1][2] = 0.0
  _M[3][1][3] = -(1.0*(3715.248982235241*N+188.7935380250076*D))/(480.0*dx[2]^2*N-32.0*dx[2]^2*D)
  _M[3][1][4] = 0.0
  _M[3][1][5] = 0.0
  _M[3][1][6] = (2985.15074996222*N+221.3707297724792*D)/(480.0*dx[2]^2*N-32.0*dx[2]^2*D)
  _M[3][1][7] = 0.0
  _M[3][1][8] = 0.0
  _M[3][1][9] = 0.0
  _M[3][1][10] = -(1.0*(9167.528292838808*N+1039.780265248384*D))/(3360.0*dx[2]^2*N-224.0*dx[2]^2*D)
  _M[3][1][11] = 0.0
  _M[3][1][12] = 0.0
  _M[3][1][13] = 0.0
  _M[3][1][14] = 0.0
  _M[3][1][15] = 0.0
  _M[3][1][16] = 0.0
  _M[3][2][1] = 0.0
  _M[3][2][2] = (2625.0*N+105.0*D)/(480.0*dx[2]^2*N-32.0*dx[2]^2*D)
  _M[3][2][3] = 0.0
  _M[3][2][4] = -(1.0*(11145.74694670572*N+566.3806140750228*D))/(1440.0*dx[2]^2*N-96.0*dx[2]^2*D)
  _M[3][2][5] = 0.0
  _M[3][2][6] = 0.0
  _M[3][2][7] = 0.0
  _M[3][2][8] = (2985.150749962219*N+221.3707297724792*D)/(480.0*dx[2]^2*N-32.0*dx[2]^2*D)
  _M[3][2][9] = 0.0
  _M[3][2][10] = 0.0
  _M[3][2][11] = 0.0
  _M[3][2][12] = 0.0
  _M[3][2][13] = -(1.0*(27502.58487851641*N+3119.340795745152*D))/(10080.0*dx[2]^2*N-672.0*dx[2]^2*D)
  _M[3][2][14] = 0.0
  _M[3][2][15] = 0.0
  _M[3][2][16] = 0.0
  _M[3][3][1] = (4200.223208354527*N-732.6574916016349*D)/(480.0*dx[2]^2*N-32.0*dx[2]^2*D)
  _M[3][3][2] = 0.0
  _M[3][3][3] = -(1.0*(5931.0*N-1101.0*D))/(480.0*dx[2]^2*N-32.0*dx[2]^2*D)
  _M[3][3][4] = 0.0
  _M[3][3][5] = 0.0
  _M[3][3][6] = -(1.0*(995.3567199753062*D-4736.658632411671*N))/(480.0*dx[2]^2*N-32.0*dx[2]^2*D)
  _M[3][3][7] = 0.0
  _M[3][3][8] = 0.0
  _M[3][3][9] = 0.0
  _M[3][3][10] = -(1.0*(14338.87934951682*N-3624.81737471007*D))/(3360.0*dx[2]^2*N-224.0*dx[2]^2*D)
  _M[3][3][11] = 0.0
  _M[3][3][12] = 0.0
  _M[3][3][13] = 0.0
  _M[3][3][14] = 0.0
  _M[3][3][15] = 0.0
  _M[3][3][16] = 0.0
  _M[3][4][1] = 0.0
  _M[3][4][2] = (4200.223208354527*N-732.6574916016349*D)/(480.0*dx[2]^2*N-32.0*dx[2]^2*D)
  _M[3][4][3] = 0.0
  _M[3][4][4] = -(1.0*(5931.0*N-1101.0*D))/(480.0*dx[2]^2*N-32.0*dx[2]^2*D)
  _M[3][4][5] = 0.0
  _M[3][4][6] = 0.0
  _M[3][4][7] = 0.0
  _M[3][4][8] = (4736.658632411671*N-995.3567199753062*D)/(480.0*dx[2]^2*N-32.0*dx[2]^2*D)
  _M[3][4][9] = 0.0
  _M[3][4][10] = 0.0
  _M[3][4][11] = 0.0
  _M[3][4][12] = 0.0
  _M[3][4][13] = -(1.0*(2048.41133564526*N-517.8310535300099*D))/(480.0*dx[2]^2*N-32.0*dx[2]^2*D)
  _M[3][4][14] = 0.0
  _M[3][4][15] = 0.0
  _M[3][4][16] = 0.0
  _M[3][5][1] = 0.0
  _M[3][5][2] = 0.0
  _M[3][5][3] = 0.0
  _M[3][5][4] = 0.0
  _M[3][5][5] = (2625.0*N+105.0*D)/(480.0*dx[2]^2*N-32.0*dx[2]^2*D)
  _M[3][5][6] = 0.0
  _M[3][5][7] = -(1.0*(18576.24491117621*N+943.9676901250382*D))/(2400.0*dx[2]^2*N-160.0*dx[2]^2*D)
  _M[3][5][8] = 0.0
  _M[3][5][9] = 0.0
  _M[3][5][10] = 0.0
  _M[3][5][11] = (14925.7537498111*N+1106.853648862396*D)/(2400.0*dx[2]^2*N-160.0*dx[2]^2*D)
  _M[3][5][12] = 0.0
  _M[3][5][13] = 0.0
  _M[3][5][14] = 0.0
  _M[3][5][15] = -(1.0*(45837.64146419404*N+5198.901326241921*D))/(16800.0*dx[2]^2*N-1120.0*dx[2]^2*D)
  _M[3][5][16] = 0.0
  _M[3][6][1] = (771.4434522374276*N+449.4496634774578*D)/(480.0*dx[2]^2*N-32.0*dx[2]^2*D)
  _M[3][6][2] = 0.0
  _M[3][6][3] = ((-499.6148516607568*N)-747.4857858180314*D)/(480.0*dx[2]^2*N-32.0*dx[2]^2*D)
  _M[3][6][4] = 0.0
  _M[3][6][5] = 0.0
  _M[3][6][6] = -(1.0*(645.0*N-795.0*D))/(480.0*dx[2]^2*N-32.0*dx[2]^2*D)
  _M[3][6][7] = 0.0
  _M[3][6][8] = 0.0
  _M[3][6][9] = 0.0
  _M[3][6][10] = (7613.994680849207*N-3460.906673113276*D)/(3360.0*dx[2]^2*N-224.0*dx[2]^2*D)
  _M[3][6][11] = 0.0
  _M[3][6][12] = 0.0
  _M[3][6][13] = 0.0
  _M[3][6][14] = 0.0
  _M[3][6][15] = 0.0
  _M[3][6][16] = 0.0
  _M[3][7][1] = 0.0
  _M[3][7][2] = 0.0
  _M[3][7][3] = 0.0
  _M[3][7][4] = 0.0
  _M[3][7][5] = (21001.11604177265*N-3663.287458008176*D)/(2400.0*dx[2]^2*N-160.0*dx[2]^2*D)
  _M[3][7][6] = 0.0
  _M[3][7][7] = -(1.0*(5931.0*N-1101.0*D))/(480.0*dx[2]^2*N-32.0*dx[2]^2*D)
  _M[3][7][8] = 0.0
  _M[3][7][9] = 0.0
  _M[3][7][10] = 0.0
  _M[3][7][11] = (4736.658632411671*N-995.3567199753062*D)/(480.0*dx[2]^2*N-32.0*dx[2]^2*D)
  _M[3][7][12] = 0.0
  _M[3][7][13] = 0.0
  _M[3][7][14] = 0.0
  _M[3][7][15] = -(1.0*(10242.0566782263*N-2589.15526765005*D))/(2400.0*dx[2]^2*N-160.0*dx[2]^2*D)
  _M[3][7][16] = 0.0
  _M[3][8][1] = 0.0
  _M[3][8][2] = (771.4434522374274*N+449.4496634774577*D)/(480.0*dx[2]^2*N-32.0*dx[2]^2*D)
  _M[3][8][3] = 0.0
  _M[3][8][4] = -(1.0*(499.6148516607568*N+747.4857858180314*D))/(480.0*dx[2]^2*N-32.0*dx[2]^2*D)
  _M[3][8][5] = 0.0
  _M[3][8][6] = 0.0
  _M[3][8][7] = 0.0
  _M[3][8][8] = -(1.0*(645.0*N-795.0*D))/(480.0*dx[2]^2*N-32.0*dx[2]^2*D)
  _M[3][8][9] = 0.0
  _M[3][8][10] = 0.0
  _M[3][8][11] = 0.0
  _M[3][8][12] = 0.0
  _M[3][8][13] = (7613.994680849206*N-3460.906673113276*D)/(3360.0*dx[2]^2*N-224.0*dx[2]^2*D)
  _M[3][8][14] = 0.0
  _M[3][8][15] = 0.0
  _M[3][8][16] = 0.0
  _M[3][9][1] = 0.0
  _M[3][9][2] = 0.0
  _M[3][9][3] = 0.0
  _M[3][9][4] = 0.0
  _M[3][9][5] = 0.0
  _M[3][9][6] = 0.0
  _M[3][9][7] = 0.0
  _M[3][9][8] = 0.0
  _M[3][9][9] = (2625.0*N+105.0*D)/(480.0*dx[2]^2*N-32.0*dx[2]^2*D)
  _M[3][9][10] = 0.0
  _M[3][9][11] = 0.0
  _M[3][9][12] = -(1.0*(26006.74287564669*N+1321.554766175053*D))/(3360.0*dx[2]^2*N-224.0*dx[2]^2*D)
  _M[3][9][13] = 0.0
  _M[3][9][14] = (20896.05524973553*N+1549.595108407354*D)/(3360.0*dx[2]^2*N-224.0*dx[2]^2*D)
  _M[3][9][15] = 0.0
  _M[3][9][16] = -(1.0*(9167.528292838808*N+1039.780265248384*D))/(3360.0*dx[2]^2*N-224.0*dx[2]^2*D)
  _M[3][10][1] = (3770.195618267042*N-695.8325948099873*D)/(480.0*dx[2]^2*N-32.0*dx[2]^2*D)
  _M[3][10][2] = 0.0
  _M[3][10][3] = -(1.0*(5210.388565164791*N-1040.244682754976*D))/(480.0*dx[2]^2*N-32.0*dx[2]^2*D)
  _M[3][10][4] = 0.0
  _M[3][10][5] = 0.0
  _M[3][10][6] = -(1.0*(928.8245259466398*D-3922.360896195046*N))/(480.0*dx[2]^2*N-32.0*dx[2]^2*D)
  _M[3][10][7] = 0.0
  _M[3][10][8] = 0.0
  _M[3][10][9] = 0.0
  _M[3][10][10] = -(1.0*(1449.0*N-471.0*D))/(480.0*dx[2]^2*N-32.0*dx[2]^2*D)
  _M[3][10][11] = 0.0
  _M[3][10][12] = 0.0
  _M[3][10][13] = 0.0
  _M[3][10][14] = 0.0
  _M[3][10][15] = 0.0
  _M[3][10][16] = 0.0
  _M[3][11][1] = 0.0
  _M[3][11][2] = 0.0
  _M[3][11][3] = 0.0
  _M[3][11][4] = 0.0
  _M[3][11][5] = (771.4434522374276*N+449.4496634774578*D)/(480.0*dx[2]^2*N-32.0*dx[2]^2*D)
  _M[3][11][6] = 0.0
  _M[3][11][7] = -(1.0*(499.6148516607568*N+747.4857858180314*D))/(480.0*dx[2]^2*N-32.0*dx[2]^2*D)
  _M[3][11][8] = 0.0
  _M[3][11][9] = 0.0
  _M[3][11][10] = 0.0
  _M[3][11][11] = -(1.0*(645.0*N-795.0*D))/(480.0*dx[2]^2*N-32.0*dx[2]^2*D)
  _M[3][11][12] = 0.0
  _M[3][11][13] = 0.0
  _M[3][11][14] = 0.0
  _M[3][11][15] = (7613.994680849206*N-3460.906673113275*D)/(3360.0*dx[2]^2*N-224.0*dx[2]^2*D)
  _M[3][11][16] = 0.0
  _M[3][12][1] = 0.0
  _M[3][12][2] = 0.0
  _M[3][12][3] = 0.0
  _M[3][12][4] = 0.0
  _M[3][12][5] = 0.0
  _M[3][12][6] = 0.0
  _M[3][12][7] = 0.0
  _M[3][12][8] = 0.0
  _M[3][12][9] = (29401.56245848169*N-5128.602441211446*D)/(3360.0*dx[2]^2*N-224.0*dx[2]^2*D)
  _M[3][12][10] = 0.0
  _M[3][12][11] = 0.0
  _M[3][12][12] = -(1.0*(5931.0*N-1101.0*D))/(480.0*dx[2]^2*N-32.0*dx[2]^2*D)
  _M[3][12][13] = 0.0
  _M[3][12][14] = (33156.61042688169*N-6967.497039827143*D)/(3360.0*dx[2]^2*N-224.0*dx[2]^2*D)
  _M[3][12][15] = 0.0
  _M[3][12][16] = -(1.0*(14338.87934951682*N-3624.817374710069*D))/(3360.0*dx[2]^2*N-224.0*dx[2]^2*D)
  _M[3][13][1] = 0.0
  _M[3][13][2] = (11310.58685480112*N-2087.497784429962*D)/(1440.0*dx[2]^2*N-96.0*dx[2]^2*D)
  _M[3][13][3] = 0.0
  _M[3][13][4] = (1040.244682754976*D-5210.38856516479*N)/(480.0*dx[2]^2*N-32.0*dx[2]^2*D)
  _M[3][13][5] = 0.0
  _M[3][13][6] = 0.0
  _M[3][13][7] = 0.0
  _M[3][13][8] = (11767.08268858514*N-2786.473577839919*D)/(1440.0*dx[2]^2*N-96.0*dx[2]^2*D)
  _M[3][13][9] = 0.0
  _M[3][13][10] = 0.0
  _M[3][13][11] = 0.0
  _M[3][13][12] = 0.0
  _M[3][13][13] = -(1.0*(1449.0*N-471.0*D))/(480.0*dx[2]^2*N-32.0*dx[2]^2*D)
  _M[3][13][14] = 0.0
  _M[3][13][15] = 0.0
  _M[3][13][16] = 0.0
  _M[3][14][1] = 0.0
  _M[3][14][2] = 0.0
  _M[3][14][3] = 0.0
  _M[3][14][4] = 0.0
  _M[3][14][5] = 0.0
  _M[3][14][6] = 0.0
  _M[3][14][7] = 0.0
  _M[3][14][8] = 0.0
  _M[3][14][9] = (5400.104165661992*N+3146.147644342204*D)/(3360.0*dx[2]^2*N-224.0*dx[2]^2*D)
  _M[3][14][10] = 0.0
  _M[3][14][11] = 0.0
  _M[3][14][12] = ((-3497.303961625297*N)-5232.40050072622*D)/(3360.0*dx[2]^2*N-224.0*dx[2]^2*D)
  _M[3][14][13] = 0.0
  _M[3][14][14] = -(1.0*(645.0*N-795.0*D))/(480.0*dx[2]^2*N-32.0*dx[2]^2*D)
  _M[3][14][15] = 0.0
  _M[3][14][16] = (7613.994680849206*N-3460.906673113275*D)/(3360.0*dx[2]^2*N-224.0*dx[2]^2*D)
  _M[3][15][1] = 0.0
  _M[3][15][2] = 0.0
  _M[3][15][3] = 0.0
  _M[3][15][4] = 0.0
  _M[3][15][5] = (18850.97809133521*N-3479.162974049937*D)/(2400.0*dx[2]^2*N-160.0*dx[2]^2*D)
  _M[3][15][6] = 0.0
  _M[3][15][7] = -(1.0*(26051.94282582395*N-5201.223413774878*D))/(2400.0*dx[2]^2*N-160.0*dx[2]^2*D)
  _M[3][15][8] = 0.0
  _M[3][15][9] = 0.0
  _M[3][15][10] = 0.0
  _M[3][15][11] = (3922.360896195045*N-928.8245259466397*D)/(480.0*dx[2]^2*N-32.0*dx[2]^2*D)
  _M[3][15][12] = 0.0
  _M[3][15][13] = 0.0
  _M[3][15][14] = 0.0
  _M[3][15][15] = -(1.0*(1449.0*N-471.0*D))/(480.0*dx[2]^2*N-32.0*dx[2]^2*D)
  _M[3][15][16] = 0.0
  _M[3][16][1] = 0.0
  _M[3][16][2] = 0.0
  _M[3][16][3] = 0.0
  _M[3][16][4] = 0.0
  _M[3][16][5] = 0.0
  _M[3][16][6] = 0.0
  _M[3][16][7] = 0.0
  _M[3][16][8] = 0.0
  _M[3][16][9] = (3770.195618267042*N-695.8325948099873*D)/(480.0*dx[2]^2*N-32.0*dx[2]^2*D)
  _M[3][16][10] = 0.0
  _M[3][16][11] = 0.0
  _M[3][16][12] = (1040.244682754976*D-5210.38856516479*N)/(480.0*dx[2]^2*N-32.0*dx[2]^2*D)
  _M[3][16][13] = 0.0
  _M[3][16][14] = (3922.360896195045*N-928.8245259466397*D)/(480.0*dx[2]^2*N-32.0*dx[2]^2*D)
  _M[3][16][15] = 0.0
  _M[3][16][16] = -(1.0*(1449.0*N-471.0*D))/(480.0*dx[2]^2*N-32.0*dx[2]^2*D)
  return(_M)
end

return(stencilFn)