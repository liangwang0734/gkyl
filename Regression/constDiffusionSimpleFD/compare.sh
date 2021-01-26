
pgkyl -f 'rt-testConstDiffusionSimple-1d_q*_2.bp' plot -f0

pgkyl -f 'rt-testConstDiffusionSimple-2d_q*_2.bp' select --z0 0 plot -f0

pgkyl -f 'rt-testConstDiffusionSimple-3d_q*_2.bp' select --z0 16 select --z1 16 plot -f0
