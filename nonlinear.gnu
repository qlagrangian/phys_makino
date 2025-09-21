set terminal pos eps 22
set output "nlfig1.ps"
set logscale x 10
set logscale y 10
set zero 1e-100
set title "T=1E6"
set data style lines
set yrange [1e-15:1e-1]
set ytics(1e-14, 1e-12, 1e-10, 1e-8, 1e-6, 1e-4, 1e-2)
plot "nonlinear_o1000000rk4.dat.abs"  t "RK" with l 1, \
"nonlinear_o1000000lf.dat.abs"  t "LF" with l 3
