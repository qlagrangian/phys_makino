set terminal post eps 22
set output "harmfig1.ps"
set xlabel "t"
set ylabel "error"
set data style lines
set zero 1e-100
set title "Norb=10, Nsteps=1000, RK"
plot "harmn1000o10rk.dat" using 2:4  t "x" with l 1 , "harmn1000o10rk.dat" using 2:6 t "v" with l 3

