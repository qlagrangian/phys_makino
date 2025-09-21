set terminal post eps 22
set output "harmfig2.ps"
set xlabel "steps"
set ylabel "error"
set data style lines
set zero 1e-100
set title "RK"
set logscale xy
set yrange [1e-15:1e-2]
set xrange [50:100000]
set ytics (1e-14, 1e-12, 1e-10, 1e-8, 1e-6,  1e-4,  1e-2)

plot "harmonic4a_rk4_o1.dat" using 4:7  t "Norb=1" with l 1, \
 "harmonic4a_rk4_o100.dat" using 4:7  t "Norb=100" with l 3, \
 "harmonic4a_rk4_o10000.dat" using 4:7  t "Norb=10000" with l 5
