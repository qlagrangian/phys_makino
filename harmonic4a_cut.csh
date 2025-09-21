#!/bin/csh -f
foreach itype (rk4 lf)
    set logfile = harmonic4a_${itype}.dat
    foreach norbits (1 10 100 1000 10000)
	grep "norb= $norbits " $logfile > harmonic4a_${itype}_o${norbits}.dat
    end
    foreach nsteps (64 128 256 512 1024 2048 4096 8192 16384 32768 65536 )
	grep "steps= $nsteps " $logfile > harmonic4a_${itype}_s${nsteps}.dat
    end
end
