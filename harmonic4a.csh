#!/bin/csh -f
foreach itype (rk4 lf)
if ($itype == rk4) then
   set nitype = 0
else
   set nitype = 1
endif   
set logfile = harmonic4a_${itype}.dat
rm  $logfile
foreach norbits (1 10 100 1000 10000)
foreach nsteps (64 128 256 512 1024 2048 4096 8192 16384 32768 65536 )
echo $norbits $nsteps $nitype
echo  $nsteps  $norbits $nitype  | harmonic4a >> $logfile
end
end
end




