#!/bin/csh -f 
#foreach norbits (1 10 100 1000 10000)
set k1 = 0.1
set k2 = 1
#foreach norbits ( 10 100 1000 10000 100000 1000000 )
foreach norbits (10 100 1000 10000 )
#foreach type (0 1 2 3 4 )
foreach type (5 6)
if ($type == 0) then
  set typename = rk4
else if  ($type == 1) then
  set typename = lf
else if  ($type == 2) then
  set typename = y4
else if  ($type == 3) then
  set typename = y6a
else if  ($type == 4) then
  set typename = y6b
else if  ($type == 5) then
  set typename = trp
else if  ($type == 6) then
  set typename = y4i
endif
set logfile = nonlinear_o${norbits}${typename}.dat
rm  $logfile
foreach nsteps (16 32 64 128 256 512 1024 2048  )
set interval = `echo $nsteps $norbits | awk '{print $1*$2/10}'`
echo $norbits $nsteps $interval
echo $k1 $k2 $nsteps $norbits  $interval $type | nonlinear5 >> $logfile
end
end
end




