#!/bin/csh -f 
foreach f (nonlinear*dat)
   nawk -f nonlinear_postp.awk $f > ${f}.abs
   nawk -f calc_diff.awk $f > ${f}.diff
end
