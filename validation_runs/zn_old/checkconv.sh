#!/bin/bash
for ff in [5-9]/*.out 
 do echo $ff 
 grep "convergence has" $ff
 grep JOB $ff 
 echo "" 
done
