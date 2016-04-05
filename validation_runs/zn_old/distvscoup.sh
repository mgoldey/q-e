#!/bin/bash
for ff in [5-9]/coupling.out
 do 
 energy1=`awk '/ Ha/{print $3}' $ff`
 dist=${ff%%/c*}
 echo "$dist $energy1"
done
