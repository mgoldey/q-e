#!/bin/bash
for ff in [5-9]/left.out 
 do 
 energy1=`sh gete.sh $ff`
 cor1=`sh getcor.sh $ff`
 dir=${ff%%/l*}
 sed -i "s/thef1/$energy1/" $dir/coupling.in
 sed -i "s/thec1/$cor1/" $dir/coupling.in
 echo "$ff $energy1 $cor1"
 echo "" 
done

for ff in [5-9]/right.out 
 do  
 energy2=`sh gete.sh $ff`
 cor2=`sh getcor.sh $ff`
 dir=${ff%%/r*}
 sed -i "s/thef2/$energy2/" $dir/coupling.in
 sed -i "s/thec2/$cor2/" $dir/coupling.in
 echo "$ff $energy2 $cor2"
 echo "" 
done
