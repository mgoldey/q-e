#!/bin/bash
#
# creates an input file for epcdft_coupling.x
# from 2 cdft runs left.out & right.out
#
# sh scrpt.sh left.out right.out > coupling.in
#
#
e1=`awk '/! /{print $5}' $1 | tail -n 1`
e2=`awk '/! /{print $5}' $2 | tail -n 1`
echo $1 $e1
echo $2 $e2
#
c1=`awk '/correction /{print $4}' $1 | tail -n 1`
c2=`awk '/correction /{print $4}' $2 | tail -n 1`
echo $1 $c1
echo $2 $c2
#
f1=`echo "$e1 - $c1" | bc `
f2=`echo "$e2 - $c2" | bc `
echo $1 $f1
echo $2 $f2
#
spinup=`awk '/\(up:  /{print $7}' $1`
spindown=`awk '/\(up:  /{print $9}' $1`
#
spinup=${spinup%%.*}
spindown=${spindown%%.*}
#
cat>txt<< EOF
&INPUTPP
  prefix       = "left",
  outdir       = "./out",
  prefix2      = "right",
  outdir2      = "./out",
  occup1       = $spinup
  occup2       = $spinup
  occdown1     = $spindown
  occdown2     = $spindown
  debug        = .true.
  eig_of_w     = .false.
  s_spin       = .true.
  free1 = $f1
  free2 = $f2
  cor1 =  $c1 
  cor2 =  $c2
/
EOF
#
cat txt
