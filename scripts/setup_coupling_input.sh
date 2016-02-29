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
#
i1=`echo $1 | sed 's/\.out/\.in'`
i2=`echo $2 | sed 's/\.out/\.in'`
#
#
c1=`awk '/CDFT correction /{print $4}' $1 | tail -n 1`
c2=`awk '/CDFT correction /{print $4}' $2 | tail -n 1`
#
f1=`echo "$e1 - $c1" | bc `
f2=`echo "$e2 - $c2" | bc `
#
spinup=`awk '/\(up:  /{print $7}' $1`
spindown=`awk '/\(up:  /{print $9}' $1`
#
spinup=${spinup%%.*}
spindown=${spindown%%.*}
#
cat>txt<< EOF
&INPUTPP
EOF
grep prefix $i1 >>txt
grep outdir $i1 >>txt
grep prefix $i2 | sed 's/prefix/prefix2/' >>txt
grep outdir $i2 | sed 's/outdir/outdir2/' >>txt
cat >>txt <<EOF
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
