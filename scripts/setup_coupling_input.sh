#!/bin/bash
#
# creates an input file for epcdft_coupling.x
# from 2 cdft runs left.out & right.out
#
# sh scrpt.sh left.out right.out > coupling.in
#
#
thef1=`awk '/! /{print $5}' $1 | tail -n 1`
thec1=`awk '/correction /{print $5}' $1 | tail -n 1`
#
thef2=`awk '/! /{print $5}' $2 | tail -n 1`
thec2=`awk '/correction /{print $5}' $2 | tail -n 1`
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
  free1 = $thef1
  free2 = $thef2
  cor1 =  $thec1 
  cor2 =  $thec2
/
EOF
#
cat txt
