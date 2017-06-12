#!/bin/bash
#
# creates an input file for epcdft_coupling.x
# from 2 cdft runs left.out & right.out
#
# sh scrpt.sh left.out right.out > coupling.in
#
i1=`echo $1 | sed 's/\.out/\.in/'`
i2=`echo $2 | sed 's/\.out/\.in/'`
echo "&INPUTPP"

echo `grep prefix $i1 `
echo `grep prefix $i2 | sed 's/prefix/prefix2/' `

if [ `grep -c outdir $i1` == 1 ] 
then
echo `grep outdir $i1 `
else
echo `outdir='./' `
fi

if [ `grep -c outdir $i2` == 1 ] 
then
echo `grep outdir $i2 | sed 's:outdir:outdir2:g' `
else
echo `outdir2='./' `
fi
echo "/" 
