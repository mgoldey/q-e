#!/bin/bash
#
# creates an input file for pp.x
# from a QE run
#
# sh script.sh out1 out2 out3 out4
#

for i in $@
do
i1=`echo $i | sed 's/\.out/\.in/'`
stem=`echo $i1 | sed 's/\.in//'`
sdd=${stem}_sdd.in
#
cat>$sdd<< EOF
&INPUTPP
EOF

grep prefix $i1 >>$sdd
if [ `grep -c outdir $i1` == 1 ]
then
grep outdir $i1 >>$sdd
else
echo "outdir='./'" >>$sdd
fi
cat >> $sdd <<EOF
filplot = '${stem}.sdd'
plot_num = 6
/
&PLOT
iflag = 3
output_format = 6
fileout = '${stem}.cube'
/
EOF
#
done

