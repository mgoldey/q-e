export OMP_NUM_THREADS=8

PW=~/Programs/epcdft/espresso-5.1.1/bin/pw.x
for i in *t.in
do
export j=`echo $i | sed 's/\.in/\.out/'`
if [ -e $j ]
then
continue
fi
$PW -in $i > $j
done

SCR=~/Programs/epcdft/scripts/setup_coupling_input.sh
CUP=~/Programs/epcdft/espresso-5.1.1/bin/epcdft_coupling.x

for i in *left.out
do
export j=`echo $i | sed 's/left/right/'`
export k=`echo $i | sed 's/_left.out/_coupling.in/'`
export l=`echo $i | sed 's/_left.out/_coupling.out/'`
if [ -e $l ]
then
continue
fi
$SCR $i $j > $k
$CUP -in $k > $l
done

