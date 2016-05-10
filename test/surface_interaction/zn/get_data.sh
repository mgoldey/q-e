for ff in pw*.out;
do 
  pos=`grep 'reference position' $ff | tail -n1| awk '{print $6}'`
  en=`grep Image $ff | awk '{print $4}'| tail -n1`
  dip=`grep 'Total dipole moment' $ff | tail -n1 | awk '{print $6}'`
  echo $pos ${en##*:} $dip
done>data
