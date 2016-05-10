pym(){
python -<<EOF
print ($1+$2)
EOF
}

for ff in pw[!_]*.out;
do 
  pos=`grep 'reference position' $ff | tail -n1| awk '{print $6}'`
  eni=`grep Image $ff | awk '{print $4}'| tail -n1`
  eni=${eni%%\[Ry\]}
  en=`grep ! $ff | awk '{print $5}'| tail -n1`
  #en=`pym $en $eni`
  en=$eni
  dip=`grep 'Total dipole moment' $ff | tail -n1 | awk '{print $6}'`
  echo $pos ${en##*:} $dip
done>datafld

#for ff in pw_*.out;
#do 
#  pos=`grep 'reference position' $ff | tail -n1| awk '{print $6}'`
#  en=`grep ! $ff | awk '{print $5}'| tail -n1`
#  dip=`grep 'Total dipole moment' $ff | tail -n1 | awk '{print $6}'`
#  echo $pos ${en##*:} $dip
#done>datanofld
