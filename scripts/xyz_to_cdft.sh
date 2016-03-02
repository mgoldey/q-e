#!/bin/bash

if [ -z "$2" ];then
  echo ""
  echo "create a directory with CDFT run files given an xyz file"
  echo ""
  echo "\tsh scprt.sh system.xyz acc_start acc_end donor_start donor_end pseudo_dir"
  echo ""
  echo "example : "
  echo "sh xyz_to_cdft.sh 1.xyz 1 5 6 10 ../../pseudos/"
  echo ""
  echo "or : "
  echo "sh xyz_to_cdft.sh 1.xyz ../../pseudos/"
  echo ""
  exit
fi

# get number of atoms, name of system and working directory
nat=`head -n1 $1`
name=${1%.xyz}
dir=`pwd`

# acc/donor string
# if less than 6 args assume dimer
if [ -z "$6" ];then

get_addstring(){
python - <<END
print "1 ",int($1/2.0),int($1/2.0+1),$1
END
}
adstring=$(get_addstring $nat)
psdir=$2

else

adstring="$2 $3 $4 $5"
psdir=$6

fi

# create directory
mkdir -p $name


# when called print the atomic spec card with all the pp data
get_atom_spec_card(){
python - <<END
import glob, os

olddir = os.getcwd()

print "ATOMIC_SPECIES"

fil=open('$1','r')
fil.readline()
fil.readline()
dat=[l.split() for l in fil]

# print elms only once
seen=[]
for l in dat:
  if l[0] not in seen:
    # get PP file for that elm
    os.chdir("$2")
    for file in glob.glob(l[0]+".*.UPF"):
      pp=file
      os.chdir(olddir)

    # print elm and PP
    print '{} 1.0 {}'.format(l[0], pp)
    seen.append(l[0])

END
}

# get ntyp
get_ntyp(){
python - <<END
import glob, os

olddir = os.getcwd()


fil=open('$1','r')
fil.readline()
fil.readline()
dat=[l.split() for l in fil]

# print elms only once
seen=[]
ntyp=0
for l in dat:
  if l[0] not in seen:

    # add new elms to ntyp
    ntyp=ntyp+1
    seen.append(l[0])

print ntyp
END
}
ntyp=$( get_ntyp $1 )

## get num of bands
cat > tmp_ecount.sh <<EOF
#!/bin/bash
# sh script.sh xyz.in
#

#vars
xyzfil=\$1
mydir=\$(cd \`dirname \$0\` && pwd)
ppdir="${dir}/$psdir"
etot=0

#functions
getz(){
atomz=\`awk "/z_valence/{print}" \$pp\` ; atomz=\${atomz#*\\"} ; atomz=\${atomz%%\\.*} 
}

#main

#remove tabs
sed -i "s/\\t/ /g" \$xyzfil

while read line;
do

  atom="\${line%%[0-9]*}"
  atom="\${atom%%\\t*}"
  atom="\${atom%%\\ *}"

  # if atom is there get pp file name
  if [ ! -z "\$atom" ]
  then
    pp=\`ls \$ppdir/\${atom}.*UPF\`
  fi

  # if pp is there add z electrons to etot
  if [ -e "\$pp" ]
  then
    getz
    etot=\$((\$etot + \$atomz))
  fi

done < \$xyzfil

echo \$etot
EOF

numelec=`sh tmp_ecount.sh $1`

# add 6 unocc states
numbnd(){
python - <<END
print int($1/2.0 + 6)
END
}
nbnd=$(numbnd $numelec)
rm tmp_ecount.sh

# get ntyp

# get system size
# size = 4*D+7ang
# where D is greatest dist 
# from center of mass
getsize(){
python - <<END
f=open('$1','r')
nat=float(f.readline())
f.readline()
lines=f.readlines()

cm = [0,0,0]
for line in lines:
  at,x,y,z = line.split()
  fx=float(x)
  fy=float(y)
  fz=float(z)
  cm[0] = cm[0] + fx
  cm[1] = cm[1] + fy
  cm[2] = cm[2] + fz

cm = [ c / nat for c in cm ]

rmax = 0.0
for line in lines:
  at,x,y,z = line.split()
  fx=float(x)-cm[0]
  fy=float(y)-cm[1]
  fz=float(z)-cm[2]
  r = (fx**2+fy**2+fz**2)**0.5
  if r > rmax :
    rmax = r

diam=rmax*2.0
print int(diam*2+7)
END
}
size=$(getsize $1)

#move sys to center of cell
movetocent(){
python - <<END
f=open('$1','r')
nat=float(f.readline())
f.readline()
lines=f.readlines()

cm = [0,0,0]
for line in lines:
  at,x,y,z = line.split()
  fx=float(x)
  fy=float(y)
  fz=float(z)
  cm[0] = cm[0] + fx
  cm[1] = cm[1] + fy
  cm[2] = cm[2] + fz

cm = [ c / nat for c in cm ]

cc = $2/2.0

#move coords to center of cell
rmax = 0.0
for line in lines:
  at,x,y,z = line.split()
  fx=float(x)+(cc-cm[0]) 
  fy=float(y)+(cc-cm[1])
  fz=float(z)+(cc-cm[2])
  print at,fx,fy,fz
END
}

# copy xyz to new dir and go to new dir
cp $1 $name
cd $name

# create left input
cat > left.in << EOF
&CONTROL 
  outdir = 'out'  
  prefix = 'left' 
  calculation  = "scf",
  pseudo_dir   = "./",
  wf_collect   = .TRUE.
  verbosity    = "high"
/
&SYSTEM
  ibrav     = 1,
  a         = $size,
  nat       = $nat,
  ntyp      = $ntyp,
  nbnd      = $nbnd,
  ecutwfc   = 90,
  assume_isolated = 'mt',
  tot_charge = +1,
  nspin     = 2,
  tot_magnetization = 1
  nosym = .true.
/
&ELECTRONS 
  electron_maxstep = 300
  conv_thr    = 1.D-10,
  diago_full_acc = .true.
! mixing_beta = 0.50,
/
&IONS
/
&CELL
/
EPCDFT 
 1 1.0E-06 0.05D0 50
 delta_charge $adstring 1.0 0.10 
K_POINTS {Gamma}
EOF

get_atom_spec_card ${1##*/} $dir/$psdir >> left.in
echo "ATOMIC_POSITIONS (angstrom)" >> left.in
movetocent $1 $size >> left.in # move to center of cell
#tail -n${nat} ${1##*/} >> left.in

# create right input
cat > right.in << EOF
&CONTROL 
  outdir = 'out'  
  prefix = 'right' 
  calculation  = "scf",
  pseudo_dir   = "./",
  wf_collect   = .TRUE.
  verbosity    = "high"
/
&SYSTEM
  ibrav     = 1,
  a         = $size,
  nat       = $nat,
  ntyp      = $ntyp,
  nbnd      = $nbnd,
  ecutwfc   = 90,
  assume_isolated = 'mt',
  tot_charge = +1,
  nspin     = 2,
  tot_magnetization = 1
  nosym = .true.
/
&ELECTRONS 
  electron_maxstep = 300
  conv_thr    = 1.D-10,
  diago_full_acc = .true.
! mixing_beta = 0.50,
/
&IONS
/
&CELL
/
EPCDFT 
 1 1.0E-06 0.05D0 50 
 delta_charge $adstring -1.0 -0.10 
K_POINTS {Gamma}
EOF

get_atom_spec_card ${1##*/} $dir/$psdir >> right.in
echo "ATOMIC_POSITIONS (angstrom)" >> right.in
movetocent $1 $size >> right.in # move to center of cell
#tail -n${nat} ${1##*/} >> right.in # dont move to center

# create edison run file
cat > runqe << EOF
#!/bin/bash -l
#SBATCH -p regular
#SBATCH --qos=premium
#SBATCH -N 8
#SBATCH -t 02:00:00
#SBATCH -J $name
#
# sbatch myscript.slurm
#---------------------------------
#
pwrun='/global/homes/n/nbrawand/src/epcdft/espresso-5.1.1/bin/pw.x' 
pprun='/global/homes/n/nbrawand/src/epcdft/espresso-5.1.1/bin/epcdft_coupling.x'

srun -n 192 \$pwrun < left.in >& left.out
srun -n 192 \$pwrun < right.in >& right.out
sh /global/homes/n/nbrawand/src/epcdft/scripts/setup_coupling_input.sh left.out right.out >& coupling.in
srun -n 24 \$pprun < coupling.in >& coupling.out
EOF

# copy pseudos to dir
get_atom_spec_card ${1##*/} $dir/$psdir > tmp_ps_list
while read line
do
 pp=`echo $line | awk '{print $3}'`
 if [ ! -z "$pp" ];then
 cp $dir/$psdir/$pp ./
 fi
done<tmp_ps_list
rm tmp_ps_list

# go back to working directory
cd $dir
