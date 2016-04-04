#!/bin/bash
# Check status of cdft run and resubmit if not done 
#
# bash script.sh -l left.out -r right.out -c coupling.out
#
# Nicholas Brawand nicholasbrawand@gmail.com

# parse args
usage() { echo "Usage: $0 [-l <left.out>] [-r <right.out>] [-c <coupling.out>] [-s <run script>]" 1>&2; exit 1; }

while getopts ":l:r:c:s:" o; do
    case "${o}" in
        l)
            lf=${OPTARG}
            #if [ ! -f $lf ] ; then 
			#	echo "left out file does not exist"; exit 1 
			#fi
            ;;
        r)
            rf=${OPTARG}
            #if [ ! -f $rf ] ; then 
			#	echo "right out file does not exist"; exit 1 
			#fi
            ;;
        c)
            cf=${OPTARG}
            #if [ ! -f $cf ] ; then 
			#	echo "coupling out file does not exist"; exit 1 
			#fi
            ;;
        s)
            sf=${OPTARG}
            #if [ ! -f $sf ] ; then 
			#	echo "run script does not exist"; exit 1 
			#fi
            ;;
        *)
            usage
            ;;
    esac
done

if [ -z "${lf}" ] || [ -z "${rf}" ] || [ -z "${cf}" ] || [ -z "${sf}" ]; then
    usage
	exit 1
fi


# check status of cdft calculation and set flags
# done flags
ldone=false
rdone=false
cdone=false
newsf="${sf}_current"

# run file info this section needs to be modded based on runscrpt
resubjob="qsub $newsf"

nodes=`awk '/COBALT -n/{print $3}' $sf`
time=`awk '/COBALT -t/{print $3}' $sf`
qoe=`awk '/COBALT -A/{print $3}' $sf`
name=`grep 'COBALT --jobname' $sf`
name=${name##*=}

head="
#!/bin/bash\n
#COBALT -n $nodes\n
#COBALT -t $time\n
#COBALT -A $qoe\n
#COBALT --jobname=$name\n
#\n
\n
pwpream='--block \$COBALT_PARTNAME -p 16 --envs OMP_NUM_THREADS=1'\n
root='/home/nbrawand/src/epcdft'\n
pwrun='/home/nbrawand/src/epcdft/espresso-5.1.1/bin/pw.x' \n
pprun='/home/nbrawand/src/epcdft/espresso-5.1.1/bin/epcdft_coupling.x'\n
"
sed -i "s/ #/#/g" $newsf # removing extra white space

mkcin="sh \${root}/epcdft/scripts/setup_coupling_input.sh left.out right.out >& coupling.in"
runc="runjob \$pwpream : \$pprun < coupling.in >& coupling.out"
runright="runjob \$pwpream : \$pwrun < right.in >& right.out"
runleft="runjob \$pwpream : \$pwrun < left.in >& left.out"


# commands to check if complete 
if [ ! -z "`grep 'convergence has' $lf`" ] && [ ! -z "`grep 'JOB DONE' $lf`" ]; then
	ldone=true
fi

if [ ! -z "`grep 'convergence has' $rf`" ] && [ ! -z "`grep 'JOB DONE' $rf`" ]; then
	rdone=true
fi

if [ ! -z "`grep 'JOB DONE' $cf`" ]; then
	cdone=true
fi

# resubmit or complete job

# if job is done
if $ldone && $rdone && $cdone ; then

		exit 1 # DONE

else # job NOT DONE create new run file

		echo $head > $newsf

		# left not done
		if ! $ldone ; then
				echo $runleft >> $newsf
		fi
		
		# right is done
		if ! $rdone ; then
				echo $runright >> $newsf
		fi
		
		# coupling not done
		if ! $cdone ; then
				echo $mkcin >> $newsf
				echo $runc >> $newsf
		fi

		echo "sh $0 $@" >> $newsf
		eval $resubjob

fi
