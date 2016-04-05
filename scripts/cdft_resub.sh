#!/bin/bash
#
# Check status of cdft run and resubmit if not done 
# RUN SCRIPT section needs to be modded before use
#
# Nicholas Brawand nicholasbrawand@gmail.com

#
# parse args
#
usage() { 
	echo ""
	echo "Usage: $0 [-L <left.in>] [-l <left.out>] [-R <right.in>] [-r <right.out>] [-c <coupling.out>] [-s <run script>]" 1>&2
	echo "RUN SCRIPT section needs to be modded before use."
	echo ""
	exit 1
}

while getopts ":l:L:r:R:c:s:" o; do
    case "${o}" in
        l)# left output file
            lf=${OPTARG}
            #if [ ! -f $lf ] ; then 
			#	echo "left out file does not exist"; exit 1 
			#fi
            ;;
        L)# left input file
            lif=${OPTARG}
            if [ ! -f $lif ] ; then 
				echo "run script does not exist"; exit 1 
			fi
            ;;
        r)# right output file
            rf=${OPTARG}
            #if [ ! -f $rf ] ; then 
			#	echo "right out file does not exist"; exit 1 
			#fi
            ;;
        R)# right input file
            rif=${OPTARG}
            if [ ! -f $rif ] ; then 
				echo "run script does not exist"; exit 1 
			fi
            ;;
        c)# coupling output file
            cf=${OPTARG}
            #if [ ! -f $cf ] ; then 
			#	echo "coupling out file does not exist"; exit 1 
			#fi
            ;;
        s)# run script
            sf=${OPTARG}
            if [ ! -f $sf ] ; then 
				echo "run script does not exist"; exit 1 
			fi
            ;;
        *)
            usage
            ;;
    esac
done

# if any required input is missing print usage and exit
if [ -z "${lf}" ] || [ -z "${lif}" ] || [ -z "${rif}" ] || [ -z "${rf}" ] || [ -z "${cf}" ] || [ -z "${sf}" ]; then
    usage
	exit 1
fi


# if these are false then change intput and
# rerun calc
ldone=false
rdone=false
cdone=false
newsf="${sf}_current"
resubrunfil="run_resub"


#
#===== RUN SCRIPT =====
#
# run file info this section needs to be 
# modded based on runscrpt
resubjob="qsub $newsf"
runresub="qsub $resubrunfil"

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

mkcin="sh \${root}/epcdft/scripts/setup_coupling_input.sh left.out right.out >& coupling.in"
runc="runjob \$pwpream : \$pprun < coupling.in >& coupling.out"
runright="runjob \$pwpream : \$pwrun < right.in >& right.out"
runleft="runjob \$pwpream : \$pwrun < left.in >& left.out"

resubhead="
#!/bin/bash\n
#COBALT -n 1\n
#COBALT -t 00:05:00\n
#COBALT -A $qoe\n
#COBALT --jobname=r_$name\n
#COBALT --dependencies=wjobid
\n
"
#
#===== END RUN SCRIPT =====
#

# commands to check if calcs are complete 
if [ ! -z "`grep 'convergence has' $lf`" ] && [ ! -z "`grep 'JOB DONE' $lf`" ]; then
	ldone=true
fi

if [ ! -z "`grep 'convergence has' $rf`" ] && [ ! -z "`grep 'JOB DONE' $rf`" ]; then
	rdone=true
fi

if [ ! -z "`grep 'JOB DONE' $cf`" ]; then
	cdone=true
fi

cdft_rep(){
# this function replaces input V's
# with new V's from output files
# of incomplete CDFT runs
# $1 input file
# $2 output file

	line=`grep -a2 EPCDFT $1 | tail -n1 | awk '{$NF=""; print $0}'` 
	vold=`grep -a2 EPCDFT $1 | tail -n1 | awk '{print $NF}'`
	vnew=`grep -a1 New $2 | tail -n1 | awk '{print $NF}'`

	if [ ! -z "$line" ] && [ ! -z "$vold" ] && [ ! -z "$vnew" ]; then
		sed -i "s/$line *$vold/$line $vnew/" $1
	fi

}

#
# resubmit or complete job
#

if $ldone && $rdone && $cdone ; then
	# if job is done
	exit 1 # DONE

else # job NOT DONE create new run file

		echo -e $head > $newsf

		if ! $ldone ; then
		# left not done
				echo $runleft >> $newsf
				if [ -e $lif ]; then
						# replace EPCDFT input pot
						cdft_rep $lif $lf
				fi
		fi
		
		if ! $rdone ; then
		# right is done
				echo $runright >> $newsf
				if [ -e $rif ]; then
						cdft_rep $rif $rf
				fi
		fi
		
		if ! $cdone ; then
		# coupling not done
				echo $mkcin >> $newsf
				echo $runc >> $newsf
		fi

		# submit job
		sed -i "s/ #/#/g" $newsf # removing extra white space
		chmod a+x $newsf
		jobid=`eval $resubjob`
		echo $jobid

		# create resub script and submit
		echo -e $resubhead > $resubrunfil
		chmod a+x $resubrunfil
		echo "sh $0 $@" >> $resubrunfil
		sed -i "s/wjobid/$jobid/" $resubrunfil
		sed -i "s/ #/#/g" $resubrunfil # removing extra white space
		eval $runresub
fi
