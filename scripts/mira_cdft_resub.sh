#!/bin/bash
#
# Check status of cdft run and resubmit if not done 
# RUN SCRIPT section needs to be modded before use
#
# The RUN SCRIPT section is setup for mira/cetus
#
# To start job run this script in directory with
# input files 
#
# this script will submit the cdft job and a 
# resubmit script which will wait for the 
# cdft job to exit and then resubmit this script
#
# Nicholas Brawand nicholasbrawand@gmail.com

#
# parse args
#
usage() { 
	echo ""
	echo "Usage: $0 [-L <left.in>] [-l <left.out>] [-R <right.in>] [-r <right.out>] [-c <coupling.out>] [-n <jobname>]" 1>&2
	echo ""
	echo "If no options given the default values above are used."
	echo ""
	echo "RUN SCRIPT section needs to be modded before use."
	echo ""
	exit 1
}

# defualt values
lf="left.out"
lif="left.in"
rf="right.out"
rif="right.in"
cf="coupling.out"
dir=`pwd`
name="${dir##*/}"


while getopts ":l:L:r:R:c:n:" o; do
    case "${o}" in
        l) # left output file
            lf=${OPTARG}
            ;;
        L) # left input file
            lif=${OPTARG}
            if [ ! -f $lif ] ; then 
		echo "run script does not exist"; exit 1 
	    fi
            ;;
        r) # right output file
            rf=${OPTARG}
            ;;
        R) # right input file
            rif=${OPTARG}
            if [ ! -f $rif ] ; then 
		echo "run script does not exist"; exit 1 
 	    fi
            ;;
        c) # coupling output file
            cf=${OPTARG}
            ;;
        n) # jobname
            name=${OPTARG}
            ;;
        *)
            usage
            ;;
    esac
done

# if any required input is missing print usage and exit
if [ -z "${lf}" ] || [ -z "${lif}" ] || [ -z "${rif}" ] || [ -z "${rf}" ] || [ -z "${cf}" ] ; then
	usage
	exit 1
fi

# if name is blank
if [ -z ${name// } ] ; then 
    name="cdft_job"
fi

# if these are false then change intput and
# rerun calc
ldone=false
rdone=false
cdone=false
newrun="run_${name}"
resubrunfil="run_resub"


#
#===== RUN SCRIPT =====
#
# run file info this section needs to be 
# modded based on runscrpt
resubjob="qsub $newrun"
runresub="qsub $resubrunfil"
#qoe="ECSInterface"
qoe="CompSpectro"

# this time must be less than walltime\n
resubtime=" sleep 57m ; exit 0 ; \n"


head="
#!/bin/bash -x\n
#COBALT -n 128\n
#COBALT -t 01:00:00\n
#COBALT -A $qoe\n
#COBALT --jobname=$name\n
#\n
\n
pwpream=\"--block \$COBALT_PARTNAME -p 16 --envs OMP_NUM_THREADS=1\"\n
root='/home/nbrawand/src/epcdft'\n
pwrun='/home/nbrawand/src/epcdft/espresso-5.1.1/bin/pw.x' \n
pprun='/home/nbrawand/src/epcdft/espresso-5.1.1/bin/epcdft_coupling.x'\n

{ \n
"

mkcin="bash \${root}/scripts/setup_coupling_input.sh $lf $rf >& coupling.in"
runc="runjob \$pwpream : \$pprun < coupling.in >& $cf"
runright="runjob \$pwpream : \$pwrun < $rif >& $rf"
runleft="runjob \$pwpream : \$pwrun < $lif >& $lf"

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
elif [ ! -z "`grep 'JOB DONE' $lf`" ];then #job rdy to restart
        sed -i 's/restart_mode/!restart_mode/g' $lif
        sed -i "s/CONTROL/CONTROL \n restart_mode = 'restart'/g" $lif
else #job crashed 
        sed -i 's/restart_mode/!restart_mode/g' $lif
fi


if [ ! -z "`grep 'convergence has' $rf`" ] && [ ! -z "`grep 'JOB DONE' $rf`" ]; then
	rdone=true
elif [ ! -z "`grep 'JOB DONE' $rf`" ];then #job rdy to restart
        sed -i 's/restart_mode/!restart_mode/g' $rif
        sed -i "s/CONTROL/CONTROL \n restart_mode = 'restart'/g" $rif
else #job crashed 
        sed -i 's/restart_mode/!restart_mode/g' $rif
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

		echo -e $head > $newrun

		if ! $ldone ; then
		# left not done
				echo $runleft >> $newrun
				if [ -e $lif ]; then
						# replace EPCDFT input pot
						cdft_rep $lif $lf
				fi
		fi
		
		if ! $rdone ; then
		# right is done
				echo $runright >> $newrun
				if [ -e $rif ]; then
						cdft_rep $rif $rf
				fi
		fi
		
		if ! $cdone ; then
		# coupling not done
				echo $mkcin >> $newrun
				echo $runc >> $newrun
		fi

		# submit job
		sed -i "s/ #/#/g" $newrun # removing extra white space
		chmod a+x $newrun
		echo "exit 0 ;" >> $newrun
		echo "} &" >> $newrun
		echo -e $resubtime >> $newrun
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
