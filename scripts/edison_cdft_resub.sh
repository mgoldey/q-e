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
	echo "RUN SCRIPT section needs to be modded before use."
	echo ""
	exit 1
}

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
resubjob="sbatch $newrun"
runresub="sbatch $resubrunfil"

# this time must be less than walltime\n
resubtime=" sleep 4m ; exit 0 ; \n"


head="
#!/bin/bash -l\n
#SBATCH -p debug\n
#SBATCH --qos=premium\n
#SBATCH -N 2\n
#SBATCH -t 00:10:00\n
#SBATCH -J $name \n
#\n
# sbatch myscript.slurm\n
#---------------------------------\n
#\n
pwrun='/global/homes/n/nbrawand/src/epcdft/espresso-5.1.1/bin/pw.x' \n
pprun='/global/homes/n/nbrawand/src/epcdft/espresso-5.1.1/bin/epcdft_coupling.x'\n
\n
{ \n
"

mkcin="sh /global/homes/n/nbrawand/src/epcdft/scripts/setup_coupling_input.sh $lif $rif >& coupling.in"
runc="srun -n 24 \$pprun < coupling.in >& $cf"
runright="srun -n 48 \$pwrun < $rif >& $rf"
runleft="srun -n 48 \$pwrun < $lif >& $lf"

resubhead="
#!/bin/bash -l\n
#SBATCH -p regular\n
#SBATCH --qos=normal\n
#SBATCH -N 1\n
#SBATCH -t 00:02:00\n
#SBATCH -J $name  \n
#SBATCH -d afterany:wjobid \n
#\n
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
		jobid=`echo $jobid | awk '{print $4}'`
		echo $jobid

		# create resub script and submit
		echo -e $resubhead > $resubrunfil
		chmod a+x $resubrunfil
		echo "sh $0 $@" >> $resubrunfil
		sed -i "s/wjobid/$jobid/" $resubrunfil
		sed -i "s/ #/#/g" $resubrunfil # removing extra white space
		eval $runresub
fi
