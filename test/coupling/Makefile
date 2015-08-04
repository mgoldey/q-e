include ../Makefile.inc
all: left.out right.out coupling.out
clean: 
	rm -rf left.out right.out out

left.out: left.in
	$(PW_CMD) < $^ > $@

right.out: right.in
	$(PW_CMD) < $^ > $@

coupling.out: coupling.in
	touch $@
#	export OMP_NUM_THREADS=1 ; mpirun -np 8 ../../espresso-5.1.1/bin//epcdft_coupling.x < $^ > $@
