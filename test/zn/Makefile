include ../Makefile.inc
all: hirshfeld.out voronoi.out zero.out 
clean: 
	rm -rf out hirshfeld.out voronoi.out zero.out well.out

hirshfeld.out: hirshfeld.in 
	$(QE_CMD) < $^ > $@

voronoi.out: voronoi.in
	$(QE_CMD) < $^ > $@

zero.out: zero.in
	$(QE_CMD) < $^ > $@

well.out: well.in
	$(QE_CMD) < $^ > $@

#coupling.out: coupling.in
#	export OMP_NUM_THREADS=1 ; mpirun -np 8 ../../espresso-5.1.1/bin//epcdft_coupling.x < $^ > $@

#pp.out: pp.in
#	export OMP_NUM_THREADS=1 ; mpirun -np 8 ../../espresso-5.1.1/bin/pp.x < $^ > $@
