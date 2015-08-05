include ../Makefile.inc
all: hirshfeld.out voronoi.out zero.out hirsh_homo.out voro_homo.out homo.out hirsh_poten.out voro_poten.out
clean: 
	rm -rf out *.out tmp.pp *.cube

# the meat
hirshfeld.out: hirshfeld.in 
	$(PW_CMD) < $^ > $@
voronoi.out: voronoi.in
	$(PW_CMD) < $^ > $@
zero.out: zero.in
	$(PW_CMD) < $^ > $@

# presently broken - fixable in 20 minutes
well.out: well.in
	$(PW_CMD) < $^ > $@


# plots
hirsh_homo.out: hirsh_homo.in
	$(PP_CMD) < $^ > $@
voro_homo.out: voro_homo.in
	$(PP_CMD) < $^ > $@
homo.out: homo.in
	$(PP_CMD) < $^ > $@
hirsh_poten.out: hirsh_poten.in
	$(PP_CMD) < $^ > $@
voro_poten.out: voro_poten.in
	$(PP_CMD) < $^ > $@

#coupling.out: coupling.in
#	export OMP_NUM_THREADS=1 ; mpirun -np 8 ../../espresso-5.1.1/bin//epcdft_coupling.x < $^ > $@
