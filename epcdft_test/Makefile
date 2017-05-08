# EXES AND SCRIPTS
BINDIR=../bin/
PW_EXE=${BINDIR}/pw.x
PP_EXE=${BINDIR}/pp.x
CUP_EXE=${BINDIR}epcdft_coupling.x
NP:=$(shell grep -c ^processor /proc/cpuinfo)
PW_CMD=export OMP_NUM_THREADS=1; export MKL_NUM_THREADS=1; mpirun -np $(NP) $(PW_EXE)
PP_CMD=export OMP_NUM_THREADS=1; export MKL_NUM_THREADS=1; mpirun -np $(NP) $(PP_EXE)
CUP_CMD=export OMP_NUM_THREADS=1; export MKL_NUM_THREADS=1; mpirun -np $(NP) $(CUP_EXE)
CUP_SH=./setup_coupling_input.sh 
SDD_SH=./make_sdd.sh

all: left.out right.out coupling.in coupling.out left_sdd.in left_sdd.out right_sdd.in right_sdd.out
clean: 
	rm -rf out *.out coupling.in *t_sdd.* *.cube txt *.sdd

# the meat
left.out: left.in
	$(PW_CMD) -in $^ > $@

right.out: right.in
	$(PW_CMD) -in $^ > $@

coupling.in: left.out right.out
	bash $(CUP_SH) $^ > $@

coupling.out: coupling.in
	$(CUP_CMD) -in $^ > $@

left_sdd.in: left.out
	$(SDD_SH) $^

left_sdd.out: left_sdd.in
	$(PP_CMD) -in $^ >$@

right_sdd.in: right.out
	$(SDD_SH) $^

right_sdd.out: right_sdd.in
	$(PP_CMD) -in $^ >$@
