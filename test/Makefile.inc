PW_EXE=../../espresso-5.1.1/bin/pw.x
PP_EXE=../../espresso-5.1.1/bin/pp.x
NP=8
PW_CMD=export OMP_NUM_THREADS=1; export MKL_NUM_THREADS=1; mpirun -np $(NP) $(PW_EXE)
PP_CMD=export OMP_NUM_THREADS=1; export MKL_NUM_THREADS=1; mpirun -np $(NP) $(PP_EXE)

