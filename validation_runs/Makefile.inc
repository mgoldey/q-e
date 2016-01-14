PW_EXE=../../espresso-5.1.1/bin/pw.x
PP_EXE=../../espresso-5.1.1/bin/pp.x
PP_EXE=../../espresso-5.1.1/bin/epcdft_coupling.x
NP:=$(shell grep -c ^processor /proc/cpuinfo)
PW_CMD=export OMP_NUM_THREADS=1; export MKL_NUM_THREADS=1; nohup mpirun -np $(NP) $(PW_EXE)
PP_CMD=export OMP_NUM_THREADS=1; export MKL_NUM_THREADS=1; nohup mpirun -np $(NP) $(PP_EXE)
CUP_CMD=export OMP_NUM_THREADS=1; export MKL_NUM_THREADS=1; nohup mpirun -np $(NP) $(PP_EXE)

