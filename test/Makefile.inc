QE_EXE=../../espresso-5.1.1/bin/pw.x
NP=8
QE_CMD=export OMP_NUM_THREADS=1; export MKL_NUM_THREADS=1; mpirun -np $(NP) $(QE_EXE)

