all: pw.out
clean: 
	rm -rf out pw.out

pw.out: pw.in
	export OMP_NUM_THREADS=1 ; mpirun -np 8 ../../bin/pw.x < $^ > $@

