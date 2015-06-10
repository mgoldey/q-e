all: pw.out pw2.out pw_serial.out
clean: 
	rm -rf out pw.out pw2.out pw_serial.out

pw.out: pw.in
	export OMP_NUM_THREADS=1 ; mpirun -np 8 ../../espresso-5.1.1/bin/pw.x < $^ > $@
pw_serial.out: pw.in
	export OMP_NUM_THREADS=8 ; mpirun -np 1 ../../espresso-5.1.1/bin/pw.x < $^ > $@
pw2.out: pw2.in
	export OMP_NUM_THREADS=1 ; mpirun -np 8 ../../espresso-5.1.1/bin/pw.x < $^ > $@

