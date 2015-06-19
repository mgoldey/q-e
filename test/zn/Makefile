all: pw.out pw_rerun.out pw2.out coupling.out 
clean: 
	rm -rf out out2 pw.out pw2.out coupling.out

pw.out: pw.in
	export OMP_NUM_THREADS=1 ; mpirun -np 8 ../../espresso-5.1.1/bin/pw.x < $^ > $@

pw_rerun.out: pw_rerun.in
	export OMP_NUM_THREADS=1 ; mpirun -np 8 ../../espresso-5.1.1/bin/pw.x < $^ > $@

pw2.out: pw2.in
	export OMP_NUM_THREADS=1 ; mpirun -np 8 ../../espresso-5.1.1/bin/pw.x < $^ > $@

coupling.out: coupling.in
	export OMP_NUM_THREADS=1 ; mpirun -np 8 ../../espresso-5.1.1/bin//epcdft_coupling.x < $^ > $@
