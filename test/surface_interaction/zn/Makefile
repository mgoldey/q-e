#include ../../Makefile.inc
PW_CMD = mpirun -np 4 ../../../espresso-5.1.1/bin/pw.x

all: pw2.out pw2.25.out pw2.50.out pw2.75.out pw3.out pw4.out pw5.out pw6.out

clean: 
	rm -rf *.out out *cub

pw2.out: pw2.in
	$(PW_CMD) < $^ > $@

pw2.25.out: pw2.25.in
	$(PW_CMD) < $^ > $@

pw2.50.out: pw2.50.in
	$(PW_CMD) < $^ > $@

pw2.75.out: pw2.75.in
	$(PW_CMD) < $^ > $@

pw3.out: pw3.in
	$(PW_CMD) < $^ > $@

pw4.out: pw4.in
	$(PW_CMD) < $^ > $@

pw5.out: pw5.in
	$(PW_CMD) < $^ > $@

pw6.out: pw6.in
	$(PW_CMD) < $^ > $@
