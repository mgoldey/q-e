#include ../../Makefile.inc
PW_CMD = mpirun -np 4 ../../../espresso-5.1.1/bin/pw.x

all: pw.out 

clean: 
	rm -rf *.out out 

pw.out: pw.in
	$(PW_CMD) < $^ > $@
