include ../Makefile.inc
all: left.out right.out coupling.in coupling.out
clean: 
	rm -rf *.out out coupling.in 

left.out: left.in
	$(PW_CMD) < $^ > $@

right.out: right.in
	$(PW_CMD) < $^ > $@

coupling.in: left.out right.out
	sh ../../scripts/setup_coupling_input.sh left.out right.out > coupling.in

coupling.out: coupling.in
	$(CUP_CMD) < $^ > $@
