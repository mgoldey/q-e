include ../Makefile.inc
all: left.out right.out coupling.in coupling.out
clean: 
	rm -rf *.out out 

left.out: left.in
	$(PW_CMD) < $^ > $@

right.out: right.in
	$(PW_CMD) < $^ > $@

coupling.in: ../../scripts/setup_coupling_input.sh
	sh ../../scripts/setup_coupling_input.sh left.out right.out > coupling.in

coupling.out: coupling.in
	$(CUP_CMD) < $^ > $@
