include ../Makefile.inc
all: left.out right.out coupling.in coupling.out left_pp.out right_pp.out
clean: 
	rm -rf *.out out 

left.out: left.in
	$(PW_CMD) < $^ > $@

right.out: right.in
	$(PW_CMD) < $^ > $@

coupling.in: left.out right.out
	sh ../../scripts/setup_coupling_input.sh left.out right.out > coupling.in

coupling.out: coupling.in
	$(CUP_CMD) < $^ > $@

left_pp.out: left_pp.in
	$(PP_CMD)  < $^ > $@

right_pp.out: right_pp.in
	$(PP_CMD)  < $^ > $@
